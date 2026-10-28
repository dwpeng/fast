//! Fast "first N records" extraction.
//!
//! `head` never re-serializes records: it locates the end of the N-th record
//! and copies the raw prefix.  Record/newline counts are computed in parallel
//! for large blocks, and the shared block reader keeps the I/O pipeline busy.

use std::io::Write;

use anyhow::Result;
use memchr::memchr;
use rayon::prelude::*;

use crate::parser::{Format, StreamConfig, fastq_record_end_at, stream_path};
use crate::util::is_broken_pipe;

struct HeadCut {
    /// Bytes to copy from the current block.
    cut: usize,
    /// Number of complete records in `[0, cut)`.
    records: u64,
    /// Whether the requested number of records has been satisfied.
    reached: bool,
}

/// Copy up to `*remaining` records from one input, decrementing the counter.
fn head_path(
    path: &str,
    remaining: &mut u64,
    cfg: &StreamConfig,
    out: &mut (dyn Write + Send),
) -> Result<()> {
    if *remaining == 0 {
        return Ok(());
    }

    let mut wrote_any = false;
    let mut last_byte = None;

    stream_path(path, cfg, |body, format, eof| {
        if *remaining == 0 {
            return Ok(false);
        }

        let cut = match format {
            Format::Fasta => find_fasta_cut(body, *remaining),
            Format::Fastq => find_fastq_cut(body, *remaining, eof),
        };

        if cut.cut > 0 {
            out.write_all(&body[..cut.cut])?;
            wrote_any = true;
            last_byte = Some(body[cut.cut - 1]);
        }

        *remaining = (*remaining).saturating_sub(cut.records);

        Ok(!cut.reached && !eof)
    })?;

    // Keep concatenated output valid when the source's final line has no
    // trailing newline.
    if wrote_any && last_byte != Some(b'\n') {
        out.write_all(b"\n")?;
    }

    Ok(())
}

/// Run `head` over every input.  `num` is the total number of records to print
/// across all inputs, matching `seqkit head`.
pub fn run(files: &[String], num: u64, cfg: &StreamConfig, out: &mut (dyn Write + Send)) -> Result<bool> {
    let mut remaining = num;
    let mut failed = false;

    for path in files {
        if remaining == 0 {
            break;
        }

        if let Err(e) = head_path(path, &mut remaining, cfg, out) {
            if is_broken_pipe(&e) {
                return Err(e);
            }
            failed = true;
            eprintln!("fast head: {path}: {e:#}");
        }
    }

    Ok(failed)
}

fn find_fasta_cut(data: &[u8], remaining: u64) -> HeadCut {
    if remaining == 0 {
        return HeadCut {
            cut: 0,
            records: 0,
            reached: true,
        };
    }

    if data.is_empty() {
        return HeadCut {
            cut: 0,
            records: 0,
            reached: false,
        };
    }

    // The first record starts at byte 0 and is not counted by
    // `scan_fasta`.
    let scan = scan_fasta(data, remaining);

    match scan.pos {
        Some(cut) => HeadCut {
            cut,
            records: remaining,
            reached: true,
        },
        None => {
            let total = scan.total + 1;
            HeadCut {
                cut: data.len(),
                records: total,
                reached: total == remaining,
            }
        }
    }
}

struct Scan {
    /// Position of the target-th record start/newline, when found.
    pos: Option<usize>,
    /// Total count in the block.  When `pos` is `Some` the scan may have
    /// stopped early, so this is only a lower bound; callers only use it when
    /// `pos` is `None`.
    total: u64,
}

/// Count `>` bytes that start a line (and therefore a FASTA record), excluding
/// the initial record at byte 0.
fn count_fasta_extra_starts(data: &[u8], start: usize, end: usize) -> u64 {
    let mut count = 0;
    let mut i = start.max(1);

    while i < end {
        match memchr(b'>', &data[i..end]) {
            Some(rel) => {
                let pos = i + rel;
                if data[pos - 1] == b'\n' {
                    count += 1;
                }
                i = pos + 1;
            }
            None => break,
        }
    }

    count
}

fn find_fasta_extra_start(data: &[u8], start: usize, end: usize, n: u64) -> Option<usize> {
    let mut count = 0;
    let mut i = start.max(1);

    while i < end {
        match memchr(b'>', &data[i..end]) {
            Some(rel) => {
                let pos = i + rel;
                if data[pos - 1] == b'\n' {
                    count += 1;
                    if count == n {
                        return Some(pos);
                    }
                }
                i = pos + 1;
            }
            None => break,
        }
    }

    None
}

fn scan_fasta(data: &[u8], target: u64) -> Scan {
    let ranges = parallel_ranges(data.len());

    if ranges.len() <= 1 || target <= 1024 || data.len() < (2 << 20) {
        let mut count = 0u64;
        let mut i = 1;

        while i < data.len() {
            match memchr(b'>', &data[i..]) {
                Some(rel) => {
                    let pos = i + rel;
                    if data[pos - 1] == b'\n' {
                        count += 1;
                        if count == target {
                            return Scan {
                                pos: Some(pos),
                                total: count,
                            };
                        }
                    }
                    i = pos + 1;
                }
                None => break,
            }
        }

        return Scan {
            pos: None,
            total: count,
        };
    }

    let counts: Vec<u64> = ranges
        .par_iter()
        .map(|&(start, end)| count_fasta_extra_starts(data, start, end))
        .collect();

    let total: u64 = counts.iter().sum();
    if total < target {
        return Scan { pos: None, total };
    }

    let mut before = 0u64;
    for (i, &count) in counts.iter().enumerate() {
        if before + count >= target {
            let (start, end) = ranges[i];
            let need = target - before;
            let pos = find_fasta_extra_start(data, start, end, need)
                .or_else(|| find_fasta_extra_start(data, start, data.len(), need));
            return Scan { pos, total };
        }
        before += count;
    }

    Scan { pos: None, total }
}

fn find_fastq_cut(data: &[u8], remaining: u64, eof: bool) -> HeadCut {
    if remaining == 0 {
        return HeadCut {
            cut: 0,
            records: 0,
            reached: true,
        };
    }

    // Walk complete records from the start of the block.  The block is
    // record-aligned, and wrapped sequence/quality lines mean a fixed
    // four-lines-per-record count is not valid.
    let mut pos = 0usize;
    let mut records = 0u64;

    while pos < data.len() && records < remaining {
        match fastq_record_end_at(data, pos) {
            Some(end) if end > pos => {
                pos = end;
                records += 1;
            }
            _ => break,
        }
    }

    if records == remaining {
        return HeadCut {
            cut: pos,
            records,
            reached: true,
        };
    }

    // The block contains fewer than `remaining` records.  At EOF the final
    // record may be missing its trailing newline, and the whole remainder
    // belongs to the output; otherwise every byte in the block is part of a
    // complete record.
    if eof {
        let extra = u64::from(pos < data.len());
        HeadCut {
            cut: data.len(),
            records: records + extra,
            reached: true,
        }
    } else {
        HeadCut {
            cut: data.len(),
            records,
            reached: false,
        }
    }
}

fn parallel_ranges(len: usize) -> Vec<(usize, usize)> {
    if len == 0 {
        return Vec::new();
    }

    let threads = rayon::current_num_threads().max(1);
    let chunk_size = (len / (threads * 4).max(1)).max(1 << 20);

    if chunk_size >= len {
        return vec![(0, len)];
    }

    (0..len)
        .step_by(chunk_size)
        .map(|start| (start, (start + chunk_size).min(len)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fasta(records: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for i in 0..records {
            out.extend_from_slice(format!(">r{i}\nACGT\n").as_bytes());
        }
        out
    }

    fn fastq(records: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for i in 0..records {
            out.extend_from_slice(format!("@r{i}\nACGT\n+\nIIII\n").as_bytes());
        }
        out
    }

    fn wrapped_fastq(records: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for i in 0..records {
            // The first record is four lines; the rest wrap at 60 columns.
            let len = if i == 0 { 4 } else { 100 };
            let seq: String = std::iter::repeat_n('A', len).collect();
            let qual: String = std::iter::repeat_n('I', len).collect();

            out.extend_from_slice(format!("@r{i}\n").as_bytes());
            for line in seq.as_bytes().chunks(60) {
                out.extend_from_slice(line);
                out.push(b'\n');
            }
            out.extend_from_slice(b"+\n");
            for line in qual.as_bytes().chunks(60) {
                out.extend_from_slice(line);
                out.push(b'\n');
            }
        }
        out
    }

    #[test]
    fn fasta_cut_selects_prefix() {
        let data = fasta(10);
        let cut = find_fasta_cut(&data, 3);
        assert!(cut.reached);
        assert_eq!(cut.records, 3);
        assert_eq!(&data[..cut.cut], &fasta(3)[..]);

        let cut = find_fasta_cut(&data, 10);
        assert!(cut.reached);
        assert_eq!(cut.cut, data.len());

        let cut = find_fasta_cut(&data, 20);
        assert!(!cut.reached);
        assert_eq!(cut.records, 10);
        assert_eq!(cut.cut, data.len());
    }

    #[test]
    fn fastq_cut_selects_prefix() {
        let data = fastq(10);
        let cut = find_fastq_cut(&data, 3, false);
        assert!(cut.reached);
        assert_eq!(cut.records, 3);
        assert_eq!(&data[..cut.cut], &fastq(3)[..]);

        let cut = find_fastq_cut(&data, 10, false);
        assert!(cut.reached);
        assert_eq!(cut.cut, data.len());

        let cut = find_fastq_cut(&data, 20, false);
        assert!(!cut.reached);
        assert_eq!(cut.cut, data.len());
    }

    #[test]
    fn fastq_without_final_newline_at_eof() {
        let mut data = fastq(1);
        data.pop();

        let cut = find_fastq_cut(&data, 1, true);
        assert!(cut.reached);
        assert_eq!(cut.cut, data.len());
    }

    #[test]
    fn wrapped_fastq_cut_selects_prefix() {
        let data = wrapped_fastq(10);
        let cut = find_fastq_cut(&data, 3, false);
        assert!(cut.reached);
        assert_eq!(cut.records, 3);
        assert_eq!(&data[..cut.cut], &wrapped_fastq(3)[..]);

        let cut = find_fastq_cut(&data, 10, false);
        assert!(cut.reached);
        assert_eq!(cut.cut, data.len());

        let cut = find_fastq_cut(&data, 20, false);
        assert!(!cut.reached);
        assert_eq!(cut.records, 10);
        assert_eq!(cut.cut, data.len());
    }

    #[test]
    fn fasta_cut_large_parallel_path() {
        let data = fasta(200_000);
        assert!(data.len() > (2 << 20));

        let cut = find_fasta_cut(&data, 150_000);
        assert!(cut.reached);
        assert_eq!(&data[..cut.cut], &fasta(150_000)[..]);
    }

    #[test]
    fn fastq_cut_large_parallel_path() {
        let data = fastq(200_000);
        assert!(data.len() > (2 << 20));

        let cut = find_fastq_cut(&data, 150_000, false);
        assert!(cut.reached);
        assert_eq!(&data[..cut.cut], &fastq(150_000)[..]);
    }
}
