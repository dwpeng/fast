//! FASTA <-> FASTQ conversion.
//!
//! Conversion is done on raw bytes rather than through a record model: line
//! ends are located with SIMD `memchr`, sequence lines are bulk-copied, and
//! each record-aligned block is split into chunks that are converted in
//! parallel.  The result is written in input order.

use std::io::Write;

use anyhow::{Context, Result, bail};
use memchr::{memchr, memmem};
use rayon::prelude::*;

use crate::parser::{Format, StreamConfig, plan_chunks, stream_path};

/// Convert FASTA to FASTQ, using `qual` (a Phred+33 quality byte) for every
/// generated quality score.
pub fn fa2fq(path: &str, qual: u8, cfg: &StreamConfig, out: &mut dyn Write) -> Result<()> {
    stream_path(path, cfg, |body, format, _eof| {
        if format != Format::Fasta {
            bail!("expected FASTA input, found {}", format.label());
        }

        let chunks: Vec<Result<Vec<u8>>> = plan_chunks(body, Format::Fasta)
            .par_iter()
            .map(|&(start, end)| {
                convert_fasta_chunk(&body[start..end], qual)
                    .with_context(|| format!("records in byte range {start}..{end}"))
            })
            .collect();

        for chunk in chunks {
            out.write_all(&chunk?)?;
        }

        Ok(true)
    })?;

    Ok(())
}

/// Convert FASTQ to FASTA, dropping quality scores.
pub fn fq2fa(path: &str, cfg: &StreamConfig, out: &mut dyn Write) -> Result<()> {
    stream_path(path, cfg, |body, format, _eof| {
        if format != Format::Fastq {
            bail!("expected FASTQ input, found {}", format.label());
        }

        let chunks: Vec<Result<Vec<u8>>> = plan_chunks(body, Format::Fastq)
            .par_iter()
            .map(|&(start, end)| {
                convert_fastq_chunk(&body[start..end])
                    .with_context(|| format!("records in byte range {start}..{end}"))
            })
            .collect();

        for chunk in chunks {
            out.write_all(&chunk?)?;
        }

        Ok(true)
    })?;

    Ok(())
}

/// Convert one record-aligned FASTA range.  `data` must start with `>` and
/// every `\n>` inside it starts a new record.
fn convert_fasta_chunk(data: &[u8], qual: u8) -> Result<Vec<u8>> {
    // FASTQ output is roughly sequence + quality plus a few extra newlines.
    let mut out = Vec::with_capacity(data.len().saturating_mul(2).saturating_add(64));
    let mut pos = 0;

    while pos < data.len() {
        if data[pos] != b'>' {
            bail!("expected '>' at offset {pos}");
        }

        let header_start = pos + 1;
        let header_end = memchr(b'\n', &data[header_start..])
            .map(|i| header_start + i)
            .unwrap_or(data.len());
        let header = strip_cr(&data[header_start..header_end]);

        let seq_start = if header_end < data.len() {
            header_end + 1
        } else {
            data.len()
        };
        // An empty sequence is followed immediately by the next header, so
        // the `>` is at `seq_start` rather than after a `\n` found by the
        // search below.
        let next = if data.get(seq_start) == Some(&b'>') {
            seq_start
        } else {
            memmem::find(&data[seq_start..], b"\n>")
                .map(|i| seq_start + i + 1)
                .unwrap_or(data.len())
        };

        out.push(b'@');
        out.extend_from_slice(header);
        out.push(b'\n');

        let seq_out_start = out.len();
        copy_fasta_sequence(&data[seq_start..next], &mut out);
        let seq_len = out.len() - seq_out_start;

        out.push(b'\n');
        out.push(b'+');
        out.push(b'\n');
        out.resize(out.len() + seq_len, qual);
        out.push(b'\n');

        pos = next;
    }

    Ok(out)
}

/// Copy FASTA sequence lines without line terminators.
#[inline]
fn copy_fasta_sequence(mut data: &[u8], out: &mut Vec<u8>) {
    while !data.is_empty() {
        match memchr(b'\n', data) {
            Some(i) => {
                out.extend_from_slice(strip_cr(&data[..i]));
                data = &data[i + 1..];
            }
            None => {
                out.extend_from_slice(strip_cr(data));
                break;
            }
        }
    }
}

/// Convert one record-aligned four-line FASTQ range.
fn convert_fastq_chunk(data: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(data.len() / 2 + 64);
    let mut pos = 0;

    while pos < data.len() {
        if data[pos] != b'@' {
            bail!("expected '@' at offset {pos}");
        }

        let header_start = pos + 1;
        let header_end = memchr(b'\n', &data[header_start..])
            .map(|i| header_start + i)
            .with_context(|| format!("unterminated FASTQ header at offset {pos}"))?;
        let header = strip_cr(&data[header_start..header_end]);

        let seq_start = header_end + 1;
        let seq_end = memchr(b'\n', &data[seq_start..])
            .map(|i| seq_start + i)
            .with_context(|| format!("unterminated FASTQ sequence at offset {pos}"))?;
        let seq = strip_cr(&data[seq_start..seq_end]);

        let plus_start = seq_end + 1;
        if data.get(plus_start) != Some(&b'+') {
            bail!("expected '+' at offset {plus_start}");
        }
        let plus_end = memchr(b'\n', &data[plus_start..])
            .map(|i| plus_start + i)
            .with_context(|| format!("unterminated FASTQ separator at offset {pos}"))?;

        let qual_start = plus_end + 1;
        let qual_end = memchr(b'\n', &data[qual_start..])
            .map(|i| qual_start + i)
            .unwrap_or(data.len());
        let qual = strip_cr(&data[qual_start..qual_end]);

        if seq.len() != qual.len() {
            bail!(
                "FASTQ record at offset {pos} has sequence length {} but quality length {}",
                seq.len(),
                qual.len()
            );
        }

        out.push(b'>');
        out.extend_from_slice(header);
        out.push(b'\n');
        out.extend_from_slice(seq);
        out.push(b'\n');

        pos = if qual_end < data.len() {
            qual_end + 1
        } else {
            data.len()
        };
    }

    Ok(out)
}

#[inline]
fn strip_cr(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\r").unwrap_or(line)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_to_fastq_multiline_and_crlf() {
        let data = b">a description\r\nACGT\r\nNN\r\n>b\nacgt\n";
        let out = convert_fasta_chunk(data, b'I').unwrap();
        assert_eq!(
            out,
            b"@a description\nACGTNN\n+\nIIIIII\n@b\nacgt\n+\nIIII\n"
        );
    }

    #[test]
    fn fastq_to_fasta_multiline() {
        let data = b"@r1 desc\nACGT\n+\nIIII\n@r2\nNN\n+r2\n##\n";
        let out = convert_fastq_chunk(data).unwrap();
        assert_eq!(out, b">r1 desc\nACGT\n>r2\nNN\n");
    }

    #[test]
    fn fastq_quality_length_mismatch_is_an_error() {
        let data = b"@r1\nACGT\n+\nIII\n";
        assert!(convert_fastq_chunk(data).is_err());
    }

    #[test]
    fn fasta_empty_sequence_before_next_header() {
        let data = b">empty\n>full\nACGT\n";
        let out = convert_fasta_chunk(data, b'I').unwrap();
        assert_eq!(out, b"@empty\n\n+\n\n@full\nACGT\n+\nIIII\n");
    }
}
