//! Streaming, block-based input and record-parallel parsing.
//!
//! Input is read into fixed-size blocks (16 MiB by default), parsed in
//! parallel with rayon, and then the next block is read.  Memory use is
//! independent of the file size: two block buffers, one incomplete record and
//! the decoder's own small buffers are kept.
//!
//! A dedicated reader thread fills the next buffer while the calling thread
//! parses the previous one, so disk latency and decompression are overlapped
//! with CPU work.  Block boundaries are chosen at record boundaries:
//!
//! * FASTA: everything up to the last `\n>` is complete; the final record is
//!   carried over to the next block.
//! * FASTQ: the last complete four-line record is found by scanning backwards
//!   and validating the record structure; the remainder is carried over.
//!
//! Inside a block, record-aligned chunking is used for parallelism.  noodles
//! handles record decoding and reuses its record buffers.

use std::fs::File;
use std::io::{self, Cursor, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::OnceLock;
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::thread;
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail};
use memchr::{memchr, memmem};
use noodles_fasta as fasta;
use noodles_fastq as fastq;
use rapidgzip_core::{Decoder, DecoderPool};
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};

use crate::stats::{ChunkStats, StatOptions, Summary, guess_alphabet};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

impl Format {
    pub fn label(self) -> &'static str {
        match self {
            Format::Fasta => "FASTA",
            Format::Fastq => "FASTQ",
        }
    }
}

/// Streaming input configuration shared by all subcommands.
#[derive(Clone, Debug)]
pub struct StreamConfig {
    /// Number of bytes read from the input per streaming iteration.
    pub block_size: usize,
    /// Upper bound for the read size when a single record does not fit in
    /// `block_size` and the reader grows its buffer.
    pub max_block_size: usize,
    /// Plain-file read concurrency per block.  `None` selects automatically;
    /// `Some(0)` forces sequential reads.
    pub read_threads: Option<usize>,
    /// Print read throughput progress to stderr.
    pub verbose: bool,
}

impl StreamConfig {
    pub fn new(block_size: usize, read_threads: Option<usize>, verbose: bool) -> Self {
        Self {
            block_size,
            max_block_size: block_size,
            read_threads,
            verbose,
        }
    }

    /// Use a smaller initial block size without lowering the growth cap.
    pub fn with_block_size(&self, block_size: usize) -> Self {
        Self {
            block_size,
            ..self.clone()
        }
    }
}

/// Configuration for the statistics command.
pub struct Config {
    pub stats: StatOptions,
    pub nx: Vec<f64>,
    pub stream: StreamConfig,
}

const DEFAULT_READ_THREADS: usize = 32;
const DEFAULT_NETWORK_READ_THREADS: usize = 32;
const MAX_READ_THREADS: usize = 32;

/// Process one file and produce display-ready statistics.
pub fn process_path(path: &str, cfg: &Config) -> Result<Summary> {
    let mut stats = ChunkStats::default();
    let format = stream_path(path, &cfg.stream, |body, fmt, _eof| {
        stats.merge(parse_body(body, fmt, cfg)?);
        Ok(true)
    })?;

    finish_summary((stats, format), cfg)
}

fn finish_summary((stats, format): (ChunkStats, Option<Format>), cfg: &Config) -> Result<Summary> {
    match format {
        Some(format) => Ok(Summary::from_chunks(
            stats,
            format.label(),
            &cfg.nx,
            100_000,
        )),
        // Empty or whitespace/comment-only input.
        None => Ok(Summary::empty("", &cfg.nx)),
    }
}

/// Stream a FASTA/FASTQ input and hand every record-aligned block to
/// `on_block`.  Returns the detected format (`None` for empty/comment-only
/// input).  Returning `false` from the callback stops reading early.
///
/// This is the shared input path for stats, conversions, head and merge.  It
/// uses the parallel `pread` reader for plain files and the parallel gzip
/// decoder for compressed input.
pub fn stream_path<G>(path: &str, cfg: &StreamConfig, on_block: G) -> Result<Option<Format>>
where
    G: FnMut(&[u8], Format, bool) -> Result<bool>,
{
    if path == "-" {
        return stream_stdin(cfg, on_block);
    }

    let mut file = File::open(path).with_context(|| format!("failed to open {path}"))?;
    let meta = file
        .metadata()
        .with_context(|| format!("failed to stat {path}"))?;

    if meta.len() == 0 {
        return Ok(None);
    }

    let mut magic = [0u8; 2];
    if meta.len() >= 2 {
        file.read_exact(&mut magic)
            .with_context(|| format!("failed to read {path}"))?;
        file.seek(SeekFrom::Start(0))
            .with_context(|| format!("failed to seek {path}"))?;
    }

    if meta.len() >= 2 && magic == [0x1f, 0x8b] {
        let decoder = gzip_decoder()?;
        let mut reader = decoder
            .open(path)
            .with_context(|| format!("failed to open {path} for decompression"))?;
        stream_blocks_with(
            cfg,
            path,
            move |buf| read_up_to(&mut reader, buf).context("failed to read input"),
            on_block,
        )
    } else {
        let (read_threads, prefetch) = read_strategy(path, cfg.read_threads);
        let mut filler = FileFiller {
            file: &file,
            offset: 0,
            len: meta.len(),
            read_threads,
            prefetch,
        };

        stream_blocks_with(cfg, path, move |buf| filler.fill(buf), on_block)
    }
}

fn stream_stdin<G>(cfg: &StreamConfig, on_block: G) -> Result<Option<Format>>
where
    G: FnMut(&[u8], Format, bool) -> Result<bool>,
{
    let stdin = io::stdin();

    // Peek at the gzip magic without consuming more than two bytes.  The bytes
    // are chained back in front of the remaining stream.
    let mut prefix = [0u8; 2];
    let n = {
        let mut lock = stdin.lock();
        read_up_to(&mut lock, &mut prefix).context("failed to read stdin")?
    };
    let mut chained = Cursor::new(prefix[..n].to_vec()).chain(stdin);

    if n == 2 && prefix == [0x1f, 0x8b] {
        let decoder = gzip_decoder()?;
        let mut reader = decoder
            .stream_reader(chained)
            .context("failed to initialize gzip decoding for stdin")?;
        stream_blocks_with(
            cfg,
            "-",
            move |buf| read_up_to(&mut reader, buf).context("failed to read input"),
            on_block,
        )
    } else {
        stream_blocks_with(
            cfg,
            "-",
            move |buf| read_up_to(&mut chained, buf).context("failed to read input"),
            on_block,
        )
    }
}

fn read_strategy(path: &str, read_threads: Option<usize>) -> (usize, bool) {
    match read_threads {
        Some(n) => (n, n > 0),
        None => match fs_kind(Path::new(path)) {
            // High-latency network storage: more outstanding requests help,
            // while fadvise is usually a no-op or counterproductive.
            FsKind::Network => (DEFAULT_NETWORK_READ_THREADS, false),
            // Large RAID arrays and SSDs both benefit from queue depth;
            // adjacent pieces keep the access pattern sequential.
            FsKind::Rotational | FsKind::NonRotational | FsKind::Unknown => {
                (DEFAULT_READ_THREADS, true)
            }
        },
    }
}

/// Read and parse blocks from a `Read` source until EOF.
#[cfg(test)]
fn stream_stats<R: Read + Send>(
    mut reader: R,
    cfg: &Config,
    label: &str,
) -> Result<(ChunkStats, Option<Format>)> {
    stream_stats_with(cfg, label, move |buf| {
        read_up_to(&mut reader, buf).context("failed to read input")
    })
}

/// A record-aligned block handed from the reader thread to the parser.
struct Block {
    buf: Vec<u8>,
    len: usize,
    format: Option<Format>,
    eof: bool,
}

/// Read and parse blocks until EOF.
///
/// A dedicated reader thread fills the next buffer while the calling thread
/// parses the previous one, so disk/decompression latency is overlapped with
/// CPU work.  Two buffers bound memory use to roughly `2 * block_size` (plus
/// one incomplete record), independent of the input size.
#[cfg(test)]
fn stream_stats_with<F>(cfg: &Config, label: &str, fill: F) -> Result<(ChunkStats, Option<Format>)>
where
    F: FnMut(&mut [u8]) -> Result<usize> + Send,
{
    let mut stats = ChunkStats::default();
    let format = stream_blocks_with(&cfg.stream, label, fill, |body, fmt, _eof| {
        stats.merge(parse_body(body, fmt, cfg)?);
        Ok(true)
    })?;

    Ok((stats, format))
}

/// Read record-aligned blocks and invoke `on_block` for each one.
///
/// A dedicated reader thread fills the next buffer while the calling thread
/// processes the previous one, so disk/decompression latency is overlapped with
/// CPU work.  Two buffers bound memory use to roughly `2 * block_size` (plus
/// one incomplete record), independent of the input size.
fn stream_blocks_with<F, G>(
    cfg: &StreamConfig,
    label: &str,
    mut fill: F,
    mut on_block: G,
) -> Result<Option<Format>>
where
    F: FnMut(&mut [u8]) -> Result<usize> + Send,
    G: FnMut(&[u8], Format, bool) -> Result<bool>,
{
    let block_size = cfg.block_size.max(1);
    let max_block_size = cfg.max_block_size.max(block_size);
    let (full_tx, full_rx) = sync_channel::<Block>(2);
    let (free_tx, free_rx) = sync_channel::<Vec<u8>>(2);

    for _ in 0..2 {
        let _ = free_tx.send(vec![0u8; block_size]);
    }

    thread::scope(|scope| -> Result<Option<Format>> {
        let reader = scope
            .spawn(move || reader_loop(&mut fill, block_size, max_block_size, full_tx, free_rx));

        let mut format = None;
        let started = Instant::now();
        let mut last_report = started;
        let mut bytes_read = 0u64;
        let mut stopped = false;

        while let Ok(block) = full_rx.recv() {
            if block.format.is_some() {
                format = block.format;
            }

            if block.len > 0 {
                let fmt = format.expect("a non-empty block has a detected format");
                if !on_block(&block.buf[..block.len], fmt, block.eof)? {
                    stopped = true;
                    bytes_read += block.len as u64;
                    break;
                }
            }

            bytes_read += block.len as u64;
            if cfg.verbose && !block.eof && last_report.elapsed() >= Duration::from_secs(5) {
                report_progress(label, bytes_read, started.elapsed());
                last_report = Instant::now();
            }

            let _ = free_tx.send(block.buf);

            if block.eof {
                break;
            }
        }

        if stopped {
            // Unblock the reader: it may be waiting to hand us a block (the
            // full channel is dropped) or to get a free buffer back.
            drop(full_rx);
            drop(free_tx);
        }

        match reader.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_) => return Err(anyhow::anyhow!("reader thread panicked")),
        }

        if cfg.verbose {
            report_progress(label, bytes_read, started.elapsed());
        }

        Ok(format)
    })
}

fn reader_loop<F>(
    fill: &mut F,
    initial_block_size: usize,
    max_block_size: usize,
    full_tx: SyncSender<Block>,
    free_rx: Receiver<Vec<u8>>,
) -> Result<()>
where
    F: FnMut(&mut [u8]) -> Result<usize>,
{
    let mut prefix: Vec<u8> = Vec::new();
    let mut format: Option<Format> = None;
    let mut buf = match free_rx.recv() {
        Ok(buf) => buf,
        Err(_) => return Ok(()),
    };
    let mut filled = 0usize;
    let mut block_size = initial_block_size;

    loop {
        // A block was handed to the parser; prepend the bytes belonging to the
        // incomplete trailing record to the next buffer.
        if !prefix.is_empty() {
            ensure_len(&mut buf, prefix.len() + block_size);
            buf[..prefix.len()].copy_from_slice(&prefix);
            filled = prefix.len();
            prefix.clear();
        }

        ensure_len(&mut buf, filled + block_size);

        let n = fill(&mut buf[filled..filled + block_size])?;
        let mut total = filled + n;
        let eof = n == 0;

        if total == 0 {
            let _ = full_tx.send(Block {
                buf,
                len: 0,
                format,
                eof: true,
            });
            return Ok(());
        }

        if format.is_none() {
            match detect_format(&buf[..total])? {
                (Some(f), offset) => {
                    format = Some(f);
                    if offset > 0 {
                        buf.copy_within(offset..total, 0);
                        total -= offset;
                    }
                }
                (None, _) => {
                    if eof {
                        let _ = full_tx.send(Block {
                            buf,
                            len: 0,
                            format: None,
                            eof: true,
                        });
                        return Ok(());
                    }
                    // Only whitespace/comments so far; keep reading.
                    filled = total;
                    continue;
                }
            }
        }

        let fmt = format.expect("format was just detected");

        if eof {
            let _ = full_tx.send(Block {
                buf,
                len: total,
                format,
                eof: true,
            });
            return Ok(());
        }

        let cut = match fmt {
            Format::Fasta => fasta_complete_prefix(&buf[..total]),
            Format::Fastq => fastq_complete_prefix(&buf[..total]),
        };

        if cut > 0 {
            prefix.clear();
            prefix.extend_from_slice(&buf[cut..total]);

            if full_tx
                .send(Block {
                    buf,
                    len: cut,
                    format,
                    eof: false,
                })
                .is_err()
            {
                return Ok(());
            }

            buf = match free_rx.recv() {
                Ok(buf) => buf,
                Err(_) => return Ok(()),
            };
            filled = 0;
            block_size = initial_block_size;
        } else {
            // One record (or a comment block) larger than the buffer: keep the
            // data in place and read more after it.  Grow the read size
            // geometrically up to the configured maximum so that
            // chromosome-sized records do not turn into thousands of tiny
            // reads (and, for FASTQ, thousands of full-buffer scans).
            filled = total;
            block_size = block_size.saturating_mul(2).min(max_block_size);
        }
    }
}

/// Grow a read buffer when a single record spans multiple fills.
///
/// Doubling avoids the quadratic copying that would result from growing by
/// the block size on every iteration for a chromosome-sized record.
fn ensure_len(buf: &mut Vec<u8>, needed: usize) {
    if buf.len() < needed {
        let new_len = needed.max(buf.len().saturating_mul(2)).max(1);
        buf.resize(new_len, 0);
    }
}

/// Fill `buf` unless EOF is reached.  A return value of 0 means EOF.
fn read_up_to<R: Read>(reader: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;

    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }

    Ok(filled)
}

/// Split a record-aligned block into chunks and merge the per-chunk results.
fn parse_body(body: &[u8], format: Format, cfg: &Config) -> Result<ChunkStats> {
    let chunks = plan_chunks(body, format);

    let chunk_stats: Vec<Result<ChunkStats>> = chunks
        .par_iter()
        .enumerate()
        .map(|(i, &(start, end))| {
            let slice = &body[start..end];
            let result = match format {
                Format::Fasta => parse_fasta(slice, &cfg.stats, i == 0),
                Format::Fastq => parse_fastq(slice, &cfg.stats, i == 0),
            };
            result.with_context(|| format!("records in byte range {start}..{end}"))
        })
        .collect();

    let mut merged = ChunkStats::default();
    for result in chunk_stats {
        merged.merge(result?);
    }

    Ok(merged)
}

/// Byte ranges covering `[0, len)`, sized for good load balancing.
fn chunk_ranges(len: usize, threads: usize) -> Vec<(usize, usize)> {
    if len == 0 {
        return Vec::new();
    }

    let target_chunks = threads.saturating_mul(8).max(1);
    let chunk_size = (len / target_chunks).clamp(1 << 20, 64 << 20);

    if chunk_size >= len {
        return vec![(0, len)];
    }

    (0..len)
        .step_by(chunk_size)
        .map(|start| (start, (start + chunk_size).min(len)))
        .collect()
}

pub(crate) fn plan_chunks(data: &[u8], format: Format) -> Vec<(usize, usize)> {
    let threads = rayon::current_num_threads();
    let ranges = chunk_ranges(data.len(), threads);

    let mut boundaries: Vec<usize> = match format {
        Format::Fasta => ranges
            .par_iter()
            .map(|&(start, _)| align_fasta_start(data, start))
            .collect(),
        Format::Fastq => ranges
            .par_iter()
            .map(|&(start, _)| align_fastq_start(data, start))
            .collect(),
    };

    boundaries.push(data.len());
    boundaries.sort_unstable();
    boundaries.dedup();

    boundaries
        .windows(2)
        .filter_map(|w| (w[0] < w[1]).then_some((w[0], w[1])))
        .collect()
}

/// Next `>` definition line start at or after `start` (byte offset).
#[inline]
fn align_fasta_start(data: &[u8], start: usize) -> usize {
    if start == 0 {
        return 0;
    }
    if data.get(start) == Some(&b'>') && data.get(start.wrapping_sub(1)) == Some(&b'\n') {
        return start;
    }
    match memmem::find(&data[start..], b"\n>") {
        Some(i) => start + i + 1,
        None => data.len(),
    }
}

/// Find the next complete FASTQ record start at or after `start`.
///
/// The block starts on a record boundary, so scanning forward from each chunk
/// boundary needs to inspect at most one record; no global line counting is
/// required.
#[inline]
fn align_fastq_start(data: &[u8], start: usize) -> usize {
    let mut pos = start;

    // Move to the next line start if `pos` is in the middle of a line.
    if pos > 0 && data[pos - 1] != b'\n' {
        match memchr(b'\n', &data[pos..]) {
            Some(i) => pos += i + 1,
            None => return data.len(),
        }
    }

    while pos < data.len() {
        if data[pos] == b'@' && is_fastq_record_start(data, pos) {
            return pos;
        }
        match memchr(b'\n', &data[pos..]) {
            Some(i) => pos += i + 1,
            None => return data.len(),
        }
    }

    data.len()
}

/// Byte offset of the last complete FASTA record in the block.
///
/// The last record may continue in the next block, so it is always carried
/// over (0 is returned when the block contains only one record).
fn fasta_complete_prefix(data: &[u8]) -> usize {
    memmem::rfind(data, b"\n>").map_or(0, |i| i + 1)
}

/// Byte offset of the last complete FASTQ record in the block.
///
/// Instead of counting every line, scan backwards for a candidate `@` line and
/// validate the four-line record structure.  Normally this only examines the
/// final record (or the previous one when the block ends mid-record).
fn fastq_complete_prefix(data: &[u8]) -> usize {
    let mut end = data.len();

    while let Some(pos) = memmem::rfind(&data[..end], b"\n@") {
        let start = pos + 1;
        if is_fastq_record_start(data, start) {
            return start;
        }
        end = pos;
    }

    0
}

/// Check whether a four-line FASTQ record starting at `start` is complete.
fn is_fastq_record_start(data: &[u8], start: usize) -> bool {
    if data.get(start) != Some(&b'@') {
        return false;
    }

    // Definition line.
    let Some(i) = memchr(b'\n', &data[start..]) else {
        return false;
    };
    let seq_start = start + i + 1;

    // Sequence line.
    let Some(i) = memchr(b'\n', &data[seq_start..]) else {
        return false;
    };
    let plus_start = seq_start + i + 1;
    if data.get(plus_start) != Some(&b'+') {
        return false;
    }
    let seq_len = plus_start - 1 - seq_start;

    // Separator line.
    let Some(i) = memchr(b'\n', &data[plus_start..]) else {
        return false;
    };
    let qual_start = plus_start + i + 1;

    // Quality line (may be unterminated at EOF).
    let qual_end = match memchr(b'\n', &data[qual_start..]) {
        Some(i) => qual_start + i,
        None => data.len(),
    };
    let qual_len = qual_end - qual_start;

    seq_len == qual_len
}

fn parse_fasta(slice: &[u8], opts: &StatOptions, capture_alphabet: bool) -> Result<ChunkStats> {
    let mut reader = fasta::io::Reader::new(slice);
    let mut definition = fasta::record::Definition::default();
    let mut sequence = Vec::new();
    let mut stats = ChunkStats::default();

    loop {
        let n = reader
            .read_definition(&mut definition)
            .context("failed to read FASTA definition")?;
        if n == 0 {
            break;
        }

        sequence.clear();
        reader
            .read_sequence(&mut sequence)
            .context("failed to read FASTA sequence")?;

        if capture_alphabet && stats.alphabet.is_none() {
            stats.alphabet = Some(guess_alphabet(&sequence));
        }

        stats
            .add(&sequence, None, opts)
            .context("failed to accumulate FASTA record")?;
    }

    Ok(stats)
}

fn parse_fastq(slice: &[u8], opts: &StatOptions, capture_alphabet: bool) -> Result<ChunkStats> {
    let mut reader = fastq::io::Reader::new(slice);
    let mut record = fastq::Record::default();
    let mut stats = ChunkStats::default();

    loop {
        let n = reader
            .read_record(&mut record)
            .context("failed to read FASTQ record")?;
        if n == 0 {
            break;
        }

        if capture_alphabet && stats.alphabet.is_none() {
            stats.alphabet = Some(guess_alphabet(record.sequence()));
        }

        stats
            .add(record.sequence(), Some(record.quality_scores()), opts)
            .with_context(|| {
                format!(
                    "failed to accumulate FASTQ record {:?}",
                    String::from_utf8_lossy(record.name())
                )
            })?;
    }

    Ok(stats)
}

/// Detect the format and return the offset of the first record marker.
///
/// Leading blank lines and comment lines (`#...` at the start of a line) are
/// skipped, matching seqkit's tolerance for leading comments.
fn detect_format(data: &[u8]) -> Result<(Option<Format>, usize)> {
    let mut i = if data.starts_with(&[0xef, 0xbb, 0xbf]) {
        3
    } else {
        0
    };

    while i < data.len() {
        match data[i] {
            b'\n' | b'\r' | b'\t' | b' ' => i += 1,
            b'>' => return Ok((Some(Format::Fasta), i)),
            b'@' => return Ok((Some(Format::Fastq), i)),
            b'#' if i == 0 || data[i - 1] == b'\n' => {
                while i < data.len() && data[i] != b'\n' {
                    i += 1;
                }
            }
            b => bail!(
                "not a FASTA/FASTQ file (unexpected byte {:?} at offset {i})",
                b as char
            ),
        }
    }

    Ok((None, i))
}

/// Process-wide decoder pool: decompression threads are shared by all files
/// that are being processed concurrently, avoiding thread oversubscription when
/// many gzip files are given at once.
static DECODER_POOL: OnceLock<DecoderPool> = OnceLock::new();

fn gzip_decoder() -> Result<Decoder> {
    let threads = rayon::current_num_threads().max(1);

    let pool = DECODER_POOL.get_or_init(|| {
        DecoderPool::builder()
            .workers(threads)
            .build()
            .expect("decoder pool workers must be non-zero")
    });

    Decoder::builder()
        .decoder_threads(threads)
        .decoder_pool(pool.clone())
        // Keep the decoder's internal buffers small: this is a streaming
        // reader, not a decode-to-memory path.
        .decoded_chunk_size(1 << 20)
        .input_page_size(1 << 18)
        .in_flight_chunks(threads.clamp(2, 8))
        .build()
        .context("failed to configure the parallel gzip decoder")
}

/// Dedicated pool for `pread`-based block reads.  Keeping it separate from the
/// parsing pool avoids a deadlock: the parsing thread may be blocked waiting
/// for the reader, so reader tasks must not depend on the parsing workers.
static READ_POOL: OnceLock<ThreadPool> = OnceLock::new();

fn read_pool() -> &'static ThreadPool {
    READ_POOL.get_or_init(|| {
        let threads = rayon::current_num_threads().clamp(1, MAX_READ_THREADS);
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .thread_name(|i| format!("fast-read-{i}"))
            .build()
            .expect("failed to build the read thread pool")
    })
}

/// Fills blocks from a plain file.
struct FileFiller<'a> {
    file: &'a File,
    offset: u64,
    len: u64,
    read_threads: usize,
    prefetch: bool,
}

impl FileFiller<'_> {
    fn fill(&mut self, buf: &mut [u8]) -> Result<usize> {
        let n = read_block_parallel(
            self.file,
            &mut self.offset,
            self.len,
            buf,
            self.read_threads,
        )?;

        // Ask the kernel to prefetch the next block while the parser works on
        // the current one.  This keeps I/O in flight without another buffer.
        #[cfg(target_os = "linux")]
        if self.prefetch && self.offset < self.len {
            use std::os::unix::io::AsRawFd;

            let bytes = (self.len - self.offset).min(buf.len() as u64) as libc::off_t;
            unsafe {
                libc::posix_fadvise(
                    self.file.as_raw_fd(),
                    self.offset as libc::off_t,
                    bytes,
                    libc::POSIX_FADV_WILLNEED,
                );
            }
        }

        Ok(n)
    }
}

/// Filesystem classes relevant to the read strategy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FsKind {
    Rotational,
    NonRotational,
    Network,
    Unknown,
}

#[cfg(target_os = "linux")]
fn fs_kind(path: &Path) -> FsKind {
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;

    const NFS: i64 = 0x0000_6969;
    const SMB: i64 = 0x0000_517b;
    const SMB2: i64 = 0xfe53_4d42;
    const CIFS: i64 = 0xff53_4d42;
    const FUSE: i64 = 0x6573_5546;
    const AFS: i64 = 0x5346_414f;
    const LUSTRE: i64 = 0x0bd0_0bd0;
    const CEPH: i64 = 0x00c3_6400;
    const GPFS: i64 = 0x4750_4653;

    let Ok(c_path) = CString::new(path.as_os_str().as_bytes()) else {
        return FsKind::Unknown;
    };

    let mut stat: libc::statfs = unsafe { std::mem::zeroed() };
    if unsafe { libc::statfs(c_path.as_ptr(), &mut stat) } != 0 {
        return FsKind::Unknown;
    }

    let f_type = stat.f_type as i64;
    if matches!(
        f_type,
        NFS | SMB | SMB2 | CIFS | FUSE | AFS | LUSTRE | CEPH | GPFS
    ) {
        return FsKind::Network;
    }

    match is_rotational(path) {
        Some(true) => FsKind::Rotational,
        Some(false) => FsKind::NonRotational,
        None => FsKind::Unknown,
    }
}

#[cfg(not(target_os = "linux"))]
fn fs_kind(path: &Path) -> FsKind {
    match is_rotational(path) {
        Some(true) => FsKind::Rotational,
        Some(false) => FsKind::NonRotational,
        None => FsKind::Unknown,
    }
}

fn report_progress(label: &str, bytes: u64, elapsed: Duration) {
    let seconds = elapsed.as_secs_f64().max(1e-9);
    let gib = bytes as f64 / (1024.0 * 1024.0 * 1024.0);
    let mib_per_second = bytes as f64 / (1024.0 * 1024.0) / seconds;
    eprintln!("{label}: {gib:.2} GiB read, {mib_per_second:.0} MiB/s");
}

/// Report whether the filesystem holding `path` lives on a rotational device.
///
/// Returns `None` when this cannot be determined (e.g. network filesystems),
/// in which case the caller keeps the parallel read strategy.
#[cfg(target_os = "linux")]
fn is_rotational(path: &Path) -> Option<bool> {
    use std::os::unix::fs::MetadataExt;

    let metadata = std::fs::metadata(path).ok()?;
    let device = metadata.dev();

    let major = libc::major(device);
    let minor = libc::minor(device);
    let link = std::fs::read_link(format!("/sys/dev/block/{major}:{minor}")).ok()?;

    let link = link.to_string_lossy();
    let name = link.rsplit("block/").next()?.split('/').next()?.to_string();

    let rotational = std::fs::read_to_string(format!("/sys/block/{name}/queue/rotational")).ok()?;

    Some(rotational.trim() == "1")
}

#[cfg(not(target_os = "linux"))]
fn is_rotational(_path: &Path) -> Option<bool> {
    None
}

/// Read one block using positional reads.
///
/// With `max_pieces > 1` the block is split into that many pieces and filled
/// with concurrent `pread` calls (best for SSDs and parallel filesystems).
/// With `max_pieces <= 1` a single sequential positional read is used, which
/// lets the kernel's readahead handle mechanical disks optimally.
#[cfg(unix)]
fn read_block_parallel(
    file: &File,
    offset: &mut u64,
    file_len: u64,
    buf: &mut [u8],
    max_pieces: usize,
) -> Result<usize> {
    let remaining = file_len.saturating_sub(*offset) as usize;
    let want = buf.len().min(remaining);

    if want == 0 {
        return Ok(0);
    }

    let dst = &mut buf[..want];
    let base = *offset;
    let pieces = max_pieces.min(dst.len() / (256 << 10)).max(1);

    if pieces <= 1 {
        read_exact_at(file, dst, base)?;
    } else {
        let piece_len = dst.len().div_ceil(pieces);
        read_pool().install(|| -> Result<()> {
            let results: Vec<io::Result<()>> = dst
                .par_chunks_mut(piece_len)
                .enumerate()
                .map(|(i, chunk)| read_exact_at(file, chunk, base + (i * piece_len) as u64))
                .collect();

            for result in results {
                result?;
            }

            Ok(())
        })?;
    }

    *offset += want as u64;
    Ok(want)
}

#[cfg(unix)]
fn read_exact_at(file: &File, mut buf: &mut [u8], mut offset: u64) -> io::Result<()> {
    use std::os::unix::fs::FileExt;

    while !buf.is_empty() {
        match file.read_at(buf, offset) {
            Ok(0) => {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "file truncated while reading",
                ));
            }
            Ok(n) => {
                buf = &mut buf[n..];
                offset += n as u64;
            }
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }

    Ok(())
}

#[cfg(not(unix))]
fn read_block_parallel(
    file: &File,
    offset: &mut u64,
    file_len: u64,
    buf: &mut [u8],
    _max_pieces: usize,
) -> Result<usize> {
    use std::io::{Read, Seek};

    let remaining = file_len.saturating_sub(*offset) as usize;
    let want = buf.len().min(remaining);
    if want == 0 {
        return Ok(0);
    }

    let mut f = file;
    f.seek(SeekFrom::Start(*offset))?;
    let n = read_up_to(&mut f, &mut buf[..want])?;
    *offset += n as u64;
    Ok(n)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::StatOptions;

    fn config(all: bool) -> Config {
        Config {
            stats: StatOptions::new(all, 33, b"- ."),
            nx: vec![10.0, 50.0, 90.0],
            stream: StreamConfig::new(32 << 20, None, false),
        }
    }

    /// Multi-line FASTA with mixed case, N, gaps and variable lengths.
    pub(super) fn fasta_body(records: usize) -> Vec<u8> {
        let mut out = String::new();
        for i in 0..records {
            out.push_str(&format!(">seq{i} description {i}\n"));
            let len = 60 + (i * 37) % 400;
            let seq: String = (0..len)
                .map(|j| match (i + j) % 13 {
                    0 => 'N',
                    1 => 'n',
                    2 => 'g',
                    3 => 'c',
                    4 => '-',
                    5 => '.',
                    6 => 't',
                    _ => 'A',
                })
                .collect();
            for line in seq.as_bytes().chunks(60) {
                out.push_str(std::str::from_utf8(line).unwrap());
                out.push('\n');
            }
        }
        out.into_bytes()
    }

    pub(super) fn fastq_body(records: usize) -> Vec<u8> {
        let mut out = String::new();
        out.reserve(records * 300);
        for i in 0..records {
            let len = 50 + (i * 29) % 200;
            let seq: String = (0..len)
                .map(|j| match (i + j) % 5 {
                    0 => 'N',
                    1 => 'G',
                    2 => 'C',
                    _ => 'A',
                })
                .collect();
            let qual: String = (0..len)
                .map(|j| (b'!' + ((i + j) % 41) as u8) as char)
                .collect();
            out.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
        }
        out.into_bytes()
    }

    pub(super) fn assert_same(whole: &ChunkStats, chunked: &ChunkStats) {
        assert_eq!(whole.count, chunked.count);
        assert_eq!(whole.sum_len, chunked.sum_len);
        assert_eq!(whole.min_len, chunked.min_len);
        assert_eq!(whole.max_len, chunked.max_len);
        assert_eq!(whole.gc, chunked.gc);
        assert_eq!(whole.n_nuc, chunked.n_nuc);
        assert_eq!(whole.n_prot, chunked.n_prot);
        assert_eq!(whole.gaps, chunked.gaps);
        assert_eq!(whole.qual_hist, chunked.qual_hist);
        assert_eq!(whole.alphabet, chunked.alphabet);

        let mut a = whole.lengths.clone();
        let mut b = chunked.lengths.clone();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }

    #[test]
    fn fasta_chunking_matches_single_chunk() {
        let body = fasta_body(20_000);
        let opts = config(true);
        let whole = parse_fasta(&body, &opts.stats, true).unwrap();
        let chunked = parse_body(&body, Format::Fasta, &opts).unwrap();
        assert!(plan_chunks(&body, Format::Fasta).len() > 1);
        assert_same(&whole, &chunked);
    }

    #[test]
    fn fastq_chunking_matches_single_chunk() {
        let body = fastq_body(25_000);
        let opts = config(true);
        let whole = parse_fastq(&body, &opts.stats, true).unwrap();
        let chunked = parse_body(&body, Format::Fastq, &opts).unwrap();
        assert!(plan_chunks(&body, Format::Fastq).len() > 1);
        assert_same(&whole, &chunked);
    }

    #[test]
    fn fastq_chunking_with_crlf() {
        let lf = String::from_utf8(fastq_body(10_000)).unwrap();
        let crlf = lf.replace('\n', "\r\n").into_bytes();

        let opts = config(true);
        let whole = parse_fastq(&crlf, &opts.stats, true).unwrap();
        let chunked = parse_body(&crlf, Format::Fastq, &opts).unwrap();
        assert_same(&whole, &chunked);
    }

    #[test]
    fn fasta_boundaries_start_at_definitions() {
        let body = fasta_body(5_000);
        let chunks = plan_chunks(&body, Format::Fasta);
        assert!(chunks.len() > 1);
        for &(start, end) in &chunks {
            assert_eq!(body[start], b'>');
            assert!(start == 0 || body[start - 1] == b'\n');
            assert!(end > start);
        }
        assert_eq!(chunks.last().unwrap().1, body.len());
    }

    #[test]
    fn fastq_boundaries_start_at_records() {
        let body = fastq_body(10_000);
        let chunks = plan_chunks(&body, Format::Fastq);
        assert!(chunks.len() > 1);
        for &(start, end) in &chunks {
            assert_eq!(body[start], b'@');
            assert!(start == 0 || body[start - 1] == b'\n');
            assert!(
                body[start..end]
                    .iter()
                    .filter(|&&b| b == b'\n')
                    .count()
                    .is_multiple_of(4)
            );
        }
        assert_eq!(chunks.last().unwrap().1, body.len());
    }

    #[test]
    fn detects_format_with_leading_comments() {
        let (format, offset) = detect_format(b"#comment\n\n>seq\nACGT\n").unwrap();
        assert_eq!(format, Some(Format::Fasta));
        assert_eq!(&b"#comment\n\n>seq\nACGT\n"[offset..offset + 1], b">");

        assert_eq!(detect_format(b"").unwrap(), (None, 0));
        assert_eq!(detect_format(b"\n \n").unwrap(), (None, 3));
        assert!(detect_format(b"ACGT\n").is_err());
    }
}

#[cfg(test)]
mod stream_tests {
    use super::tests::{assert_same, fasta_body, fastq_body};
    use super::*;
    use crate::stats::StatOptions;
    use std::io::Cursor;

    fn config(block_size: usize) -> Config {
        Config {
            stats: StatOptions::new(true, 33, b"- ."),
            nx: vec![],
            stream: StreamConfig::new(block_size, None, false),
        }
    }

    #[test]
    fn streaming_fasta_matches_bulk() {
        let body = fasta_body(5_000);
        let whole = parse_body(&body, Format::Fasta, &config(32 << 20)).unwrap();
        let (streamed, format) = stream_stats(Cursor::new(&body), &config(512), "test").unwrap();
        assert_eq!(format, Some(Format::Fasta));
        assert_same(&whole, &streamed);
    }

    #[test]
    fn streaming_fastq_matches_bulk() {
        let body = fastq_body(5_000);
        let whole = parse_body(&body, Format::Fastq, &config(32 << 20)).unwrap();
        let (streamed, format) = stream_stats(Cursor::new(&body), &config(700), "test").unwrap();
        assert_eq!(format, Some(Format::Fastq));
        assert_same(&whole, &streamed);
    }

    #[test]
    fn streaming_crlf_matches_bulk() {
        let lf = String::from_utf8(fastq_body(2_000)).unwrap();
        let crlf = lf.replace('\n', "\r\n").into_bytes();
        let whole = parse_body(&crlf, Format::Fastq, &config(32 << 20)).unwrap();
        let (streamed, _) = stream_stats(Cursor::new(&crlf), &config(333), "test").unwrap();
        assert_same(&whole, &streamed);
    }

    #[test]
    fn streaming_handles_tiny_blocks_and_leading_comments() {
        let mut body = b"# leading comment\n\n".to_vec();
        body.extend_from_slice(&fasta_body(100));

        let whole = parse_body(
            &body[body.iter().position(|&b| b == b'>').unwrap()..],
            Format::Fasta,
            &config(32 << 20),
        )
        .unwrap();
        let (streamed, format) = stream_stats(Cursor::new(&body), &config(64), "test").unwrap();
        assert_eq!(format, Some(Format::Fasta));
        assert_same(&whole, &streamed);
    }

    #[test]
    fn streaming_empty_and_comments_only() {
        let (stats, format) = stream_stats(Cursor::new(b""), &config(64), "test").unwrap();
        assert_eq!(format, None);
        assert_eq!(stats.count, 0);

        let (stats, format) = stream_stats(Cursor::new(b"# hi\n\n"), &config(64), "test").unwrap();
        assert_eq!(format, None);
        assert_eq!(stats.count, 0);
    }
}
