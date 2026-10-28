mod convert;
mod head;
mod merge;
mod parser;
mod stats;
mod util;

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand};
use rayon::ThreadPoolBuilder;
use unicode_width::UnicodeWidthStr;

use crate::parser::{Config, StreamConfig, process_path};
use crate::stats::{StatOptions, Summary};
use crate::util::is_broken_pipe;

/// Fast FASTA/FASTQ toolkit.
///
/// `fast stats` prints summary statistics (compatible with `seqkit stats`);
/// `fast head`, `fast fa2fq`, `fast fq2fa` and `fast merge` provide fast
/// record-level operations.
#[derive(Parser, Debug)]
#[command(name = "fast", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    #[command(flatten)]
    global: GlobalOpts,
}

/// Options shared by every subcommand.
#[derive(Args, Debug, Clone)]
struct GlobalOpts {
    /// Number of threads (default: all logical CPUs)
    #[arg(short = 't', long, value_name = "N", global = true)]
    threads: Option<usize>,

    /// Streaming read block size (e.g. 8M, 16M, 32M, 64M; minimum 1M)
    #[arg(
        short = 'B',
        long = "block-size",
        value_name = "SIZE",
        default_value = "16M",
        value_parser = parse_size,
        global = true
    )]
    block_size: usize,

    /// Plain-file read concurrency per block (0 = sequential; default: auto)
    #[arg(long = "read-threads", value_name = "N", global = true)]
    read_threads: Option<usize>,

    /// Print read throughput progress to stderr every 5 seconds
    #[arg(short = 'v', long, global = true)]
    verbose: bool,
}

#[derive(Args, Debug)]
struct StatsArgs {
    /// Input FASTA/FASTQ files (plain or gzip; "-" for stdin)
    #[arg(value_name = "FILE")]
    files: Vec<String>,

    /// All statistics: quartiles, gaps, N50_num, quality, GC and N count
    #[arg(short = 'a', long)]
    all: bool,

    /// Aligned tabular output (tab-separated, padded columns)
    #[arg(short = 'T', long)]
    tabular: bool,

    /// Only print the basename of each file
    #[arg(short = 'b', long)]
    basename: bool,

    /// Additional Nx statistics, e.g. -N 10,75,90
    #[arg(short = 'N', long = "N", value_name = "LIST", value_delimiter = ',')]
    nx: Vec<f64>,

    /// FASTQ quality offset (33 for Sanger/Illumina 1.8+, 64 for Illumina 1.3+)
    #[arg(long, value_name = "OFFSET", default_value_t = 33)]
    phred_offset: u8,

    /// Gap letters used by -a
    #[arg(
        short = 'G',
        long = "gap-letters",
        value_name = "CHARS",
        default_value = "- ."
    )]
    gap_letters: String,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Print summary statistics
    Stats(StatsArgs),

    /// Convert FASTA to FASTQ
    #[command(name = "fa2fq")]
    Fa2fq(Fa2fqArgs),

    /// Convert FASTQ to FASTA
    #[command(name = "fq2fa")]
    Fq2fa(Fq2faArgs),

    /// Merge (concatenate) FASTA/FASTQ files in input order
    Merge(MergeArgs),

    /// Print the first N records in total
    Head(HeadArgs),
}

#[derive(Args, Debug)]
struct Fa2fqArgs {
    /// Quality character for generated FASTQ (default 'I' = Phred 40)
    #[arg(short = 'q', long = "qual", value_name = "CHAR", default_value = "I")]
    qual: String,

    /// Output file (default: stdout; "-" for stdout)
    #[arg(short = 'o', long = "out-file", value_name = "FILE")]
    out_file: Option<String>,

    /// Input FASTA files (plain or gzip; "-" for stdin)
    #[arg(value_name = "FILE")]
    files: Vec<String>,
}

#[derive(Args, Debug)]
struct Fq2faArgs {
    /// Output file (default: stdout; "-" for stdout)
    #[arg(short = 'o', long = "out-file", value_name = "FILE")]
    out_file: Option<String>,

    /// Input FASTQ files (plain or gzip; "-" for stdin)
    #[arg(value_name = "FILE")]
    files: Vec<String>,
}

#[derive(Args, Debug)]
struct MergeArgs {
    /// Output file (default: stdout; "-" for stdout)
    #[arg(short = 'o', long = "out-file", value_name = "FILE")]
    out_file: Option<String>,

    /// Input FASTA/FASTQ files (plain or gzip; "-" for stdin)
    #[arg(value_name = "FILE")]
    files: Vec<String>,
}

#[derive(Args, Debug)]
struct HeadArgs {
    /// Number of records to print in total
    #[arg(
        short = 'n',
        long = "num",
        alias = "number",
        value_name = "N",
        default_value_t = 10
    )]
    num: u64,

    /// Output file (default: stdout; "-" for stdout)
    #[arg(short = 'o', long = "out-file", value_name = "FILE")]
    out_file: Option<String>,

    /// Input FASTA/FASTQ files (plain or gzip; "-" for stdin)
    #[arg(value_name = "FILE")]
    files: Vec<String>,
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    if let Some(n) = cli.global.threads {
        if n == 0 {
            eprintln!("fast: --threads must be at least 1");
            return ExitCode::FAILURE;
        }
        if let Err(e) = ThreadPoolBuilder::new().num_threads(n).build_global() {
            eprintln!("fast: failed to configure the thread pool: {e}");
            return ExitCode::FAILURE;
        }
    }

    let stream = StreamConfig::new(
        cli.global.block_size,
        cli.global.read_threads,
        cli.global.verbose,
    );

    match cli.command {
        Command::Stats(stats) => run_stats(stats, &stream),
        Command::Fa2fq(args) => run_fa2fq(args, &stream),
        Command::Fq2fa(args) => run_fq2fa(args, &stream),
        Command::Merge(args) => run_merge(args, &stream),
        Command::Head(args) => run_head(args, &stream),
    }
}

fn run_stats(args: StatsArgs, stream: &StreamConfig) -> ExitCode {
    if args.phred_offset != 33 && args.phred_offset != 64 {
        eprintln!("fast stats: --phred-offset must be 33 or 64");
        return ExitCode::FAILURE;
    }

    if args
        .nx
        .iter()
        .any(|x| !x.is_finite() || !(0.0..=100.0).contains(x))
    {
        eprintln!("fast stats: -N values must be in the range [0, 100]");
        return ExitCode::FAILURE;
    }

    if !args.gap_letters.is_ascii() {
        eprintln!("fast stats: --gap-letters must be ASCII");
        return ExitCode::FAILURE;
    }

    let cfg = Config {
        stats: StatOptions::new(args.all, args.phred_offset, args.gap_letters.as_bytes()),
        nx: args.nx.clone(),
        stream: stream.clone(),
    };

    let files: Vec<String> = if args.files.is_empty() {
        vec!["-".to_string()]
    } else {
        args.files.clone()
    };

    // Bound the number of files processed at once so that file-level
    // parallelism cannot blow up memory.  Each input keeps at most three
    // full-size pipeline buffers for plain input; gzip inputs use smaller 4 MiB
    // blocks and therefore stay below this estimate.
    let budget = 256usize << 20;
    let per_file = stream.block_size.saturating_mul(3).max(1);
    let max_concurrent = (budget / per_file)
        .clamp(1, 32)
        .min(rayon::current_num_threads().max(1));

    let mut results: Vec<(String, anyhow::Result<Summary>)> = Vec::with_capacity(files.len());
    for batch in files.chunks(max_concurrent) {
        let cfg = &cfg;
        // Use OS threads for file-level parallelism: each file already uses the
        // rayon pool for block/chunk processing, and blocking a rayon worker on
        // a file's pipeline could starve the very tasks it waits for.
        let batch_results: Vec<(String, anyhow::Result<Summary>)> = std::thread::scope(|scope| {
            let handles: Vec<_> = batch
                .iter()
                .map(|path| {
                    let name = display_name(path, args.basename);
                    scope.spawn(move || (name, process_path(path, cfg)))
                })
                .collect();

            handles
                .into_iter()
                .map(|handle| handle.join().expect("stats worker panicked"))
                .collect()
        });
        results.extend(batch_results);
    }

    let mut rows: Vec<Vec<String>> = Vec::with_capacity(results.len());
    let mut failed = false;

    for (path, result) in &results {
        match result {
            Ok(summary) => rows.push(build_row(path, summary, &args)),
            Err(e) => {
                failed = true;
                eprintln!("fast stats: {path}: {e:#}");
            }
        }
    }

    print_table(&args, &rows);

    if failed {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
}

fn run_fa2fq(args: Fa2fqArgs, stream: &StreamConfig) -> ExitCode {
    let qual = match parse_qual(&args.qual) {
        Ok(q) => q,
        Err(e) => {
            eprintln!("fast fa2fq: {e}");
            return ExitCode::FAILURE;
        }
    };

    let files = input_files(&args.files);
    let mut out = match open_output(args.out_file.as_deref()) {
        Ok(out) => out,
        Err(e) => {
            eprintln!("fast fa2fq: {e:#}");
            return ExitCode::FAILURE;
        }
    };

    let mut failed = false;
    for path in &files {
        match convert::fa2fq(path, qual, stream, &mut *out) {
            Ok(()) => {}
            Err(e) if is_broken_pipe(&e) => return finish_output(out, false),
            Err(e) => {
                failed = true;
                eprintln!("fast fa2fq: {path}: {e:#}");
            }
        }
    }

    finish_output(out, failed)
}

fn run_fq2fa(args: Fq2faArgs, stream: &StreamConfig) -> ExitCode {
    let files = input_files(&args.files);
    let mut out = match open_output(args.out_file.as_deref()) {
        Ok(out) => out,
        Err(e) => {
            eprintln!("fast fq2fa: {e:#}");
            return ExitCode::FAILURE;
        }
    };

    let mut failed = false;
    for path in &files {
        match convert::fq2fa(path, stream, &mut *out) {
            Ok(()) => {}
            Err(e) if is_broken_pipe(&e) => return finish_output(out, false),
            Err(e) => {
                failed = true;
                eprintln!("fast fq2fa: {path}: {e:#}");
            }
        }
    }

    finish_output(out, failed)
}

fn run_merge(args: MergeArgs, stream: &StreamConfig) -> ExitCode {
    let files = input_files(&args.files);
    let mut out = match open_output(args.out_file.as_deref()) {
        Ok(out) => out,
        Err(e) => {
            eprintln!("fast merge: {e:#}");
            return ExitCode::FAILURE;
        }
    };

    match merge::run(&files, stream, &mut *out) {
        Ok(failed) => finish_output(out, failed),
        Err(e) if is_broken_pipe(&e) => finish_output(out, false),
        Err(e) => {
            eprintln!("fast merge: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn run_head(args: HeadArgs, stream: &StreamConfig) -> ExitCode {
    let files = input_files(&args.files);
    // For a small -n, do not read a full 16 MiB block before stopping.
    let head_stream = stream.with_block_size(head_block_size(stream.block_size, args.num));

    let mut out = match open_output(args.out_file.as_deref()) {
        Ok(out) => out,
        Err(e) => {
            eprintln!("fast head: {e:#}");
            return ExitCode::FAILURE;
        }
    };

    match head::run(&files, args.num, &head_stream, &mut *out) {
        Ok(failed) => finish_output(out, failed),
        Err(e) if is_broken_pipe(&e) => finish_output(out, false),
        Err(e) => {
            eprintln!("fast head: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn head_block_size(configured: usize, num: u64) -> usize {
    // A rough upper bound of 1 KiB per record keeps small `-n` values from
    // reading a full block.  The 1 MiB floor avoids tiny reads when a single
    // record does not fit and the reader has to grow its buffer.  Long records
    // simply span more blocks.
    let estimate = num.saturating_mul(1024).max(1 << 20);
    let estimate = usize::try_from(estimate).unwrap_or(usize::MAX);
    configured.min(estimate)
}

fn input_files(files: &[String]) -> Vec<String> {
    if files.is_empty() {
        vec!["-".to_string()]
    } else {
        files.to_vec()
    }
}

fn open_output(path: Option<&str>) -> anyhow::Result<Box<dyn Write + Send>> {
    match path {
        Some(path) if path != "-" => {
            let file =
                File::create(path).map_err(|e| anyhow::anyhow!("failed to create {path}: {e}"))?;
            Ok(Box::new(BufWriter::with_capacity(1 << 20, file)))
        }
        _ => Ok(Box::new(BufWriter::with_capacity(1 << 20, io::stdout()))),
    }
}

fn finish_output(mut out: Box<dyn Write + Send>, failed: bool) -> ExitCode {
    match out.flush() {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => return ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("fast: failed to write output: {e}");
            return ExitCode::FAILURE;
        }
    }

    if failed {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
}

fn parse_qual(s: &str) -> Result<u8, String> {
    let bytes = s.as_bytes();
    if bytes.len() != 1 {
        return Err("--qual must be a single ASCII character".to_string());
    }

    let b = bytes[0];
    if !(33..=126).contains(&b) {
        return Err("--qual must be a printable ASCII character (33..=126)".to_string());
    }

    Ok(b)
}

/// Parse a human-readable byte size, e.g. `32M`, `8MiB`, `65536`.
fn parse_size(s: &str) -> Result<usize, String> {
    let t = s.trim().to_ascii_lowercase();

    let (digits, multiplier) = if let Some(x) = t.strip_suffix("kib") {
        (x, 1usize << 10)
    } else if let Some(x) = t.strip_suffix("kb") {
        (x, 1usize << 10)
    } else if let Some(x) = t.strip_suffix('k') {
        (x, 1usize << 10)
    } else if let Some(x) = t.strip_suffix("mib") {
        (x, 1usize << 20)
    } else if let Some(x) = t.strip_suffix("mb") {
        (x, 1usize << 20)
    } else if let Some(x) = t.strip_suffix('m') {
        (x, 1usize << 20)
    } else if let Some(x) = t.strip_suffix("gib") {
        (x, 1usize << 30)
    } else if let Some(x) = t.strip_suffix("gb") {
        (x, 1usize << 30)
    } else if let Some(x) = t.strip_suffix('g') {
        (x, 1usize << 30)
    } else if let Some(x) = t.strip_suffix('b') {
        (x, 1usize)
    } else {
        (t.as_str(), 1usize)
    };

    let value: usize = digits
        .trim()
        .parse()
        .map_err(|_| format!("invalid block size: {s:?}"))?;
    let bytes = value
        .checked_mul(multiplier)
        .ok_or_else(|| format!("block size is too large: {s:?}"))?;

    if bytes < (1 << 20) {
        return Err("block size must be at least 1M".to_string());
    }

    Ok(bytes)
}

fn display_name(path: &str, basename: bool) -> String {
    if !basename || path == "-" {
        return path.to_string();
    }

    Path::new(path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string())
}

fn build_row(file: &str, s: &Summary, cli: &StatsArgs) -> Vec<String> {
    let mut row = Vec::with_capacity(if cli.all {
        19 + cli.nx.len()
    } else {
        9 + cli.nx.len()
    });

    row.push(file.to_string());
    row.push(s.format.to_string());
    row.push(s.seq_type.clone());

    if cli.all {
        push_int(&mut row, s.num_seqs, cli.tabular);
        push_int(&mut row, s.sum_len, cli.tabular);
        push_int(&mut row, s.min_len, cli.tabular);
        row.push(format_avg_len(s.avg_len, cli.tabular));
        push_int(&mut row, s.max_len, cli.tabular);
        row.push(format_quartile(s.q1, cli.tabular));
        row.push(format_quartile(s.q2, cli.tabular));
        row.push(format_quartile(s.q3, cli.tabular));
        push_int(&mut row, s.sum_gap, cli.tabular);
        push_int(&mut row, s.n50, cli.tabular);
        push_int(&mut row, s.l50, cli.tabular);
        row.push(format_pct(s.q20, cli.tabular));
        row.push(format_pct(s.q30, cli.tabular));
        row.push(format_float(s.avg_qual, 2, cli.tabular));
        row.push(format_float(s.gc, 2, cli.tabular));
        push_int(&mut row, s.sum_n, cli.tabular);
    } else {
        push_int(&mut row, s.num_seqs, cli.tabular);
        push_int(&mut row, s.sum_len, cli.tabular);
        push_int(&mut row, s.min_len, cli.tabular);
        row.push(format_avg_len(s.avg_len, cli.tabular));
        push_int(&mut row, s.max_len, cli.tabular);
        push_int(&mut row, s.n50, cli.tabular);
    }

    for (_, value) in &s.nx {
        push_int(&mut row, *value, cli.tabular);
    }

    row
}

fn push_int(row: &mut Vec<String>, value: u64, tabular: bool) {
    if tabular {
        row.push(value.to_string());
    } else {
        row.push(add_commas(&value.to_string()));
    }
}

fn format_avg_len(value: f64, tabular: bool) -> String {
    if tabular {
        format!("{value:.1}")
    } else {
        pretty_float(value, 1)
    }
}

fn format_quartile(value: f64, tabular: bool) -> String {
    if tabular {
        // seqkit uses printf("%.0f"), i.e. round half to even.
        format!("{value:.0}")
    } else {
        pretty_float(value, 1)
    }
}

fn format_pct(value: f64, tabular: bool) -> String {
    if tabular {
        // seqkit stores two decimals but prints Q20/Q30 with "%.0f".
        format!("{value:.0}")
    } else {
        pretty_float(value, 2)
    }
}

fn format_float(value: f64, decimals: usize, tabular: bool) -> String {
    if tabular {
        format!("{value:.decimals$}")
    } else {
        pretty_float(value, decimals)
    }
}

/// Format a float the way seqkit's humanized table does: at most `decimals`
/// fractional digits, trailing zeros removed, thousands separators added.
fn pretty_float(value: f64, decimals: usize) -> String {
    if !value.is_finite() {
        return "0".to_string();
    }

    let mut s = format!("{value:.decimals$}");
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }

    add_commas(&s)
}

fn add_commas(s: &str) -> String {
    let (sign, rest) = match s.strip_prefix('-') {
        Some(rest) => ("-", rest),
        None => ("", s),
    };

    let (int_part, frac_part) = match rest.split_once('.') {
        Some((int_part, frac_part)) => (int_part, Some(frac_part)),
        None => (rest, None),
    };

    let mut out = String::with_capacity(s.len() + int_part.len() / 3 + 1);
    out.push_str(sign);

    let n = int_part.len();
    for (i, ch) in int_part.char_indices() {
        if i > 0 && (n - i) % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }

    if let Some(frac) = frac_part {
        out.push('.');
        out.push_str(frac);
    }

    out
}

fn nx_header(x: f64) -> String {
    if x.fract() == 0.0 {
        format!("N{}", x as u64)
    } else {
        let mut s = format!("{x}");
        if s.contains('.') {
            while s.ends_with('0') {
                s.pop();
            }
            if s.ends_with('.') {
                s.pop();
            }
        }
        format!("N{s}")
    }
}

fn headers(cli: &StatsArgs) -> Vec<String> {
    let mut headers: Vec<String> = if cli.all {
        [
            "file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len", "Q1",
            "Q2", "Q3", "sum_gap", "N50", "N50_num", "Q20(%)", "Q30(%)", "AvgQual", "GC(%)",
            "sum_n",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect()
    } else {
        [
            "file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len", "N50",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect()
    };

    for &x in &cli.nx {
        headers.push(nx_header(x));
    }

    headers
}

fn print_table(cli: &StatsArgs, rows: &[Vec<String>]) {
    let headers = headers(cli);

    let left_aligned = [true, true, true];

    let mut widths: Vec<usize> = headers.iter().map(|h| h.as_str().width()).collect();

    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            widths[i] = widths[i].max(cell.as_str().width());
        }
    }

    // In tabular mode the separator is a tab; because every field is padded to
    // the same width in every row, each tab lands on the same column stop.
    let separator = if cli.tabular { "\t" } else { "  " };

    println!(
        "{}",
        render_line(&headers, &widths, &left_aligned, separator)
    );

    for row in rows {
        println!("{}", render_line(row, &widths, &left_aligned, separator));
    }
}

fn render_line(
    cells: &[String],
    widths: &[usize],
    left_aligned: &[bool],
    separator: &str,
) -> String {
    let mut line = String::new();

    for (i, cell) in cells.iter().enumerate() {
        if i > 0 {
            line.push_str(separator);
        }

        let width = widths[i];
        let padding = width.saturating_sub(cell.as_str().width());
        let align_left = left_aligned.get(i).copied().unwrap_or(false);

        if align_left {
            line.push_str(cell);
            for _ in 0..padding {
                line.push(' ');
            }
        } else {
            for _ in 0..padding {
                line.push(' ');
            }
            line.push_str(cell);
        }
    }

    line.trim_end().to_string()
}
