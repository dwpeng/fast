mod parser;
mod stats;

use std::path::Path;
use std::process::ExitCode;

use clap::Parser;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use unicode_width::UnicodeWidthStr;

use crate::parser::{Config, process_path};
use crate::stats::{StatOptions, Summary};

/// Fast FASTA/FASTQ summary statistics.
///
/// Output columns follow `seqkit stats`; N50 is included by default and the
/// full seqkit column set is available with `-a`.
#[derive(Parser, Debug)]
#[command(name = "fast", version, about, long_about = None)]
struct Cli {
    /// Input FASTA/FASTQ files (plain or gzip; "-" for stdin)
    #[arg(value_name = "FILE")]
    files: Vec<String>,

    /// All statistics: quartiles, gaps, N50_num, quality, GC and N count
    #[arg(short = 'a', long)]
    all: bool,

    /// Aligned tabular output (tab-separated, padded columns)
    #[arg(short = 'T', long)]
    tabular: bool,

    /// Plain-file read concurrency per block (0 = sequential; default: auto)
    #[arg(long = "read-threads", value_name = "N")]
    read_threads: Option<usize>,

    /// Print read throughput progress to stderr every 5 seconds
    #[arg(short = 'v', long)]
    verbose: bool,

    /// Only print the basename of each file
    #[arg(short = 'b', long)]
    basename: bool,

    /// Number of threads (default: all logical CPUs)
    #[arg(short = 't', long, value_name = "N")]
    threads: Option<usize>,

    /// Streaming read block size (e.g. 8M, 16M, 32M, 64M; minimum 1M)
    #[arg(
        short = 'B',
        long = "block-size",
        value_name = "SIZE",
        default_value = "16M",
        value_parser = parse_size
    )]
    block_size: usize,

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

fn main() -> ExitCode {
    let cli = Cli::parse();

    if let Some(n) = cli.threads {
        if n == 0 {
            eprintln!("fast: --threads must be at least 1");
            return ExitCode::FAILURE;
        }
        if let Err(e) = ThreadPoolBuilder::new().num_threads(n).build_global() {
            eprintln!("fast: failed to configure the thread pool: {e}");
            return ExitCode::FAILURE;
        }
    }

    if cli.phred_offset != 33 && cli.phred_offset != 64 {
        eprintln!("fast: --phred-offset must be 33 or 64");
        return ExitCode::FAILURE;
    }

    if cli
        .nx
        .iter()
        .any(|x| !x.is_finite() || !(0.0..=100.0).contains(x))
    {
        eprintln!("fast: -N values must be in the range [0, 100]");
        return ExitCode::FAILURE;
    }

    if !cli.gap_letters.is_ascii() {
        eprintln!("fast: --gap-letters must be ASCII");
        return ExitCode::FAILURE;
    }

    let cfg = Config {
        stats: StatOptions::new(cli.all, cli.phred_offset, cli.gap_letters.as_bytes()),
        nx: cli.nx.clone(),
        block_size: cli.block_size,
        read_threads: cli.read_threads,
        verbose: cli.verbose,
    };

    let files: Vec<String> = if cli.files.is_empty() {
        vec!["-".to_string()]
    } else {
        cli.files.clone()
    };

    // Bound the number of files processed at once so that file-level
    // parallelism cannot blow up memory (each file keeps up to two blocks).
    let budget = 256usize << 20;
    let per_file = cli.block_size.saturating_mul(2).max(1);
    let max_concurrent = (budget / per_file)
        .clamp(1, 32)
        .min(rayon::current_num_threads().max(1));

    let mut results: Vec<(String, anyhow::Result<Summary>)> = Vec::with_capacity(files.len());
    for batch in files.chunks(max_concurrent) {
        let batch_results: Vec<(String, anyhow::Result<Summary>)> = batch
            .par_iter()
            .map(|path| {
                let name = display_name(path, cli.basename);
                (name, process_path(path, &cfg))
            })
            .collect();
        results.extend(batch_results);
    }

    let mut rows: Vec<Vec<String>> = Vec::with_capacity(results.len());
    let mut failed = false;

    for (path, result) in &results {
        match result {
            Ok(summary) => rows.push(build_row(path, summary, &cli)),
            Err(e) => {
                failed = true;
                eprintln!("fast: {path}: {e:#}");
            }
        }
    }

    print_table(&cli, &rows);

    if failed {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
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

fn build_row(file: &str, s: &Summary, cli: &Cli) -> Vec<String> {
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

fn headers(cli: &Cli) -> Vec<String> {
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

fn print_table(cli: &Cli, rows: &[Vec<String>]) {
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
