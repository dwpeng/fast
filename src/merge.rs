//! Merge (concatenate) FASTA/FASTQ inputs in command-line order.
//!
//! Inputs are streamed through the same parallel block reader / gzip decoder
//! used by `stats`, so plain files get parallel `pread` and compressed files
//! get parallel decompression.  Output is plain text; pipe it through `gzip`
//! when compressed output is needed.

use std::io::Write;

use anyhow::{Result, bail};

use crate::parser::{Format, StreamConfig, stream_path};
use crate::util::is_broken_pipe;

pub fn run(files: &[String], cfg: &StreamConfig, out: &mut dyn Write) -> Result<bool> {
    let mut expected: Option<Format> = None;
    let mut failed = false;

    for path in files {
        match merge_path(path, cfg, out, &mut expected) {
            Ok(()) => {}
            Err(e) => {
                if is_broken_pipe(&e) {
                    return Err(e);
                }
                failed = true;
                eprintln!("fast merge: {path}: {e:#}");
                // A format mismatch makes all following output meaningless.
                if e.downcast_ref::<FormatMismatch>().is_some() {
                    break;
                }
            }
        }
    }

    Ok(failed)
}

#[derive(Debug)]
struct FormatMismatch;

impl std::fmt::Display for FormatMismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("mixed FASTA and FASTQ inputs cannot be merged")
    }
}

impl std::error::Error for FormatMismatch {}

fn merge_path(
    path: &str,
    cfg: &StreamConfig,
    out: &mut dyn Write,
    expected: &mut Option<Format>,
) -> Result<()> {
    let mut file_format: Option<Format> = None;
    let mut last_byte: Option<u8> = None;

    stream_path(path, cfg, |body, format, _eof| {
        match file_format {
            None => {
                file_format = Some(format);
                if let Some(want) = *expected {
                    if want != format {
                        return Err(FormatMismatch.into());
                    }
                } else {
                    *expected = Some(format);
                }
            }
            Some(have) if have != format => {
                bail!(
                    "input format changed mid-file ({} -> {})",
                    have.label(),
                    format.label()
                );
            }
            _ => {}
        }

        out.write_all(body)?;
        if let Some(&b) = body.last() {
            last_byte = Some(b);
        }
        Ok(true)
    })?;

    if last_byte.is_some() && last_byte != Some(b'\n') {
        out.write_all(b"\n")?;
    }

    Ok(())
}
