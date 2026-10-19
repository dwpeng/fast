//! Sequence length / base / quality statistics.
//!
//! The definitions follow `seqkit stats` so that the two tools can be compared
//! column by column:
//!
//! * sequence lengths include gaps (unless they are removed before input),
//! * `N50` uses half of `sum_len` as the threshold,
//! * quartiles use the same element-index convention as seqkit (the lower/upper
//!   halves include the median element for odd record counts),
//! * `GC(%)` = (G + C) / sum_len,
//! * `AvgQual` = -10 * log10(mean error probability).

use anyhow::{Result, bail};
use rayon::prelude::*;

pub const PHRED_MAX: usize = 94;
pub const ALPHABET_GUESS_LEN: usize = 10_000;

/// Options shared by every parsing thread.
pub struct StatOptions {
    /// Compute the extra statistics (`-a`).
    pub all: bool,
    pub phred_offset: u8,
    /// Gap byte lookup.
    pub gap: [bool; 256],
}

impl StatOptions {
    pub fn new(all: bool, phred_offset: u8, gap_letters: &[u8]) -> Self {
        let mut gap = [false; 256];
        for &b in gap_letters {
            gap[b as usize] = true;
        }

        Self {
            all,
            phred_offset,
            gap,
        }
    }
}

/// Alphabet classification, mirroring shenwei356/bio.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Alphabet {
    Dna,
    Rna,
    DnaRedundant,
    RnaRedundant,
    Protein,
    Unlimit,
}

impl Alphabet {
    /// Label used in the `type` column.
    pub fn label(self) -> &'static str {
        match self {
            Alphabet::Dna | Alphabet::DnaRedundant => "DNA",
            Alphabet::Rna | Alphabet::RnaRedundant => "RNA",
            Alphabet::Protein => "Protein",
            Alphabet::Unlimit => "Unlimit",
        }
    }

    pub fn is_nucleotide(self) -> bool {
        matches!(
            self,
            Alphabet::Dna | Alphabet::Rna | Alphabet::DnaRedundant | Alphabet::RnaRedundant
        )
    }
}

#[inline]
fn is_dna(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' | b'N' | b'n' | b' ' | b'-' | b'.'
    )
}

#[inline]
fn is_rna(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C' | b'G' | b'U' | b'a' | b'c' | b'g' | b'u' | b'N' | b'n' | b' ' | b'-' | b'.'
    )
}

#[inline]
fn is_ambiguity(b: u8) -> bool {
    matches!(
        b,
        b'R' | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'r'
            | b'y'
            | b's'
            | b'w'
            | b'k'
            | b'm'
            | b'b'
            | b'd'
            | b'h'
            | b'v'
    )
}

#[inline]
fn is_dna_redundant(b: u8) -> bool {
    is_dna(b) || is_ambiguity(b)
}

#[inline]
fn is_rna_redundant(b: u8) -> bool {
    is_rna(b) || is_ambiguity(b)
}

/// Protein alphabet used by shenwei356/bio: all letters plus `* _ . -` and space.
#[inline]
fn is_protein(b: u8) -> bool {
    b.is_ascii_alphabetic() || matches!(b, b'*' | b'_' | b'.' | b'-' | b' ')
}

/// Guess the alphabet from (a prefix of) one sequence.
pub fn guess_alphabet(seq: &[u8]) -> Alphabet {
    let s = if seq.len() > ALPHABET_GUESS_LEN {
        &seq[..ALPHABET_GUESS_LEN]
    } else {
        seq
    };

    if s.is_empty() {
        return Alphabet::Unlimit;
    }
    if s.iter().all(|&b| is_dna(b)) {
        return Alphabet::Dna;
    }
    if s.iter().all(|&b| is_rna(b)) {
        return Alphabet::Rna;
    }
    if s.iter().all(|&b| is_dna_redundant(b)) {
        return Alphabet::DnaRedundant;
    }
    if s.iter().all(|&b| is_rna_redundant(b)) {
        return Alphabet::RnaRedundant;
    }
    if s.iter().all(|&b| is_protein(b)) {
        return Alphabet::Protein;
    }
    Alphabet::Unlimit
}

/// Per-thread (per-chunk) accumulators.
pub struct ChunkStats {
    pub count: u64,
    pub sum_len: u64,
    pub min_len: u64,
    pub max_len: u64,
    pub gc: u64,
    pub n_nuc: u64,
    pub n_prot: u64,
    pub gaps: u64,
    /// Quality-score histogram (index = decoded phred score).
    pub qual_hist: [u64; PHRED_MAX],
    pub lengths: Vec<u32>,
    /// Only the first chunk sets this (alphabet of the first record).
    pub alphabet: Option<Alphabet>,
}

impl Default for ChunkStats {
    fn default() -> Self {
        Self {
            count: 0,
            sum_len: 0,
            min_len: u64::MAX,
            max_len: 0,
            gc: 0,
            n_nuc: 0,
            n_prot: 0,
            gaps: 0,
            qual_hist: [0; PHRED_MAX],
            lengths: Vec::new(),
            alphabet: None,
        }
    }
}

impl ChunkStats {
    #[inline]
    pub fn add(&mut self, seq: &[u8], qual: Option<&[u8]>, opts: &StatOptions) -> Result<()> {
        let len = seq.len();
        if len > u32::MAX as usize {
            bail!("sequences longer than 4 GiB are not supported");
        }

        self.count += 1;
        self.sum_len += len as u64;
        self.min_len = self.min_len.min(len as u64);
        self.max_len = self.max_len.max(len as u64);
        self.lengths.push(len as u32);

        if opts.all {
            // One pass over the sequence: GC, N/n, X/x and gaps.
            let mut gc = 0u64;
            let mut n_nuc = 0u64;
            let mut n_prot = 0u64;
            let mut gaps = 0u64;
            for &b in seq {
                let f = BASE_CLASS[b as usize];
                gc += (f & 1) as u64;
                n_nuc += ((f >> 1) & 1) as u64;
                n_prot += ((f >> 2) & 1) as u64;
                gaps += opts.gap[b as usize] as u64;
            }
            self.gc += gc;
            self.n_nuc += n_nuc;
            self.n_prot += n_prot;
            self.gaps += gaps;

            if let Some(q) = qual {
                if q.len() != len {
                    bail!(
                        "quality length ({}) does not match sequence length ({len})",
                        q.len()
                    );
                }

                let offset = opts.phred_offset;

                for (i, &b) in q.iter().enumerate() {
                    if b < offset {
                        bail!(
                            "invalid quality byte {:?} at position {} (below ASCII base {offset})",
                            b as char,
                            i + 1
                        );
                    }
                    let v = (b - offset) as usize;
                    if v >= PHRED_MAX {
                        bail!(
                            "invalid quality byte {:?} at position {} (score {v})",
                            b as char,
                            i + 1
                        );
                    }
                    self.qual_hist[v] += 1;
                }
            }
        }

        Ok(())
    }

    pub fn merge(&mut self, other: ChunkStats) {
        if other.count == 0 {
            if self.alphabet.is_none() {
                self.alphabet = other.alphabet;
            }
            return;
        }

        self.count += other.count;
        self.sum_len += other.sum_len;
        self.min_len = self.min_len.min(other.min_len);
        self.max_len = self.max_len.max(other.max_len);
        self.gc += other.gc;
        self.n_nuc += other.n_nuc;
        self.n_prot += other.n_prot;
        self.gaps += other.gaps;
        for (dst, src) in self.qual_hist.iter_mut().zip(other.qual_hist) {
            *dst += src;
        }
        self.lengths.extend(other.lengths);
        if self.alphabet.is_none() {
            self.alphabet = other.alphabet;
        }
    }
}

/// Byte classes: bit 0 = GC, bit 1 = N/n, bit 2 = X/x.
static BASE_CLASS: [u8; 256] = build_base_class();

const fn build_base_class() -> [u8; 256] {
    let mut table = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        let mut flags = 0u8;
        if matches!(b, b'G' | b'C' | b'g' | b'c') {
            flags |= 1;
        }
        if matches!(b, b'N' | b'n') {
            flags |= 2;
        }
        if matches!(b, b'X' | b'x') {
            flags |= 4;
        }
        table[i] = flags;
        i += 1;
    }
    table
}

/// Final, display-ready summary for one file.
pub struct Summary {
    pub format: &'static str,
    pub seq_type: String,
    pub num_seqs: u64,
    pub sum_len: u64,
    pub min_len: u64,
    pub avg_len: f64,
    pub max_len: u64,
    pub q1: f64,
    pub q2: f64,
    pub q3: f64,
    pub sum_gap: u64,
    pub n50: u64,
    pub l50: u64,
    pub q20: f64,
    pub q30: f64,
    pub avg_qual: f64,
    pub gc: f64,
    pub sum_n: u64,
    pub nx: Vec<(f64, u64)>,
}

impl Summary {
    pub fn empty(format: &'static str, nx_values: &[f64]) -> Self {
        Self {
            format,
            seq_type: String::new(),
            num_seqs: 0,
            sum_len: 0,
            min_len: 0,
            avg_len: 0.0,
            max_len: 0,
            q1: 0.0,
            q2: 0.0,
            q3: 0.0,
            sum_gap: 0,
            n50: 0,
            l50: 0,
            q20: 0.0,
            q30: 0.0,
            avg_qual: 0.0,
            gc: 0.0,
            sum_n: 0,
            nx: nx_values.iter().map(|&x| (x, 0)).collect(),
        }
    }

    pub fn from_chunks(
        mut st: ChunkStats,
        format: &'static str,
        nx_values: &[f64],
        parallel_sort_threshold: usize,
    ) -> Self {
        if st.count == 0 {
            return Self::empty(format, nx_values);
        }

        if st.lengths.len() > parallel_sort_threshold {
            st.lengths.par_sort_unstable();
        } else {
            st.lengths.sort_unstable();
        }

        let alphabet = st.alphabet.unwrap_or(Alphabet::Unlimit);
        let seq_type = alphabet.label().to_string();

        let num_seqs = st.count;
        let sum_len = st.sum_len;
        let min_len = st.min_len;
        let max_len = st.max_len;
        let avg_len = round_half_away(sum_len as f64 / num_seqs as f64, 1);

        let (q1, q2, q3) = quartiles(&st.lengths);
        let (n50, l50) = n50_and_l50(&st.lengths, sum_len);
        let nx = nx_values
            .iter()
            .map(|&x| (x, nx_value(&st.lengths, sum_len, x)))
            .collect();

        let (q20, q30, avg_qual, gc) = if sum_len == 0 {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let q20_count: u64 = st.qual_hist[20..].iter().sum();
            let q30_count: u64 = st.qual_hist[30..].iter().sum();
            let err_sum: f64 = st
                .qual_hist
                .iter()
                .enumerate()
                .map(|(q, &count)| count as f64 * 10f64.powf(-(q as f64) / 10.0))
                .sum();

            let q20 = round_half_away(q20_count as f64 / sum_len as f64 * 100.0, 2);
            let q30 = round_half_away(q30_count as f64 / sum_len as f64 * 100.0, 2);
            let avg_qual = if err_sum == 0.0 {
                0.0
            } else {
                round_half_away(-10.0 * (err_sum / sum_len as f64).log10(), 2)
            };
            // seqkit only reports GC content for nucleotide alphabets.
            let gc = if alphabet.is_nucleotide() {
                round_half_away(st.gc as f64 / sum_len as f64 * 100.0, 2)
            } else {
                0.0
            };
            (q20, q30, avg_qual, gc)
        };

        let sum_n = if alphabet.is_nucleotide() {
            st.n_nuc
        } else {
            st.n_prot
        };

        Self {
            format,
            seq_type,
            num_seqs,
            sum_len,
            min_len,
            avg_len,
            max_len,
            q1,
            q2,
            q3,
            sum_gap: st.gaps,
            n50,
            l50,
            q20,
            q30,
            avg_qual,
            gc,
            sum_n,
            nx,
        }
    }
}

/// Round half away from zero (matches Go's `math.Round`).
#[inline]
pub fn round_half_away(x: f64, digits: i32) -> f64 {
    let p = 10f64.powi(digits);
    (x * p).round() / p
}

#[inline]
fn median(sorted: &[u32]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n.is_multiple_of(2) {
        (sorted[n / 2 - 1] as f64 + sorted[n / 2] as f64) / 2.0
    } else {
        sorted[n / 2] as f64
    }
}

/// Seqkit-compatible quartiles.
///
/// For an odd number of records the median record belongs to both halves; for
/// an even number the two halves are split exactly in the middle.
pub fn quartiles(sorted: &[u32]) -> (f64, f64, f64) {
    let n = sorted.len();
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }
    let lower_end = n.div_ceil(2);
    let upper_start = n / 2;
    (
        median(&sorted[..lower_end]),
        median(sorted),
        median(&sorted[upper_start..]),
    )
}

/// N50 and L50 (number of records needed to reach the N50 threshold).
pub fn n50_and_l50(sorted: &[u32], sum_len: u64) -> (u64, u64) {
    if sorted.is_empty() || sum_len == 0 {
        return (0, 0);
    }

    let mut cum = 0u64;
    for (i, &len) in sorted.iter().rev().enumerate() {
        cum += len as u64;
        if cum.saturating_mul(2) >= sum_len {
            return (len as u64, (i + 1) as u64);
        }
    }

    (0, 0)
}

/// Nx value (x in [0, 100]).
pub fn nx_value(sorted: &[u32], sum_len: u64, x: f64) -> u64 {
    if sorted.is_empty() || sum_len == 0 {
        return 0;
    }

    let boundary = sum_len as f64 * x / 100.0;
    let mut cum = 0.0f64;
    for &len in sorted.iter().rev() {
        cum += len as f64;
        if cum >= boundary {
            return len as u64;
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quartiles_match_seqkit() {
        assert_eq!(quartiles(&[2, 10, 10]), (6.0, 10.0, 10.0));
        assert_eq!(quartiles(&[6, 16]), (6.0, 11.0, 16.0));
        assert_eq!(quartiles(&[1, 2, 3, 4]), (1.5, 2.5, 3.5));
        assert_eq!(quartiles(&[1, 2, 3, 4, 5]), (2.0, 3.0, 4.0));
        assert_eq!(quartiles(&[1, 2, 3, 4, 5, 6]), (2.0, 3.5, 5.0));
    }

    #[test]
    fn n50_matches_seqkit() {
        assert_eq!(n50_and_l50(&[2, 10, 10], 22), (10, 2));
        assert_eq!(n50_and_l50(&[1, 2, 3, 4], 10), (3, 2));
        assert_eq!(n50_and_l50(&[6, 16], 22), (16, 1));
        assert_eq!(n50_and_l50(&[], 0), (0, 0));
    }

    #[test]
    fn alphabet_guessing() {
        assert_eq!(guess_alphabet(b"ACGTNN"), Alphabet::Dna);
        assert_eq!(guess_alphabet(b"ACGUUU"), Alphabet::Rna);
        assert_eq!(guess_alphabet(b"ACGTRYSWKMBDHVN"), Alphabet::DnaRedundant);
        assert_eq!(guess_alphabet(b"ACGTU"), Alphabet::Protein);
        assert_eq!(guess_alphabet(b"MEEPQSDPSV"), Alphabet::Protein);
        assert_eq!(guess_alphabet(b"acgtacgt"), Alphabet::Dna);
        assert_eq!(guess_alphabet(b"ACGT*ACGT"), Alphabet::Protein);
        assert_eq!(guess_alphabet(b""), Alphabet::Unlimit);
    }
}
