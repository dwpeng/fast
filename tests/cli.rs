use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fast")
}

fn tempdir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fast-cli-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn processes_multiple_files_in_input_order() {
    let dir = tempdir("multi");

    fs::write(dir.join("a.fa"), ">a1\nACGT\n>a2\nACGTAC\n").unwrap();
    fs::write(dir.join("b.fa"), ">b1\nAAAA\n").unwrap();
    fs::write(dir.join("c.fq"), "@r1\nACGT\n+\nIIII\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["stats", "-T", "a.fa", "b.fa", "c.fq"])
        .output()
        .unwrap();

    assert!(out.status.success());

    let text = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 4);

    // Header plus one row per input file, in the same order as the arguments.
    let fields =
        |line: &str| -> Vec<String> { line.split('\t').map(|f| f.trim().to_string()).collect() };
    assert!(
        lines[1].starts_with("a.fa")
            && fields(lines[1])[1..8].join(" ") == "FASTA DNA 2 10 4 5.0 6"
    );
    assert!(
        lines[2].starts_with("b.fa") && fields(lines[2])[1..8].join(" ") == "FASTA DNA 1 4 4 4.0 4"
    );
    assert!(
        lines[3].starts_with("c.fq") && fields(lines[3])[1..8].join(" ") == "FASTQ DNA 1 4 4 4.0 4"
    );

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn continues_after_a_missing_file() {
    let dir = tempdir("errors");

    fs::write(dir.join("ok.fa"), ">s\nACGT\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["stats", "-T", "missing.fa", "ok.fa"])
        .output()
        .unwrap();

    // The valid file is still reported, but the exit status is non-zero.
    assert!(!out.status.success());

    let stdout = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[1].starts_with("ok.fa"));

    let stderr = String::from_utf8(out.stderr).unwrap();
    assert!(stderr.contains("missing.fa"));

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn basename_applies_to_all_files() {
    let dir = tempdir("basename");

    fs::create_dir_all(dir.join("sub")).unwrap();
    fs::write(dir.join("sub/x.fa"), ">x\nACGT\n").unwrap();
    fs::write(dir.join("y.fa"), ">y\nACGT\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["stats", "-b", "-T", "sub/x.fa", "y.fa"])
        .output()
        .unwrap();

    assert!(out.status.success());

    let text = String::from_utf8(out.stdout).unwrap();
    let names: Vec<&str> = text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(names, ["x.fa", "y.fa"]);

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn tabular_output_is_aligned() {
    let dir = tempdir("align");

    fs::write(dir.join("a.fa"), ">a\nACGT\n").unwrap();
    fs::write(dir.join("a-much-longer-file-name.fa"), ">b\nACGTACGTAC\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["stats", "-T", "a.fa", "a-much-longer-file-name.fa"])
        .output()
        .unwrap();

    assert!(out.status.success());

    let text = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 3);

    // Every field is padded to the column width, so the byte offsets of the
    // tab-separated columns are identical on all rows.
    let column_offsets =
        |line: &str| -> Vec<usize> { line.match_indices('\t').map(|(i, _)| i).collect() };

    let header_offsets = column_offsets(lines[0]);
    for line in &lines[1..] {
        assert_eq!(column_offsets(line), header_offsets);
    }
}

#[test]
fn a_subcommand_is_required() {
    let dir = tempdir("subcommand-required");
    fs::write(dir.join("a.fa"), ">a\nACGT\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .arg("a.fa")
        .output()
        .unwrap();

    assert!(!out.status.success());
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert!(stderr.contains("subcommand") || stderr.contains("Usage"));

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn fa2fq_and_fq2fa_convert_records() {
    let dir = tempdir("convert");
    fs::write(dir.join("a.fa"), ">a desc\nACGT\nNN\n>b\nacgt\n").unwrap();
    fs::write(dir.join("a.fq"), "@r desc\nACGT\n+\nIIII\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["fa2fq", "a.fa"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "@a desc\nACGTNN\n+\nIIIIII\n@b\nacgt\n+\nIIII\n"
    );

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["fq2fa", "a.fq"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(String::from_utf8(out.stdout).unwrap(), ">r desc\nACGT\n");

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn fa2fq_handles_an_empty_sequence_before_another_record() {
    let dir = tempdir("fa2fq-empty");
    fs::write(dir.join("a.fa"), ">empty\n>full\nACGT\n>last\nTT\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["fa2fq", "a.fa"])
        .output()
        .unwrap();

    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "@empty\n\n+\n\n@full\nACGT\n+\nIIII\n@last\nTT\n+\nII\n"
    );

    let _ = fs::remove_dir_all(&dir);
}

/// Three wrapped FASTQ records: 100 bp sequence/quality split over two lines.
fn wrapped_fastq() -> String {
    let mut out = String::new();
    for i in 0..3 {
        let seq = "A".repeat(100);
        let qual = "I".repeat(100);
        out.push_str(&format!("@r{i} desc\n"));
        for line in seq.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
        out.push_str("+\n");
        for line in qual.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
    }
    out
}

#[test]
fn wrapped_fastq_is_supported_end_to_end() {
    let dir = tempdir("wrapped-fastq");
    let body = wrapped_fastq();
    fs::write(dir.join("a.fq"), &body).unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["stats", "-T", "a.fq"])
        .output()
        .unwrap();
    assert!(out.status.success());

    let text = String::from_utf8(out.stdout).unwrap();
    let fields: Vec<&str> = text.lines().nth(1).unwrap().split('\t').map(str::trim).collect();
    assert_eq!(fields[1..5].join(" "), "FASTQ DNA 3 300");

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["fq2fa", "a.fq"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let expected_seq = "A".repeat(100);
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        format!(">r0 desc\n{expected_seq}\n>r1 desc\n{expected_seq}\n>r2 desc\n{expected_seq}\n")
    );

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["head", "-n", "2", "a.fq"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(String::from_utf8(out.stdout).unwrap(), &body[..body.find("@r2").unwrap()]);

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn head_takes_the_first_n_records() {
    let dir = tempdir("head");
    fs::write(dir.join("a.fa"), ">a\nACGT\n>b\nACGTAC\n>c\nAA\n").unwrap();
    fs::write(
        dir.join("a.fq"),
        "@a\nACGT\n+\nIIII\n@b\nACGTAC\n+\nIIIIII\n@c\nAA\n+\nII\n",
    )
    .unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["head", "-n", "2", "a.fa"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        ">a\nACGT\n>b\nACGTAC\n"
    );

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["head", "-n", "2", "a.fq"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "@a\nACGT\n+\nIIII\n@b\nACGTAC\n+\nIIIIII\n"
    );

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn head_applies_n_across_all_inputs() {
    let dir = tempdir("head-multiple");
    fs::write(dir.join("a.fa"), ">a1\nA\n>a2\nC\n").unwrap();
    fs::write(dir.join("b.fa"), ">b1\nG\n>b2\nT\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["head", "-n", "3", "a.fa", "b.fa"])
        .output()
        .unwrap();

    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        ">a1\nA\n>a2\nC\n>b1\nG\n"
    );

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn merge_concatenates_in_order_and_rejects_mixed_formats() {
    let dir = tempdir("merge");
    fs::write(dir.join("a.fa"), ">a\nACGT\n").unwrap();
    fs::write(dir.join("b.fa"), ">b\nTTTT\n").unwrap();
    fs::write(dir.join("c.fq"), "@c\nACGT\n+\nIIII\n").unwrap();

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["merge", "a.fa", "b.fa"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        ">a\nACGT\n>b\nTTTT\n"
    );

    let out = Command::new(bin())
        .current_dir(&dir)
        .args(["merge", "a.fa", "c.fq"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert!(stderr.contains("mixed FASTA and FASTQ"));

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn global_options_work_on_either_side_of_the_subcommand() {
    let dir = tempdir("global-options");
    fs::write(dir.join("a.fa"), ">a\nACGT\n>b\nTTTT\n").unwrap();

    for args in [
        vec!["-t", "2", "head", "-n", "1", "a.fa"],
        vec!["head", "-n", "1", "-t", "2", "a.fa"],
    ] {
        let out = Command::new(bin())
            .current_dir(&dir)
            .args(args)
            .output()
            .unwrap();
        assert!(out.status.success());
        assert_eq!(String::from_utf8(out.stdout).unwrap(), ">a\nACGT\n");
    }

    let _ = fs::remove_dir_all(&dir);
}
