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
        .args(["-T", "a.fa", "b.fa", "c.fq"])
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
        .args(["-T", "missing.fa", "ok.fa"])
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
        .args(["-b", "-T", "sub/x.fa", "y.fa"])
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
        .args(["-T", "a.fa", "a-much-longer-file-name.fa"])
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
