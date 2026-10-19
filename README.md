# fast

用 Rust 实现的极速 FASTA/FASTQ 序列统计工具。
[`seqkit stats`](https://bioinf.shenwei.me/seqkit/usage/#stats)，并在此之上做了
**分块流式读取 + 记录边界对齐的并行处理**。默认每块 16 MiB，内存占用与文件大小
基本无关，同时充分利用多核 CPU。

核心依赖：

- [`noodles`](https://github.com/zaeleus/noodles) (`noodles-fasta` / `noodles-fastq`)：FASTA/FASTQ 记录解析
- [`rayon`](https://github.com/rayon-rs/rayon)：数据并行
- [`clap`](https://github.com/clap-rs/clap)：命令行界面
- [`rapidgzip-core`](https://github.com/COMBINE-lab/rapidgzip-rust)：多线程 gzip/BGZF 流式解压
- `memchr`：SIMD 字节扫描

## 特性

- 统计：序列数、总长度、最小/平均/最大长度、N50、N50_num、Q1/Q2/Q3、N10/N30/...、GC%、gap 数、N 数
- FASTQ 额外统计：Q20(%)、Q30(%)、AvgQual（与 seqkit 相同的误差率平均算法）
- 支持 `.fa/.fasta/.fq/.fastq`、gzip/BGZF 压缩文件、stdin
- 分块流式读取：每次只保留一个块和跨块尾记录，峰值内存约 20–50 MB，适合超大文件
- 普通文件的块使用多线程 `pread` 填充，避免 mmap 缺页抖动
- gzip/BGZF 使用 `rapidgzip-core` 并行解压（BGZF 块级并行；普通 gzip 使用 marker/window 算法自适应选择）
- 自动识别 FASTA/FASTQ 和 DNA/RNA/Protein/Unlimit 类型
- 多文件并行处理，输出顺序与输入一致
- 默认表格对齐输出（带千位分隔符）；`-T` 输出对齐的制表符表格

## 编译

```bash
cd fast
cargo build --release
# 或安装到 ~/.cargo/bin
cargo install --path .
```

## 用法

```bash
# 默认：序列数、长度信息 + N50
fast reads.fq.gz assembly.fa

# 全部统计（四分位、gap、N50_num、Q20/Q30/AvgQual、GC、N 数）
fast -a reads.fq.gz

# 对齐的制表符输出
fast -a -T reads.fq.gz

# 额外输出 N10/N75/N90 列
fast -a -N 10,75,90 -T reads.fq.gz

# 指定线程数（默认使用全部逻辑核心）
fast -t 8 -T reads.fq

# 调整流式块大小（默认 16M；越小内存越省，越大吞吐略高）
fast -B 8M -T reads.fq

# 多文件并发处理：文件之间并行，文件内部再按记录分块并行；
# 输出顺序与命令行参数顺序一致，单个文件失败不影响其他文件
fast -a -T sample1.fq.gz sample2.fq.gz sample3.fa

# basename、stdin
fast -b -T /data/*.fq.gz
cat reads.fq | fast -a -T
```

### 选项

```
Usage: fast [OPTIONS] [FILE]...

Arguments:
  [FILE]...  Input FASTA/FASTQ files (plain or gzip; "-" for stdin)

Options:
  -a, --all                  All statistics: quartiles, gaps, N50_num, quality, GC and N count
  -T, --tabular              Aligned tabular output (tab-separated, padded columns)
  -b, --basename             Only print the basename of each file
  -t, --threads <N>          Number of threads (default: all logical CPUs)
  -B, --block-size <SIZE>    Streaming read block size, e.g. 8M, 16M, 32M [default: 16M]
      --read-threads <N>     Plain-file read concurrency per block (0 = sequential, default: auto)
  -N, --N <LIST>             Additional Nx statistics, e.g. -N 10,75,90
      --phred-offset <OFFSET>  FASTQ quality offset (33 or 64) [default: 33]
  -G, --gap-letters <CHARS>  Gap letters used by -a [default: "- ."]
  -h, --help                 Print help
  -V, --version              Print version
```

## 输出列

默认（`fast`）：

| 列 | 说明 |
| --- | --- |
| `file` | 文件名（`-b` 时只显示 basename） |
| `format` | `FASTA` / `FASTQ`（空文件为空） |
| `type` | `DNA` / `RNA` / `Protein` / `Unlimit` |
| `num_seqs` | 序列条数 |
| `sum_len` | 总长度 |
| `min_len` | 最短序列长度 |
| `avg_len` | 平均长度（保留 1 位小数） |
| `max_len` | 最长序列长度 |
| `N50` | N50 |

`-a` 在默认列基础上增加（列顺序与 `seqkit stats -a` 一致）：

`Q1, Q2, Q3, sum_gap, N50, N50_num, Q20(%), Q30(%), AvgQual, GC(%), sum_n`

- `N50_num`（即 L50）：达到 N50 阈值所需的序列条数
- `sum_gap`：gap 字符数（默认 `-`、空格、`.`）
- `sum_n`：核苷酸序列统计 `N/n`，蛋白序列统计 `X/x`（与 seqkit 一致）
- `Q20/Q30/AvgQual`：FASTA 恒为 0

## 与 seqkit 的统计口径

已经逐项核对：

- 长度包含 gap 字符；N50 以 `sum_len / 2` 为阈值
- 四分位采用 seqkit 的取法：奇数条数时中位数同时属于上下半区；
  偶数条数时正好对半切分
- `avg_len` 四舍五入（half away from zero）到 1 位小数
- 表格模式中 `Q1/Q2/Q3` 在 `-T` 下按 `%.0f`（round half to even）显示，
  非表格模式显示精确值（如 `1.5`）
- `GC(%) = (G + C) / sum_len`，仅对核苷酸类型计算
- `AvgQual = -10 * log10(mean(10^(-Q/10)))`
- 类型识别取第一条序列的前 10,000 bp，使用与 shenwei356/bio 相同的字母表判断

在测试数据（DNA/RNA/Protein/Unlimit、大小写、N、gap、CRLF、gzip、空文件、
单条超长序列等）上，`fast -a -T` 与 `seqkit stats -a -T` 的数值列完全一致
（仅全零长度记录的极端情况除外：seqkit 会输出 NaN/整数溢出，fast 输出 0）。

## 性能

测试环境：20 vCPU（Intel i9-10900X），Linux，page cache 热数据，输出重定向到
`/dev/null`，取 2–3 次运行中的最好成绩。默认块大小 16 MiB、20 线程。

### FASTQ（951 MB，240 万条 reads，100–300 bp）

| 工具 | 默认统计 | 全统计 `-a` |
| --- | ---: | ---: |
| `seqkit stats -T` | 0.77 s | 3.60 s |
| `fast -T` | **0.24 s**（RSS ≈ 46 MB） | **0.36 s**（RSS ≈ 46 MB） |
| 加速比 | **≈ 3.2×** | **≈ 10×** |

线程扩展（默认统计）：`-t 1` 0.51 s、`-t 4` 0.20 s、`-t 8` 0.20 s、
`-t 20` 0.24 s。8 线程附近吞吐最佳（内存带宽/超线程限制）。

### FASTA（546 MB，22 万条序列，50–5000 bp）

| 工具 | 默认统计 | 全统计 `-a` |
| --- | ---: | ---: |
| `seqkit stats -T` | 0.22 s | 1.38 s |
| `fast -T` | **0.13 s**（RSS ≈ 37 MB） | **0.18 s**（RSS ≈ 37 MB） |
| 加速比 | **≈ 1.7×** | **≈ 7.7×** |

### gzip（534 MB 压缩 → 951 MB 解压，951 MB FASTQ）

| 工具 / 输入 | 时间 | 峰值 RSS |
| --- | ---: | ---: |
| `seqkit stats -T`（普通 gzip） | 5.3 s | — |
| `fast -T`（普通单 member gzip） | **2.16 s** | ≈ 62 MB |
| `fast -T`（BGZF，块级并行） | **0.38 s** | ≈ 62 MB |

gzip 使用 `rapidgzip-core` 流式解码：BGZF 按块直接并行；普通 gzip 会根据
DEFLATE 块结构自适应选择 marker/window 并行或顺序解码。解压出的数据按块交给
解析器，不再把整个解压结果放进内存。

### 超大文件（冷缓存）

20 GB 的文件可以装进内存形成页缓存，重复运行几乎是纯 CPU 工作，多线程效果明显；
58 GB 这类超出内存的文件每次都要从磁盘读，瓶颈是磁盘带宽而不是 CPU。

对 951 MB 文件用 `posix_fadvise(DONTNEED)` 清空缓存后的实测（本地盘约 216 MB/s）：

| 工具 | 时间 | CPU 利用率 |
| --- | ---: | ---: |
| `seqkit stats -T` | 4.37 s | 27% |
| `fast -T` | 4.37 s | 41% |

两者都到达磁盘上限；fast 的读线程与解析线程重叠，CPU 利用率更高，但墙钟时间
不可能低于磁盘读取时间。可以先确认磁盘吞吐：

```bash
dd if=big.fq of=/dev/null bs=1M
```

如果磁盘上限远低于预期，建议把文件放到更快的存储。默认使用 32 路 `pread`
（相邻切片），对 SSD、RAID 阵列和并行文件系统都能提供足够队列深度；实测某
HDD RAID 冷区单流只有 3 MB/s，4 路并发可达 13 MB/s。只有在并发读导致明显
寻道抖动时才需要 `--read-threads 0`（顺序读）或调小数值。
