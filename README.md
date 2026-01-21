# Mikado

A modern Python implementation of the McDonald-Kreitman test toolkit for detecting selection in molecular evolution.

## Features

- **Standard MK test**: Classic 2x2 contingency table with Fisher's exact test
- **Polarized MK test**: Uses a third outgroup to assign mutations to lineages
- **Asymptotic MK test**: Frequency-bin α estimates with exponential extrapolation (Messer & Petrov 2013)
- **Multiple output formats**: Pretty-print, TSV, and JSON
- **Batch processing**: Process multiple genes in one command

## Installation

```bash
# Clone the repository
git clone https://github.com/andrewkern/mikado.git
cd mikado

# Install with uv
uv sync

# Or install with pip
pip install .
```

## Quick Start

```bash
# Run the standard MK test
mikado test ingroup.fa outgroup.fa

# Run with JSON output
mikado test ingroup.fa outgroup.fa --format json

# Run polarized MK test with a second outgroup
mikado test ingroup.fa outgroup1.fa -p outgroup2.fa

# Run asymptotic MK test
mikado asymptotic ingroup.fa outgroup.fa

# Get file info
mikado info sequences.fa
```

## Combined File Mode

Mikado supports alignments where multiple species are in a single FASTA file, differentiated by sequence name patterns. Use `--ingroup-match` and `--outgroup-match` to filter sequences by substring matching on the header.

```bash
# Run MK test on a combined alignment file
# Sequences with "_agam_" in the name are ingroup
# Sequences with "_afun_" in the name are outgroup
mikado test combined.fa combined.fa --ingroup-match "_agam_" --outgroup-match "_afun_"

# Polarized test with three species in one file
mikado test combined.fa combined.fa \
  --ingroup-match "_agam_" \
  --outgroup-match "_afun_" \
  --polarize-match "_amin_"

# Asymptotic test with combined file
mikado asymptotic combined.fa combined.fa \
  --ingroup-match "_agam_" \
  --outgroup-match "_afun_"
```

## Example Output

```
$ mikado test kreitmanAdh.fa mauritianaAdh.fa

MK Test Results:
  Divergence:    Dn=6, Ds=8
  Polymorphism:  Pn=1, Ps=8
  Fisher's exact p-value: 0.176
  Neutrality Index (NI):  0.1667
  Alpha (α):              0.8333
```

## Commands

### `mikado test`

Run the standard or polarized McDonald-Kreitman test.

```bash
mikado test INGROUP OUTGROUP [OPTIONS]

Options:
  -p, --polarize PATH       Second outgroup for polarized MK test
  -f, --format TEXT         Output format: pretty, tsv, json
  -r, --reading-frame       Reading frame (1, 2, or 3)
  --ingroup-match TEXT      Filter ingroup by name pattern (combined file mode)
  --outgroup-match TEXT     Filter outgroup by name pattern (combined file mode)
  --polarize-match TEXT     Filter polarizing outgroup by name pattern
```

### `mikado asymptotic`

Run the asymptotic MK test (Messer & Petrov 2013).

```bash
mikado asymptotic INGROUP OUTGROUP [OPTIONS]

Options:
  -f, --format TEXT         Output format: pretty, tsv, json
  -r, --reading-frame       Reading frame (1, 2, or 3)
  -b, --bins INTEGER        Number of frequency bins (default: 10)
  --bootstrap INTEGER       Bootstrap replicates for CI (default: 100)
  --ingroup-match TEXT      Filter ingroup by name pattern (combined file mode)
  --outgroup-match TEXT     Filter outgroup by name pattern (combined file mode)
```

### `mikado batch`

Run MK test on multiple gene files. Supports standard, asymptotic, and polarized tests.

**Separate file mode** (default): Uses glob patterns to match ingroup/outgroup files.

```bash
mikado batch DIRECTORY [OPTIONS]

Options:
  --ingroup-pattern TEXT    Glob pattern for ingroup files
  --outgroup-pattern TEXT   Glob pattern for outgroup files
  -a, --asymptotic          Use asymptotic MK test
  -p, --polarize-pattern    Glob pattern for second outgroup (polarized test)
  -b, --bins INTEGER        Number of frequency bins (asymptotic only)
  --bootstrap INTEGER       Bootstrap replicates for CI (asymptotic only)
  -f, --format TEXT         Output format
  -r, --reading-frame       Reading frame (1, 2, or 3)
```

**Combined file mode**: Uses `--file-pattern` with `--ingroup-match` and `--outgroup-match` to process combined alignment files.

```bash
mikado batch DIRECTORY [OPTIONS]

Options (combined mode):
  --file-pattern TEXT       Glob pattern for combined alignment files
  --ingroup-match TEXT      Substring to match ingroup sequences
  --outgroup-match TEXT     Substring to match outgroup sequences
  --polarize-match TEXT     Substring to match polarizing outgroup
  -a, --asymptotic          Use asymptotic MK test
```

Examples:

```bash
# Standard batch MK test (separate files)
mikado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"

# Asymptotic batch MK test (separate files)
mikado batch genes/ -a --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"

# Polarized batch MK test (separate files)
mikado batch genes/ -p "*_outgroup2.fa" --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"

# Combined file batch mode
mikado batch alignments/ --file-pattern "*_aligned.fa" \
  --ingroup-match "_agam_" --outgroup-match "_afun_"

# Combined file with asymptotic test
mikado batch alignments/ -a --file-pattern "*_aligned.fa" \
  --ingroup-match "_agam_" --outgroup-match "_afun_"
```

### `mikado info`

Display information about a FASTA file.

```bash
mikado info FASTA [OPTIONS]

Options:
  -r, --reading-frame    Reading frame (1, 2, or 3)
```

## Python API

```python
from mikado import mk_test, asymptotic_mk_test, SequenceSet

# Run MK test
result = mk_test("ingroup.fa", "outgroup.fa")
print(f"Alpha: {result.alpha}")
print(f"P-value: {result.p_value}")

# Run asymptotic MK test
result = asymptotic_mk_test("ingroup.fa", "outgroup.fa")
print(f"Asymptotic Alpha: {result.alpha_asymptotic}")
print(f"95% CI: {result.ci_low} - {result.ci_high}")

# Work with sequences directly
seqs = SequenceSet.from_fasta("sequences.fa")
print(f"Polymorphic codons: {len(seqs.polymorphic_codons())}")

# Combined file mode - filter by sequence name
all_seqs = SequenceSet.from_fasta("combined.fa")
ingroup = all_seqs.filter_by_name("_gamb_")
outgroup = all_seqs.filter_by_name("_mel_")
result = mk_test(ingroup, outgroup)
```

## Interpretation

### Neutrality Index (NI)

- **NI = 1**: Neutral evolution
- **NI > 1**: Excess polymorphism (negative/purifying selection)
- **NI < 1**: Excess divergence (positive selection)


## Development

```bash
# Install dev dependencies
uv sync

# Run tests
uv run pytest

# Run linter
uv run ruff check src/

# Run formatter
uv run ruff format src/
```

## References

- McDonald, J. H., & Kreitman, M. (1991). Adaptive protein evolution at the Adh locus in Drosophila. Nature, 351(6328), 652-654.
- Messer, P. W., & Petrov, D. A. (2013). Frequent adaptation and the McDonald–Kreitman test. PNAS, 110(21), 8615-8620.

## License

MIT License
