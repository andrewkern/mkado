# Anopheles Batch Example

This example demonstrates batch processing of MK tests using aligned coding sequences from three Anopheles mosquito species.

## Data

Three genes with aligned codon sequences:
- `AGAP000010.fa` - 27 sequences (12 gamb, 12 afun, 3 amin)
- `AGAP000021.fa` - 27 sequences (12 gamb, 12 afun, 3 amin)
- `AGAP000041.fa` - 27 sequences (12 gamb, 12 afun, 3 amin)

Species:
- **gamb**: Anopheles gambiae (ingroup)
- **afun**: Anopheles funestus (outgroup)
- **amin**: Anopheles minimus (second outgroup for polarized tests)

## Example Commands

### Standard MK Test

```bash
# Batch MK test: gamb (ingroup) vs afun (outgroup)
mkado batch examples/anopheles_batch/ -i gamb -o afun

# TSV output for downstream analysis
mkado batch examples/anopheles_batch/ -i gamb -o afun -f tsv

# JSON output
mkado batch examples/anopheles_batch/ -i gamb -o afun -f json
```

### Asymptotic MK Test

```bash
# Aggregate asymptotic MK test (pools polymorphism across genes)
mkado batch examples/anopheles_batch/ -i gamb -o afun -a

# Per-gene asymptotic analysis
mkado batch examples/anopheles_batch/ -i gamb -o afun -a --per-gene

# Custom frequency bins
mkado batch examples/anopheles_batch/ -i gamb -o afun -a -b 5
```

### Polarized MK Test

```bash
# Use amin as second outgroup to polarize mutations
mkado batch examples/anopheles_batch/ -i gamb -o afun --polarize-match amin
```

### Single Gene Analysis

```bash
# Standard MK test on a single gene
mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun

# Asymptotic test on a single gene
mkado test examples/anopheles_batch/AGAP000010.fa -i gamb -o afun -a
```

## Expected Output

Running `mkado batch examples/anopheles_batch/ -i gamb -o afun` should produce output similar to:

```
Processing 3 alignment files...

Gene             Dn    Ds    Pn    Ps    p-value    NI       alpha
AGAP000010       XX    XX    XX    XX    X.XXX      X.XXXX   X.XXXX
AGAP000021       XX    XX    XX    XX    X.XXX      X.XXXX   X.XXXX
AGAP000041       XX    XX    XX    XX    X.XXX      X.XXXX   X.XXXX
```

## Notes

- Sequences are filtered by name pattern: `-i gamb` matches all sequences containing "gamb"
- The alignments are codon-aligned (in-frame); use `-r 1` (default) for reading frame
- For polarized tests, use a more distant outgroup (amin) to assign mutations to lineages
