# MKado Examples

This directory contains example data and tutorials for using mkado.

## Available Examples

### [anopheles_batch/](anopheles_batch/)

Batch processing example with Anopheles mosquito genes:
- 3 aligned coding sequences from Anopheles gambiae, A. funestus, and A. minimus
- Demonstrates batch MK test workflow with combined alignment files
- Includes both standard and asymptotic MK test examples

## Quick Start

```bash
# Install mkado
pip install mkado

# Run batch MK test (gamb vs afun)
mkado batch examples/anopheles_batch/ -i gamb -o afun

# Run asymptotic MK test with aggregated results
mkado batch examples/anopheles_batch/ -i gamb -o afun -a

# Run per-gene asymptotic analysis
mkado batch examples/anopheles_batch/ -i gamb -o afun -a --per-gene
```

## Data Sources

The example data is derived from aligned coding sequences of:
- **A. gambiae** (gamb) - African malaria mosquito
- **A. funestus** (afun) - African malaria vector
- **A. minimus** (amin) - Asian malaria vector (can be used as outgroup for polarized tests)
