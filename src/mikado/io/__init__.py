"""Input/output utilities."""

from mikado.io.fasta import read_fasta, write_fasta
from mikado.io.output import format_result, OutputFormat

__all__ = ["read_fasta", "write_fasta", "format_result", "OutputFormat"]
