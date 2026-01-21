"""Command-line interface for Mikado."""

from __future__ import annotations

from pathlib import Path

import click

from mikado import __version__
from mikado.analysis.asymptotic import asymptotic_mk_test
from mikado.analysis.mk_test import mk_test
from mikado.analysis.polarized import polarized_mk_test
from mikado.io.output import OutputFormat, format_batch_results, format_result


@click.group()
@click.version_option(version=__version__)
def main() -> None:
    """Mikado: McDonald-Kreitman test toolkit.

    A modern Python implementation for detecting selection using
    the McDonald-Kreitman test and related methods.
    """
    pass


@main.command()
@click.argument("ingroup", type=click.Path(exists=True, path_type=Path))
@click.argument("outgroup", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-p",
    "--polarize",
    type=click.Path(exists=True, path_type=Path),
    help="Second outgroup for polarized MK test",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["pretty", "tsv", "json"]),
    default="pretty",
    help="Output format",
)
@click.option(
    "-r",
    "--reading-frame",
    type=click.IntRange(1, 3),
    default=1,
    help="Reading frame (1, 2, or 3)",
)
def test(
    ingroup: Path,
    outgroup: Path,
    polarize: Path | None,
    output_format: str,
    reading_frame: int,
) -> None:
    """Run the McDonald-Kreitman test.

    INGROUP: Path to FASTA file with ingroup sequences (e.g., population sample)

    OUTGROUP: Path to FASTA file with outgroup sequences

    Examples:

        mikado test ingroup.fa outgroup.fa

        mikado test ingroup.fa outgroup.fa -p outgroup2.fa

        mikado test ingroup.fa outgroup.fa --format tsv
    """
    fmt = OutputFormat(output_format)

    if polarize:
        result = polarized_mk_test(
            ingroup=ingroup,
            outgroup1=outgroup,
            outgroup2=polarize,
            reading_frame=reading_frame,
        )
    else:
        result = mk_test(
            ingroup=ingroup,
            outgroup=outgroup,
            reading_frame=reading_frame,
        )

    click.echo(format_result(result, fmt))


@main.command()
@click.argument("ingroup", type=click.Path(exists=True, path_type=Path))
@click.argument("outgroup", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["pretty", "tsv", "json"]),
    default="pretty",
    help="Output format",
)
@click.option(
    "-r",
    "--reading-frame",
    type=click.IntRange(1, 3),
    default=1,
    help="Reading frame (1, 2, or 3)",
)
@click.option(
    "-b",
    "--bins",
    type=int,
    default=10,
    help="Number of frequency bins",
)
@click.option(
    "--bootstrap",
    type=int,
    default=100,
    help="Number of bootstrap replicates for CI",
)
def asymptotic(
    ingroup: Path,
    outgroup: Path,
    output_format: str,
    reading_frame: int,
    bins: int,
    bootstrap: int,
) -> None:
    """Run the asymptotic McDonald-Kreitman test.

    This method accounts for slightly deleterious mutations by examining
    how alpha changes across the frequency spectrum.

    Based on Messer & Petrov (2013) PNAS.

    INGROUP: Path to FASTA file with ingroup sequences

    OUTGROUP: Path to FASTA file with outgroup sequences

    Examples:

        mikado asymptotic ingroup.fa outgroup.fa

        mikado asymptotic ingroup.fa outgroup.fa --bins 20
    """
    fmt = OutputFormat(output_format)

    result = asymptotic_mk_test(
        ingroup=ingroup,
        outgroup=outgroup,
        reading_frame=reading_frame,
        num_bins=bins,
        bootstrap_replicates=bootstrap,
    )

    click.echo(format_result(result, fmt))


@main.command()
@click.argument("ingroup_dir", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--outgroup-pattern",
    default="*_outgroup.fa",
    help="Glob pattern to match outgroup files",
)
@click.option(
    "--ingroup-pattern",
    default="*_ingroup.fa",
    help="Glob pattern to match ingroup files",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["pretty", "tsv", "json"]),
    default="tsv",
    help="Output format",
)
@click.option(
    "-r",
    "--reading-frame",
    type=click.IntRange(1, 3),
    default=1,
    help="Reading frame (1, 2, or 3)",
)
@click.option(
    "-a",
    "--asymptotic",
    is_flag=True,
    help="Use asymptotic MK test instead of standard",
)
@click.option(
    "-p",
    "--polarize-pattern",
    default=None,
    help="Glob pattern for second outgroup (enables polarized MK test)",
)
@click.option(
    "-b",
    "--bins",
    type=int,
    default=10,
    help="Number of frequency bins (asymptotic only)",
)
@click.option(
    "--bootstrap",
    type=int,
    default=100,
    help="Bootstrap replicates for CI (asymptotic only)",
)
def batch(
    ingroup_dir: Path,
    outgroup_pattern: str,
    ingroup_pattern: str,
    output_format: str,
    reading_frame: int,
    asymptotic: bool,
    polarize_pattern: str | None,
    bins: int,
    bootstrap: int,
) -> None:
    """Run MK test on multiple gene files.

    Expects paired ingroup/outgroup files with matching prefixes.

    INGROUP_DIR: Directory containing FASTA files

    Examples:

        mikado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"

        mikado batch genes/ --asymptotic --bins 10

        mikado batch genes/ -p "*_outgroup2.fa"  # polarized MK test
    """
    fmt = OutputFormat(output_format)

    # Check for mutually exclusive options
    if asymptotic and polarize_pattern:
        click.echo("Error: --asymptotic and --polarize-pattern are mutually exclusive", err=True)
        return

    # Find all ingroup files
    ingroup_files = sorted(ingroup_dir.glob(ingroup_pattern))

    if not ingroup_files:
        click.echo(f"No files matching '{ingroup_pattern}' found in {ingroup_dir}", err=True)
        return

    results = []

    for ingroup_file in ingroup_files:
        # Construct outgroup filename
        # Try to match by replacing the ingroup pattern suffix
        base_name = ingroup_file.stem
        if "_ingroup" in base_name:
            outgroup_name = base_name.replace("_ingroup", "_outgroup") + ".fa"
        elif "_in" in base_name:
            outgroup_name = base_name.replace("_in", "_out") + ".fa"
        else:
            # Try finding matching outgroup in same directory
            outgroup_name = base_name + "_outgroup.fa"

        outgroup_file = ingroup_dir / outgroup_name

        if not outgroup_file.exists():
            # Try with the outgroup pattern
            potential_matches = list(ingroup_dir.glob(outgroup_pattern))
            # Find one with similar prefix
            prefix = base_name.split("_")[0]
            matches = [f for f in potential_matches if f.stem.startswith(prefix)]
            if matches:
                outgroup_file = matches[0]
            else:
                click.echo(f"Warning: No outgroup found for {ingroup_file.name}", err=True)
                continue

        # Find second outgroup for polarized test if specified
        outgroup2_file: Path | None = None
        if polarize_pattern:
            potential_matches = list(ingroup_dir.glob(polarize_pattern))
            prefix = base_name.split("_")[0]
            matches = [f for f in potential_matches if f.stem.startswith(prefix)]
            if matches:
                outgroup2_file = matches[0]
            else:
                click.echo(
                    f"Warning: No second outgroup found for {ingroup_file.name}", err=True
                )
                continue

        try:
            if asymptotic:
                result = asymptotic_mk_test(
                    ingroup=ingroup_file,
                    outgroup=outgroup_file,
                    reading_frame=reading_frame,
                    num_bins=bins,
                    bootstrap_replicates=bootstrap,
                )
            elif polarize_pattern and outgroup2_file:
                result = polarized_mk_test(
                    ingroup=ingroup_file,
                    outgroup1=outgroup_file,
                    outgroup2=outgroup2_file,
                    reading_frame=reading_frame,
                )
            else:
                result = mk_test(
                    ingroup=ingroup_file,
                    outgroup=outgroup_file,
                    reading_frame=reading_frame,
                )
            results.append((ingroup_file.stem, result))
        except Exception as e:
            click.echo(f"Error processing {ingroup_file.name}: {e}", err=True)

    if results:
        click.echo(format_batch_results(results, fmt))
    else:
        click.echo("No results to display", err=True)


@main.command()
@click.argument("fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-r",
    "--reading-frame",
    type=click.IntRange(1, 3),
    default=1,
    help="Reading frame (1, 2, or 3)",
)
def info(fasta: Path, reading_frame: int) -> None:
    """Display information about a FASTA file.

    Shows sequence count, alignment length, and codon statistics.

    FASTA: Path to FASTA file
    """
    from mikado.core.sequences import SequenceSet

    seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)

    click.echo(f"File: {fasta.name}")
    click.echo(f"Sequences: {len(seqs)}")

    if seqs.sequences:
        click.echo(f"Alignment length: {seqs.alignment_length} bp")
        click.echo(f"Codons: {seqs.num_codons}")
        click.echo(f"Reading frame: {reading_frame}")

        # Count polymorphic sites
        poly_sites = seqs.polymorphic_codons()
        click.echo(f"Polymorphic codons: {len(poly_sites)}")

        # List sequence names
        click.echo("\nSequences:")
        for seq in seqs.sequences:
            click.echo(f"  {seq.name} ({len(seq)} bp)")


if __name__ == "__main__":
    main()
