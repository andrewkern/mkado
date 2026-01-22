"""Tests for plotting functionality."""

from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from mkado.analysis.mk_test import MKResult
from mkado.analysis.polarized import PolarizedMKResult
from mkado.io.plotting import create_volcano_plot


class TestVolcanoPlot:
    """Tests for volcano plot generation."""

    def test_basic_volcano_plot(self) -> None:
        """Test basic volcano plot generation with MKResult."""
        results = [
            ("gene1", MKResult(dn=10, ds=20, pn=5, ps=10, p_value=0.5, ni=1.0, alpha=0.0)),
            ("gene2", MKResult(dn=15, ds=10, pn=3, ps=8, p_value=0.01, ni=0.5, alpha=0.5)),
            ("gene3", MKResult(dn=8, ds=25, pn=6, ps=12, p_value=0.001, ni=1.5, alpha=-0.5)),
        ]

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "volcano.png"
            create_volcano_plot(results, output_path)
            assert output_path.exists()
            assert output_path.stat().st_size > 0

    def test_volcano_plot_with_polarized_results(self) -> None:
        """Test volcano plot with PolarizedMKResult."""
        results = [
            (
                "gene1",
                PolarizedMKResult(
                    dn_ingroup=10,
                    ds_ingroup=20,
                    pn_ingroup=5,
                    ps_ingroup=10,
                    dn_outgroup=3,
                    ds_outgroup=5,
                    dn_unpolarized=2,
                    ds_unpolarized=3,
                    p_value_ingroup=0.5,
                    ni_ingroup=1.0,
                    alpha_ingroup=0.0,
                ),
            ),
            (
                "gene2",
                PolarizedMKResult(
                    dn_ingroup=15,
                    ds_ingroup=10,
                    pn_ingroup=3,
                    ps_ingroup=8,
                    dn_outgroup=5,
                    ds_outgroup=7,
                    dn_unpolarized=1,
                    ds_unpolarized=2,
                    p_value_ingroup=0.01,
                    ni_ingroup=0.5,
                    alpha_ingroup=0.5,
                ),
            ),
        ]

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "volcano_polarized.png"
            create_volcano_plot(results, output_path)
            assert output_path.exists()
            assert output_path.stat().st_size > 0

    def test_volcano_plot_skips_invalid_ni(self) -> None:
        """Test that genes with NI=0 or NI=None are skipped."""
        results = [
            ("gene1", MKResult(dn=10, ds=20, pn=5, ps=10, p_value=0.5, ni=1.0, alpha=0.0)),
            ("gene2", MKResult(dn=15, ds=10, pn=0, ps=8, p_value=1.0, ni=0.0, alpha=1.0)),
            ("gene3", MKResult(dn=8, ds=25, pn=6, ps=0, p_value=1.0, ni=None, alpha=None)),
        ]

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "volcano.png"
            create_volcano_plot(results, output_path)
            assert output_path.exists()

    def test_volcano_plot_no_valid_data(self) -> None:
        """Test that ValueError is raised when no valid data."""
        results = [
            ("gene1", MKResult(dn=10, ds=20, pn=0, ps=10, p_value=1.0, ni=0.0, alpha=1.0)),
            ("gene2", MKResult(dn=15, ds=10, pn=0, ps=8, p_value=1.0, ni=None, alpha=None)),
        ]

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "volcano.png"
            with pytest.raises(ValueError, match="No valid NI/p-value pairs"):
                create_volcano_plot(results, output_path)

    def test_volcano_plot_different_formats(self) -> None:
        """Test saving volcano plot in different formats."""
        results = [
            ("gene1", MKResult(dn=10, ds=20, pn=5, ps=10, p_value=0.5, ni=1.0, alpha=0.0)),
            ("gene2", MKResult(dn=15, ds=10, pn=3, ps=8, p_value=0.01, ni=0.5, alpha=0.5)),
        ]

        with TemporaryDirectory() as tmpdir:
            for ext in ["png", "pdf", "svg"]:
                output_path = Path(tmpdir) / f"volcano.{ext}"
                create_volcano_plot(results, output_path)
                assert output_path.exists()
                assert output_path.stat().st_size > 0
