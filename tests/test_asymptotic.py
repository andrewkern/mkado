"""Tests for asymptotic MK test."""

from pathlib import Path

import pytest

from mikado.analysis.asymptotic import AsymptoticMKResult, asymptotic_mk_test


class TestAsymptoticMKResult:
    """Tests for AsymptoticMKResult class."""

    def test_result_to_dict(self) -> None:
        """Test converting result to dictionary."""
        result = AsymptoticMKResult(
            frequency_bins=[0.1, 0.3, 0.5],
            alpha_by_freq=[0.1, 0.2, 0.3],
            alpha_asymptotic=0.4,
            ci_low=0.2,
            ci_high=0.6,
            fit_a=0.4,
            fit_b=0.3,
            fit_c=2.0,
            dn=10,
            ds=20,
        )

        d = result.to_dict()

        assert d["alpha_asymptotic"] == 0.4
        assert d["ci_low"] == 0.2
        assert d["ci_high"] == 0.6
        assert d["dn"] == 10
        assert d["ds"] == 20
        assert "fit_parameters" in d

    def test_result_str(self) -> None:
        """Test string representation."""
        result = AsymptoticMKResult(
            alpha_asymptotic=0.4,
            ci_low=0.2,
            ci_high=0.6,
            dn=10,
            ds=20,
        )

        s = str(result)
        assert "0.4" in s or "0.40" in s
        assert "Dn=10" in s


class TestAsymptoticMKTest:
    """Integration tests for asymptotic MK test."""

    def test_asymptotic_mk_test_basic(self, tmp_path: Path) -> None:
        """Test basic asymptotic MK test."""
        # Create test data with some variation
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
>seq4
ATGATGATGATGCTGATG
>seq5
ATGATGATGATGATGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        result = asymptotic_mk_test(ingroup_fa, outgroup_fa, num_bins=5)

        # Should complete without errors
        assert isinstance(result, AsymptoticMKResult)
        assert result.dn >= 0
        assert result.ds >= 0

    def test_asymptotic_mk_test_insufficient_data(self, tmp_path: Path) -> None:
        """Test asymptotic MK test with insufficient data."""
        # Create minimal test data
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATG
>seq2
ATGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGATG
""")

        result = asymptotic_mk_test(ingroup_fa, outgroup_fa)

        # Should handle gracefully
        assert isinstance(result, AsymptoticMKResult)

    def test_asymptotic_mk_test_with_kreitman_data(self) -> None:
        """Test asymptotic MK test with Kreitman data."""
        test_data = Path(__file__).parent / "data"
        ingroup = test_data / "kreitmanAdh.fa"
        outgroup = test_data / "mauritianaAdh.fa"

        if not ingroup.exists() or not outgroup.exists():
            pytest.skip("Test data files not found")

        result = asymptotic_mk_test(ingroup, outgroup, num_bins=5)

        assert isinstance(result, AsymptoticMKResult)
        # With only 4 ingroup sequences, might not have enough frequency data
        # but should still complete
