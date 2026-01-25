from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from windtunnel.utils import get_files
from windtunnel.concentration.utils import calculate_uncertainties


def test_get_files_matches_prefix_and_sorts(tmp_path: Path):
    (tmp_path / "abc_2.txt").write_text("x")
    (tmp_path / "abc_1.txt").write_text("x")
    (tmp_path / "abcd.txt").write_text("x")
    (tmp_path / "zzz.txt").write_text("x")

    files = get_files(str(tmp_path), "abc")
    assert files == ["abc_1.txt", "abc_2.txt", "abcd.txt"]


@dataclass
class _FakeTs:
    config_name: str
    c_star: np.ndarray


def test_calculate_uncertainties_output_shape_is_stable():
    # Two files for the same config => finite uncertainty values.
    conc = {
        "a": _FakeTs(config_name="cfgA", c_star=np.array([1.0, 2.0, 3.0])),
        "b": _FakeTs(config_name="cfgA", c_star=np.array([1.5, 2.5, 3.5])),
    }
    out = calculate_uncertainties(
        conc,
        verbose=False,
        metrics_to_calculate=["Mean", "Median", "Peak2Mean", "P95"],
        concentration_types=["c_star"],
        include_abs=True,
        include_pct=True,
    )
    assert "cfgA" in out

    cfg = out["cfgA"]
    for k in [
        "Avg_c_star [-]_uncertainty_abs",
        "Avg_c_star [-]_uncertainty_pct",
        "Median_c_star_uncertainty_abs",
        "Median_c_star_uncertainty_pct",
        "Peak2MeanRatio_cstar_uncertainty_abs",
        "Peak2MeanRatio_cstar_uncertainty_pct",
        "Percentiles 95_cstar_uncertainty_abs",
        "Percentiles 95_cstar_uncertainty_pct",
    ]:
        assert k in cfg


