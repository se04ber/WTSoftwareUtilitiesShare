from __future__ import annotations

from pathlib import Path

import numpy as np


def _load_turb_table() -> np.ndarray:
    root = Path(__file__).resolve().parents[3]
    turb = root / "tests" / "fixtures" / "flow" / "wtref" / "UBA_BL_Rep2025_UW_01_turb.txt"
    assert turb.exists(), f"Missing flow fixture turb table: {turb}"
    arr = np.genfromtxt(turb, comments="#")
    assert arr.ndim == 2 and arr.shape[1] == 13
    return arr


def test_flow_turb_table_has_valid_ranges_for_plot_inputs():
    """
    Guardrails for the main arrays that feed turb-intensity and Lux plots.

    We validate the *data contract* (shape/ranges/monotonicity), not figures.
    """
    arr = _load_turb_table()

    z = arr[:, 3]
    u_over_uref = arr[:, 5]
    i_u = arr[:, 10]
    i_w = arr[:, 11]
    lux = arr[:, 12]

    assert np.all(np.isfinite(z))
    assert np.all(np.diff(z) > 0), "Expected strictly increasing profile heights"

    assert np.all(np.isfinite(u_over_uref))
    assert np.min(u_over_uref) > 0
    assert np.max(u_over_uref) < 2
    # Profiles can have small local non-monotonicity due to measurement noise.
    # We only guard against major issues (e.g. inverted order, wrong units).
    diffs = np.diff(u_over_uref)
    assert (u_over_uref[-1] - u_over_uref[0]) > 0.1, "Expected overall increase from bottom to top"
    assert np.sum(diffs < -1e-3) <= 2, "Too many large drops in U/Uref profile"

    for name, v in [("I_u", i_u), ("I_w", i_w)]:
        assert np.all(np.isfinite(v))
        assert np.min(v) >= 0
        assert np.max(v) <= 1

    assert np.all(np.isfinite(lux))
    assert np.min(lux) > 0


