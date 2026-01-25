from __future__ import annotations

from pathlib import Path

import numpy as np

import windtunnel as wt


def _load_profile() -> tuple[np.ndarray, np.ndarray]:
    root = Path(__file__).resolve().parents[3]
    turb = root / "tests" / "fixtures" / "flow" / "wtref" / "UBA_BL_Rep2025_UW_01_turb.txt"
    assert turb.exists(), f"Missing flow fixture turb table: {turb}"
    arr = np.genfromtxt(turb, comments="#")

    heights = arr[:, 3]  # Z_fs [m]
    u_over_uref = arr[:, 5]  # U/Uref [-]
    return u_over_uref, heights


def test_calc_alpha_matches_current_reference_value():
    """
    Regression guardrail: alpha fitting is easy to break when refactoring.

    This test locks down the current output on the shipped example profile.
    """
    u, z = _load_profile()
    alpha, ref = wt.calc_alpha(u, z, d0=0.0, BL_height=600.0, BL=[])

    np.testing.assert_allclose(alpha, 0.18882396819546157, rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(ref, 0.02409906757461394, rtol=0.0, atol=1e-6)


def test_calc_z0_matches_current_reference_value():
    """
    Regression guardrail for z0 estimation on the shipped example profile.
    """
    u, z = _load_profile()
    z0, err, fitted_height = wt.calc_z0(u, z, d0=0.0, sfc_height=100.0, sfc_layer=[])

    np.testing.assert_allclose(z0, 0.3048512663590977, rtol=0.0, atol=1e-6)
    np.testing.assert_allclose(err, 0.03180213715934185, rtol=0.0, atol=1e-6)
    assert isinstance(fitted_height, np.ndarray)
    assert fitted_height.size >= 3


