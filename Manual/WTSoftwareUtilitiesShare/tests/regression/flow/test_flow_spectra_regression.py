from __future__ import annotations

from pathlib import Path

import numpy as np

import windtunnel as wt
from support.config import FLOW_DATA_ND, FLOW_NAME, FLOW_SCALE, FLOW_SCALE_FACTOR, FLOW_SPECTRA_GOLDEN


def _fixture_paths():
    root = Path(__file__).resolve().parents[3]
    flow = root / "tests" / "fixtures" / "flow"
    return {
        "ts_dir": flow / "Zeitserien",
        "wtref_dir": flow / "wtref",
        "spectra": flow / "spectra" / FLOW_SPECTRA_GOLDEN,
    }


def test_calc_spectra_matches_saved_reference():
    """
    Regression test for spectra inputs used by spectra plots.

    We compare computed (f_sm, S_uu_sm, S_vv_sm, S_uv_sm) against the saved
    `spectra_*.txt` reference file.
    """
    p = _fixture_paths()
    assert p["spectra"].exists()

    # Pick the exact file the saved spectra is based on.
    ts_file = p["ts_dir"] / f"{FLOW_NAME}.000024.txt"
    assert ts_file.exists()

    # Index in wtref file is derived from suffix (000024 -> 23).
    idx = 24 - 1

    ts = wt.Timeseries.from_file(str(ts_file))
    ts.get_wind_comps(str(ts_file))
    ts.get_wtref(str(p["wtref_dir"]) + "/", FLOW_NAME, index=idx, vscale=FLOW_SCALE_FACTOR)

    if FLOW_DATA_ND == 0:
        ts.nondimensionalise()

    ts.adapt_scale(FLOW_SCALE, Lref=1 / FLOW_SCALE)
    ts.mask_outliers()
    ts.calc_equidistant_timesteps()

    # The legacy script uses *dimensional* wind values for spectra.
    u_dim = ts.u_eq.dropna() * ts.wtref
    v_dim = ts.v_eq.dropna() * ts.wtref
    t_eq = ts.t_eq[~np.isnan(ts.t_eq)]

    f_sm, s_uu, s_vv, s_uv, *_ = wt.calc_spectra(u_dim, v_dim, t_eq, ts.z)

    expected = np.genfromtxt(p["spectra"], comments="#")
    assert expected.ndim == 2 and expected.shape[1] == 4

    actual = np.vstack([f_sm, s_uu, s_vv, s_uv]).T
    assert actual.shape == expected.shape

    # Saved with fmt='%.8f'; allow tiny float noise.
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=5e-6)


