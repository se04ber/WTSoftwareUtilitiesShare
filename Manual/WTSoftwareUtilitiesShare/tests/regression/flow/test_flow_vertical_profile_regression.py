from __future__ import annotations

from pathlib import Path

import numpy as np

import windtunnel as wt
from support.config import (
    FLOW_DATA_ND,
    FLOW_NAME,
    FLOW_SCALE,
    FLOW_SCALE_FACTOR,
    FLOW_TS_FILES,
    FLOW_TURB_GOLDEN,
    FLOW_WTREF_FILE,
)


def _flow_fixture_paths():
    # tests/regression/flow/... -> repo root is 3 levels up
    root = Path(__file__).resolve().parents[3]
    base = root / "tests" / "fixtures" / "flow"
    return {
        "ts_dir": base / "Zeitserien",
        "wtref_dir": base / "wtref",
        "turb_golden": base / "wtref" / FLOW_TURB_GOLDEN,
        "wtref_file": base / "wtref" / FLOW_WTREF_FILE,
        "name": FLOW_NAME,
        "scale_factor": FLOW_SCALE_FACTOR,
        "scale": FLOW_SCALE,
        "data_nd": FLOW_DATA_ND,
    }


def _compute_row(
    ts_path: Path,
    wtref_dir: Path,
    name: str,
    index: int,
    scale_factor: float,
    scale: float,
    data_nd: int,
):
    ts = wt.Timeseries.from_file(str(ts_path))
    ts.get_wind_comps(str(ts_path))
    ts.get_wtref(str(wtref_dir) + "/", name, index=index, vscale=scale_factor)

    if data_nd == 0:
        ts.nondimensionalise()

    ts.adapt_scale(scale, Lref=1 / scale)
    ts.mask_outliers()
    ts.calc_equidistant_timesteps()

    wdir = float(ts.mean_direction)
    x = float(ts.x)
    y = float(ts.y)
    z = float(ts.z)
    wtref = float(ts.wtref)
    u_mean = float(np.mean(ts.u))
    v_mean = float(np.mean(ts.v))
    u_std = float(np.std(ts.u))
    v_std = float(np.std(ts.v))

    turb = wt.calc_turb_data(ts.u.dropna(), ts.v.dropna())
    flux = float(turb[2])
    I_u = float(turb[0])
    I_v = float(turb[1])

    dt = float(ts.t_eq[1] - ts.t_eq[0])
    if data_nd == 0:
        # In the legacy script, u_eq is non-dimensional; multiply with wtref to obtain dimensional u for Lux.
        lux = float(wt.calc_lux_data(dt, ts.u_eq.dropna().values * ts.wtref))
    else:
        lux = float(wt.calc_lux_data(dt, ts.u_eq.dropna().values))

    return np.array(
        [wdir, x, y, z, wtref, u_mean, v_mean, u_std, v_std, flux, I_u, I_v, lux], dtype=float
    )


def test_flow_fixture_matches_golden_turb_file_subset():
    p = _flow_fixture_paths()
    assert p["ts_dir"].exists()
    assert p["wtref_dir"].exists()
    assert p["turb_golden"].exists()
    assert p["wtref_file"].exists()

    expected = np.genfromtxt(p["turb_golden"], comments="#")
    assert expected.ndim == 2 and expected.shape[1] == 13

    for fname in FLOW_TS_FILES:
        ts_file = p["ts_dir"] / fname
        assert ts_file.exists()

        # Index in golden/wtref file is derived from suffix (000001 -> 0, ...).
        suffix = fname.split(".")[-2]  # "000001" from "...000001.txt"
        idx = int(suffix) - 1

        row = _compute_row(
            ts_path=ts_file,
            wtref_dir=p["wtref_dir"],
            name=p["name"],
            index=idx,
            scale_factor=p["scale_factor"],
            scale=p["scale"],
            data_nd=p["data_nd"],
        )

        np.testing.assert_allclose(row, expected[idx], rtol=0.0, atol=2e-5)


