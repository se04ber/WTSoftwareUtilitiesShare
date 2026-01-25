from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np

import windtunnel as wt
from support.concentration_uba_ga import get_paths, list_uba_ga_files, process_one_file
from support.config import FLOW_NAME, FLOW_SCALE, FLOW_SCALE_FACTOR


def test_performance_smoke_budget():
    """
    Performance guardrail (CI-safe).

    This is NOT a micro-benchmark. It only protects against accidental slowdowns
    (e.g. adding an O(n^2) loop) by enforcing a generous time budget.

    You can override the threshold locally:
      WT_TEST_MAX_SECONDS=... pytest -q
    """
    max_seconds = float(os.getenv("WT_TEST_MAX_SECONDS", "60"))

    t0 = time.perf_counter()

    # ---- Concentration: process one representative file
    p = get_paths()
    file0 = list_uba_ga_files(p)[0]
    pc = process_one_file(p, file0)
    assert np.isfinite(np.nanmean(pc.c_star))

    # ---- Flow: compute spectra for one representative file
    root = Path(__file__).resolve().parents[3]
    ts_file = root / "tests" / "fixtures" / "flow" / "Zeitserien" / f"{FLOW_NAME}.000024.txt"
    wtref_dir = root / "tests" / "fixtures" / "flow" / "wtref"
    assert ts_file.exists()
    assert wtref_dir.exists()

    idx = 24 - 1
    ts = wt.Timeseries.from_file(str(ts_file))
    ts.get_wind_comps(str(ts_file))
    ts.get_wtref(str(wtref_dir) + "/", FLOW_NAME, index=idx, vscale=FLOW_SCALE_FACTOR)
    ts.nondimensionalise()
    ts.adapt_scale(FLOW_SCALE, Lref=1 / FLOW_SCALE)
    ts.mask_outliers()
    ts.calc_equidistant_timesteps()

    u_dim = ts.u_eq.dropna() * ts.wtref
    v_dim = ts.v_eq.dropna() * ts.wtref
    t_eq = ts.t_eq[~np.isnan(ts.t_eq)]
    f_sm, *_ = wt.calc_spectra(u_dim, v_dim, t_eq, ts.z)
    assert f_sm.size > 0

    elapsed = time.perf_counter() - t0
    assert elapsed <= max_seconds, f"Performance regression? took {elapsed:.2f}s (budget {max_seconds:.2f}s)"


