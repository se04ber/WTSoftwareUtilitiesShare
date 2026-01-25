from __future__ import annotations

from pathlib import Path

import numpy as np

from support.concentration_uba_ga import get_paths, list_uba_ga_files, process_one_file


def test_net_concentration_is_fast_minus_slow():
    p = get_paths()
    file = list_uba_ga_files(p)[0]
    pc = process_one_file(p, file)
    np.testing.assert_allclose(pc.net_concentration, pc.fast_FID - pc.slow_FID)


def test_c_star_formula_matches_implementation():
    p = get_paths()
    file = list_uba_ga_files(p)[0]
    pc = process_one_file(p, file)

    expected = ((pc.net_concentration / 1_000_000) * pc.wtref_mean * pc.ref_length**2) / (
        pc.mass_flow_rate / (1000 * 3600)
    )
    np.testing.assert_allclose(pc.c_star, expected)


def test_full_scale_concentration_formula_matches_implementation():
    p = get_paths()
    file = list_uba_ga_files(p)[0]
    pc = process_one_file(p, file)

    # `calc_full_scale_concentration()` mutates `full_scale_flow_rate` into [m^3/s].
    q_fs_m3s = pc.full_scale_flow_rate
    expected = pc.c_star * q_fs_m3s / (pc.full_scale_ref_length**2 * pc.full_scale_wtref) * 1_000_000
    np.testing.assert_allclose(pc.full_scale_concentration, expected)


def test_full_scale_and_non_dimensional_time_scaling():
    p = get_paths()
    file = list_uba_ga_files(p)[0]
    pc = process_one_file(p, file)

    model_time = pc.time.copy()
    fs_time = pc.calc_full_scale_time()
    nd_time = pc.calc_non_dimensional_time()

    np.testing.assert_allclose(fs_time, pc.scale * pc.wtref_mean / pc.full_scale_wtref * model_time)
    np.testing.assert_allclose(nd_time, pc.wtref_mean / pc.ref_length * model_time)


def test_save_ms_fs_nd_creates_expected_files(tmp_path: Path):
    p = get_paths()
    file = list_uba_ga_files(p)[0]

    pc_ms = process_one_file(p, file)
    pc_ms.save2file_ms(file, out_dir=str(tmp_path) + "/")
    assert (tmp_path / f"_ms_{file}").exists()

    pc_fs = process_one_file(p, file)
    pc_fs.to_full_scale()
    pc_fs.save2file_fs(file, out_dir=str(tmp_path) + "/")
    assert (tmp_path / f"_fs_{file}").exists()

    pc_nd = process_one_file(p, file)
    pc_nd.to_non_dimensional()
    pc_nd.save2file_nd(file, out_dir=str(tmp_path) + "/")
    assert (tmp_path / f"_nd_{file}").exists()


