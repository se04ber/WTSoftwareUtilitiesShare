from __future__ import annotations

from pathlib import Path

import numpy as np

import windtunnel as wt
from support.concentration_uba_ga import get_paths, list_uba_ga_files, process_one_file
from support.config import UBA_GA_PREFIX


def test_parameters_per_folder_matches_per_file_for_single_file():
    """
    Guardrail for the `parameters_PerFolder=True` configuration path.

    The shipped UBA_GA dataset has per-file varying positions, so "per folder"
    is not equivalent for the full profile. But for a single file, the per-folder
    and per-file path should produce identical results if the ambient values match.
    """
    p = get_paths()
    file0 = list_uba_ga_files(p)[0]

    # Baseline: per-file ambient conditions from the normal CSV.
    pc_per_file = process_one_file(p, file0)

    per_folder_csv = (
        Path(__file__).resolve().parents[3]
        / "tests"
        / "fixtures"
        / "concentration"
        / "params"
        / "ambient_conditions_per_folder_example.csv"
    )
    assert per_folder_csv.exists()

    ambient = wt.PointConcentration.get_ambient_conditions(
        path=".", name=UBA_GA_PREFIX, input_file=str(per_folder_csv)
    )
    assert ambient is not None

    (
        x_source,
        y_source,
        z_source,
        x_measure,
        y_measure,
        z_measure,
        pressure,
        temperature,
        calibration_curve,
        mass_flow_controller,
        calibration_factor,
        scaling_factor,
        scale,
        ref_length,
        ref_height,
        gas_name,
        mol_weight,
        gas_factor,
        full_scale_wtref,
        full_scale_flow_rate,
        full_scale_temp,
        full_scale_pressure,
        config_name,
    ) = wt.PointConcentration.read_ambient_conditions(ambient, UBA_GA_PREFIX)

    pc_per_folder = wt.PointConcentration.from_file(str(p.input_dir / file0))
    pc_per_folder.ambient_conditions(
        x_source=x_source,
        y_source=y_source,
        z_source=z_source,
        x_measure=x_measure,
        y_measure=y_measure,
        z_measure=z_measure,
        pressure=pressure,
        temperature=temperature,
        calibration_curve=calibration_curve,
        mass_flow_controller=mass_flow_controller,
        calibration_factor=calibration_factor,
        config_name=config_name,
    )
    pc_per_folder.scaling_information(
        scaling_factor=scaling_factor,
        scale=scale,
        ref_length=ref_length,
        ref_height=ref_height,
    )
    pc_per_folder.tracer_information(
        gas_name=gas_name, mol_weight=mol_weight, gas_factor=gas_factor
    )
    pc_per_folder.full_scale_information(
        full_scale_wtref=full_scale_wtref,
        full_scale_flow_rate=full_scale_flow_rate,
        full_scale_temp=full_scale_temp,
        full_scale_pressure=full_scale_pressure,
    )

    pc_per_folder.convert_temperature()
    pc_per_folder.calc_wtref_mean()
    pc_per_folder.calc_model_mass_flow_rate(usingMaxFlowRate="True", applyCalibration="False")
    pc_per_folder.calc_net_concentration()
    pc_per_folder.calc_c_star()
    pc_per_folder.calc_full_scale_concentration()

    np.testing.assert_allclose(np.nanmean(pc_per_folder.net_concentration), np.nanmean(pc_per_file.net_concentration))
    np.testing.assert_allclose(np.nanmean(pc_per_folder.c_star), np.nanmean(pc_per_file.c_star))
    np.testing.assert_allclose(
        np.nanmean(pc_per_folder.full_scale_concentration), np.nanmean(pc_per_file.full_scale_concentration)
    )


