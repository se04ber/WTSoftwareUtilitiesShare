from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import windtunnel as wt

from support.config import (
    UBA_GA_CONTROL_EXCEL,
    UBA_GA_EXPECTED_COMBINED_CSV,
    UBA_GA_FILES,
    UBA_GA_PARAMETER_CSV,
)


@dataclass(frozen=True)
class UbaGaExamplePaths:
    root: Path
    input_dir: Path
    parameter_csv: Path
    expected_combined_csv: Path
    excel_control: Path


def get_paths() -> UbaGaExamplePaths:
    root = Path(__file__).resolve().parents[2]
    base = root / "tests" / "fixtures" / "concentration"
    return UbaGaExamplePaths(
        root=root,
        input_dir=base / "input",
        parameter_csv=base / "params" / UBA_GA_PARAMETER_CSV,
        expected_combined_csv=base / "expected" / UBA_GA_EXPECTED_COMBINED_CSV,
        excel_control=base / "excel" / UBA_GA_CONTROL_EXCEL,
    )


def list_uba_ga_files(paths: UbaGaExamplePaths) -> list[str]:
    for f in UBA_GA_FILES:
        if not (paths.input_dir / f).exists():
            raise FileNotFoundError(f"Missing test input file: {paths.input_dir / f}")
    return list(UBA_GA_FILES)


def process_one_file(paths: UbaGaExamplePaths, file_name: str):
    ambient = wt.PointConcentration.get_ambient_conditions(
        path=".", name=file_name, input_file=str(paths.parameter_csv)
    )
    if ambient is None:
        raise RuntimeError("Failed to read ambient conditions CSV.")

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
    ) = wt.PointConcentration.read_ambient_conditions(ambient, file_name)

    pc = wt.PointConcentration.from_file(str(paths.input_dir / file_name))
    pc.ambient_conditions(
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
    pc.scaling_information(
        scaling_factor=scaling_factor,
        scale=scale,
        ref_length=ref_length,
        ref_height=ref_height,
    )
    pc.tracer_information(gas_name=gas_name, mol_weight=mol_weight, gas_factor=gas_factor)
    pc.full_scale_information(
        full_scale_wtref=full_scale_wtref,
        full_scale_flow_rate=full_scale_flow_rate,
        full_scale_temp=full_scale_temp,
        full_scale_pressure=full_scale_pressure,
    )

    pc.convert_temperature()
    pc.calc_wtref_mean()
    pc.calc_model_mass_flow_rate(usingMaxFlowRate="True", applyCalibration="False")
    pc.calc_net_concentration()
    pc.calc_c_star()
    pc.calc_full_scale_concentration()

    return pc


