from __future__ import annotations

from pathlib import Path

import numpy as np

import windtunnel as wt


def _process_one(input_dir: Path, parameter_csv: Path, file_name: str):
    ambient = wt.PointConcentration.get_ambient_conditions(
        path=".", name=file_name, input_file=str(parameter_csv)
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
    ) = wt.PointConcentration.read_ambient_conditions(ambient, file_name)

    pc = wt.PointConcentration.from_file(str(input_dir / file_name))
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
        scaling_factor=scaling_factor, scale=scale, ref_length=ref_length, ref_height=ref_height
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


def test_multiple_prefixes_are_processed_and_combined_without_overwrites(tmp_path: Path):
    """
    Guardrail: ensure we can process multiple `namelist` prefixes and still keep
    results distinct (no overwrites) and combinable.
    """
    root = Path(__file__).resolve().parents[3]
    fixture = root / "tests" / "fixtures" / "concentration_multi"
    input_dir = fixture / "input"
    parameter_csv = fixture / "params" / "ambient_conditions_multi.csv"
    assert input_dir.exists()
    assert parameter_csv.exists()

    prefixes = ["UBA_GA_02_04_01_000_1_001", "UBA_GA_02_04_01_000_1_002"]
    all_files: list[str] = []
    pcs: dict[str, object] = {}

    for prefix in prefixes:
        files = wt.get_files(str(input_dir) + "/", prefix)
        assert files == [f"{prefix}.txt.ts#0", f"{prefix}.txt.ts#1"]
        all_files.extend(files)
        for f in files:
            pcs[f] = _process_one(input_dir, parameter_csv, f)

    assert len(all_files) == 4
    assert len(set(all_files)) == 4
    assert set(pcs.keys()) == set(all_files)

    # Save avg outputs into per-prefix folders to mimic real runs.
    avg_dirs = {}
    for prefix in prefixes:
        d = tmp_path / "Point_Data_avg" / prefix
        d.mkdir(parents=True, exist_ok=True)
        avg_dirs[prefix] = d

    for f in all_files:
        prefix = f.split(".txt.ts#")[0]
        pcs[f].save2file_avg(f, out_dir=str(avg_dirs[prefix]) + "/")

    from windtunnel.concentration.utils import combine_to_csv

    out_csv = tmp_path / "combined_multi.csv"
    df = combine_to_csv(
        ["_avg_" + f for f in all_files],
        base_path=str(tmp_path),
        file_type="avg",
        output_filename=str(out_csv),
    )
    assert df is not None
    assert out_csv.exists()
    assert len(df) == 4

    # Basic sanity: the two prefixes produce identical numeric results (since we copied the raw files),
    # but they remain distinct rows due to distinct filenames.
    assert df["Filename"].nunique() == 4
    cstar = df["Avg_c_star [-]"].to_numpy()
    assert np.isfinite(cstar).all()


