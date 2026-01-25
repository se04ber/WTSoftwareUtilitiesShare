from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import windtunnel as wt


@dataclass(frozen=True)
class PointConcentrationRunResult:
    conc_ts: dict
    files_by_name: dict[str, list[str]]
    namelist: list[str]


def _read_ambient_for_key(ambient_conditions, key: str):
    """
    `ambient_conditions` is a pandas DataFrame created by `PointConcentration.get_ambient_conditions`.
    We keep using the existing parser to avoid silently changing behavior.
    """
    return wt.PointConcentration.read_ambient_conditions(ambient_conditions, key)


def run_point_concentration(
    *,
    path: str,
    namelist: list[str],
    csv_file: str,
    parameters_per_folder: bool = False,
    full_scale: str = "ms",  # "ms" | "fs" | "nd"
    calculate_uncertainty: bool = False,
    split_factor: float = 1.8,
    uncertainty_threshold: float = 1e-4,
    uncertainty_metrics=None,
    uncertainty_concentration_types=None,
    include_abs_uncertainty: bool = True,
    include_pct_uncertainty: bool = True,
    save_all: bool = False,
    save_ts: bool = False,
    save_avg: bool = False,
    save_stats: bool = False,
    save_combined: bool = False,
    os_type: str = "Linux",
    output_path: Optional[str] = None,
    combined_filename: str = "combined_data.csv",
) -> PointConcentrationRunResult:
    """
    Compute-only pipeline for point concentration time series.

    Notes:
    - This intentionally does not do plotting.
    - Saving is optional; paths are kept compatible with existing scripts.
    """
    if full_scale not in {"ms", "fs", "nd"}:
        raise ValueError("full_scale must be one of: 'ms', 'fs', 'nd'")

    path = str(path)
    if not path.endswith("/"):
        path += "/"

    ambient_conditions = wt.PointConcentration.get_ambient_conditions(
        path=path, name=namelist[0] if namelist else None, input_file=csv_file
    )

    conc_ts: dict[str, dict[str, object]] = {}
    files_by_name: dict[str, list[str]] = {}

    for name in namelist:
        files = wt.get_files(path, name)
        files_by_name[name] = files
        conc_ts[name] = {}

        folder_ambient = None
        if parameters_per_folder and ambient_conditions is not None:
            folder_ambient = _read_ambient_for_key(ambient_conditions, name)

        for file in files:
            if (not parameters_per_folder) and ambient_conditions is not None:
                ambient_vals = _read_ambient_for_key(ambient_conditions, file)
            else:
                ambient_vals = folder_ambient

            if ambient_vals is None:
                raise RuntimeError(
                    f"Ambient conditions missing for key '{file if not parameters_per_folder else name}' in '{csv_file}'."
                )

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
            ) = ambient_vals

            pc = wt.PointConcentration.from_file(path + file)
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

            if full_scale == "fs":
                pc.to_full_scale()
            elif full_scale == "nd":
                pc.to_non_dimensional()

            conc_ts[name][file] = pc

    # Optional: uncertainty calculation is kept out of the return for now.
    # It’s still available via `windtunnel.concentration.utils.calculate_uncertainties`.
    _ = (
        calculate_uncertainty,
        split_factor,
        uncertainty_threshold,
        uncertainty_metrics,
        uncertainty_concentration_types,
        include_abs_uncertainty,
        include_pct_uncertainty,
    )

    # Optional: saving (kept compatible with existing scripts)
    if save_all:
        save_ts = save_avg = save_stats = save_combined = True
    if save_combined:
        save_avg = True
        save_stats = True

    if any([save_ts, save_avg, save_stats, save_combined]):
        if output_path is None:
            raise ValueError("output_path is required when saving is enabled.")

        output_path = str(output_path)
        if not output_path.endswith("/"):
            output_path += "/"

        for name in namelist:
            files = files_by_name[name]
            if os_type == "Windows":
                folder = "Point_Data\\" + name[: name.find(".")] + "\\"
                folder_avg = "Point_Data_avg\\" + name[: name.find(".")] + "\\"
                folder_stats = "Point_Data_stats\\" + name[: name.find(".")] + "\\"
            else:
                folder = "Files/Point_Data/" + name[: name.find(".")] + "/"
                folder_avg = "Files/Point_Data_avg/" + name[: name.find(".")] + "/"
                folder_stats = "Files/Point_Data_stats/" + name[: name.find(".")] + "/"

            wt.check_directory(output_path + folder)
            if save_avg:
                wt.check_directory(output_path + folder_avg)
            if save_stats:
                wt.check_directory(output_path + folder_stats)

            for file in files:
                pc = conc_ts[name][file]
                if save_ts:
                    if full_scale == "ms":
                        pc.save2file_ms(file, out_dir=output_path + folder)
                    elif full_scale == "fs":
                        pc.save2file_fs(file, out_dir=output_path + folder)
                    elif full_scale == "nd":
                        pc.save2file_nd(file, out_dir=output_path + folder)
                if save_avg:
                    pc.save2file_avg(file, out_dir=output_path + folder_avg)
                if save_stats:
                    pc.save2file_fullStats(file, out_dir=output_path + folder_stats)

            if save_combined:
                from windtunnel.concentration.utils import combine_to_csv

                file_names = ["_stats_" + file for file in files]
                base_path_local = output_path + f"Files/Point_Data_stats/{name[0:-1]}/"
                combine_to_csv(
                    file_names,
                    base_path_local,
                    file_type="stats",
                    output_filename=output_path + combined_filename,
                )

    return PointConcentrationRunResult(conc_ts=conc_ts, files_by_name=files_by_name, namelist=namelist)


