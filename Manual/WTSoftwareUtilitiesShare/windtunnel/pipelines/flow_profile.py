from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

import windtunnel as wt


@dataclass(frozen=True)
class FlowProfileRunResult:
    time_series: dict
    time_series_eq: dict
    files: list[str]
    namelist: list[str]
    wind_comps: dict
    spectra_data: dict
    heights: list[float]
    u_mean: list[float]
    mean_mag: list[float]
    I_u: list[float]
    I_v: list[float]
    fluxes: list[float]
    lux: list[float]


def run_flow_profile(
    *,
    path: str,
    wtref_path: str,
    namelist: list[str],
    scale_factor: float,
    scale: float,
    Lref: Optional[float] = None,
    data_nd: int,
    mode: int,
    checker: bool,
    save_data: bool,
    outdata_path: Optional[str] = None,
    save_timeseries_txt: bool = False,
    out_dir: Optional[str] = None,
    header_information: Optional[dict] = None,
) -> FlowProfileRunResult:
    """
    Compute-only pipeline for flow profile analysis.

    This mirrors the data preparation step in `WTFlowAnalysis.py` / `example_data_analysis.py`
    and returns the arrays that feed plots.
    """
    if not path.endswith("/"):
        path = path + "/"
    if not wtref_path.endswith("/"):
        wtref_path = wtref_path + "/"
    if Lref is None:
        Lref = 1 / scale
    if save_timeseries_txt and (out_dir is None or header_information is None):
        raise ValueError("out_dir and header_information are required when save_timeseries_txt=True")

    time_series = {}
    time_series.fromkeys(namelist)
    time_series_eq = {}
    time_series_eq.fromkeys(namelist)

    files_for_first: list[str] = []

    for name in namelist:
        files = wt.get_files(path, name)
        if name == namelist[0]:
            files_for_first = files

        time_series[name] = {}
        time_series[name].fromkeys(files)
        time_series_eq[name] = {}
        time_series_eq[name].fromkeys(files)

        for i, file in enumerate(files):
            ts = wt.Timeseries.from_file(path + file)
            ts.get_wind_comps(path + file)
            ts.get_wtref(wtref_path, name, index=i, vscale=scale_factor)

            # Keep legacy behavior: data_nd==0 triggers nondimensionalise() in the scripts.
            if data_nd == 0:
                ts.nondimensionalise()

            ts.adapt_scale(scale, Lref=Lref)
            ts.mask_outliers()

            ts_eq = ts
            ts_eq.calc_equidistant_timesteps()

            ts.index = ts.t_arr
            ts.weighted_component_mean
            ts_eq.weighted_component_mean
            ts.weighted_component_variance
            ts_eq.weighted_component_variance
            ts.mean_magnitude
            ts_eq.mean_magnitude
            ts.mean_direction
            ts_eq.mean_direction

            time_series[name][file] = ts
            time_series_eq[name][file] = ts_eq

            if save_timeseries_txt:
                ts.save2file(file, out_dir=out_dir, header_information=header_information)

    if not files_for_first:
        raise Exception("No Matching File Names Found. Please check namelist and/or path!")

    if checker:
        for name in namelist:
            files = wt.get_files(path, name)
            if mode in {1, 3, 4}:
                for i in range(np.size(files) - 2):
                    if time_series[name][files[i]].x != time_series[name][files[i + 1]].x:
                        raise Exception("Positions do not match! Check data file.")
                    if time_series[name][files[i]].y != time_series[name][files[i + 1]].y:
                        raise Exception("Positions do not match! Check data file.")
            if mode == 2:
                for i in range(np.size(files) - 2):
                    if time_series[name][files[i]].x != time_series[name][files[i + 1]].x:
                        raise Exception("Positions do not match! Check data file.")
                    if time_series[name][files[i]].z != time_series[name][files[i + 1]].z:
                        raise Exception("Positions do not match! Check data file.")
            if mode == 5:
                for i in range(np.size(files) - 2):
                    if time_series[name][files[i]].y != time_series[name][files[i + 1]].y:
                        raise Exception("Positions do not match! Check data file.")
                    if time_series[name][files[i]].z != time_series[name][files[i + 1]].z:
                        raise Exception("Positions do not match! Check data file.")

    wind_comps = {}
    wind_comps.fromkeys(namelist)
    turb_data = {}
    turb_data.fromkeys(namelist)
    lux_data = {}
    lux_data.fromkeys(namelist)
    spectra_data = {}
    spectra_data.fromkeys(namelist)

    for name in namelist:
        files = wt.get_files(path, name)
        wind_comps[name] = {}
        wind_comps[name].fromkeys(files)
        turb_data[name] = {}
        turb_data[name].fromkeys(files)
        lux_data[name] = {}
        lux_data[name].fromkeys(files)
        spectra_data[name] = {}
        spectra_data[name].fromkeys(files)

        for file in files:
            wind_comps[name][file] = time_series[name][file].wind_comp1, time_series[name][file].wind_comp2
            turb_data[name][file] = wt.calc_turb_data(
                time_series[name][file].u.dropna(), time_series[name][file].v.dropna()
            )
            dt = time_series[name][file].t_eq[1] - time_series[name][file].t_eq[0]
            if data_nd == 0:
                lux_data[name][file] = wt.calc_lux_data(
                    dt, (time_series[name][file].u_eq.dropna().values * time_series[name][file].wtref)
                )
            else:
                lux_data[name][file] = wt.calc_lux_data(dt, (time_series[name][file].u_eq.dropna().values))

    # Prepare arrays for plotting (first namelist entry, matching current scripts)
    name0 = namelist[0]
    files0 = wt.get_files(path, name0)
    heights: list[float] = []
    mean_mag: list[float] = []
    u_mean: list[float] = []
    I_u: list[float] = []
    I_v: list[float] = []
    fluxes: list[float] = []
    lux: list[float] = []
    x: list[float] = []
    y: list[float] = []
    wtref: list[float] = []

    for file in files0:
        mean_mag.append(time_series[name0][file].mean_magnitude)
        u_mean.append(float(np.mean(time_series[name0][file].u)))
        wtref.append(float(time_series[name0][file].wtref))
        if mode != 4:
            x.append(float(time_series[name0][file].x))
            y.append(float(time_series[name0][file].y))
            heights.append(float(time_series[name0][file].z))
            I_u.append(float(turb_data[name0][file][0]))
            I_v.append(float(turb_data[name0][file][1]))
            fluxes.append(float(turb_data[name0][file][2]))
            lux.append(float(lux_data[name0][file]))

    if save_data:
        if outdata_path is None:
            raise ValueError("outdata_path must be set when save_data=True")
        wt.check_directory(outdata_path)
        outfile = outdata_path + name0 + ".npz"
        np.savez(
            outfile,
            x=x,
            y=y,
            heights=heights,
            mean_mag=mean_mag,
            u_mean=u_mean,
            I_u=I_u,
            I_v=I_v,
            fluxes=fluxes,
            lux=lux,
            wtref=wtref,
        )

    return FlowProfileRunResult(
        time_series=time_series,
        time_series_eq=time_series_eq,
        files=files0,
        namelist=namelist,
        wind_comps=wind_comps,
        spectra_data=spectra_data,
        heights=heights,
        u_mean=u_mean,
        mean_mag=mean_mag,
        I_u=I_u,
        I_v=I_v,
        fluxes=fluxes,
        lux=lux,
    )


