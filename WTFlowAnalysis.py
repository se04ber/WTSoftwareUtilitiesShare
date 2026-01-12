#!/usr/bin/env python3
# Wind tunnel flow analysis script
# Configuration: lines 8-120
# Add new plots: see ### ADD NEW PLOTS HERE ### in run_analysis()

import os
import sys
import numpy as np
import logging
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath("/home/sabrina/Desktop/Schreibtisch/Arbeit_2025/WTSoftwareUtilitiesShare"))
import windtunnel as wt
from windtunnel.flow import *

logger = logging.getLogger()

# Paths
path = "/home/sabrina/Desktop/Schreibtisch/Arbeit_2025/FreeCAD/UBA/vertikales Profil/Zeitserien/"
wtref_path = "/home/sabrina/Desktop/Schreibtisch/Arbeit_2025/FreeCAD/UBA/vertikales Profil/wtref/"
namelist = ['UBA_BL_Rep2025_UW_01']
txt_path = "/home/sabrina/Desktop/Schreibtisch/Arbeit_2025/FreeCAD/UBA/vertikales Profil/Zeitserien/"
ref_path = "/home/sabrina/48/USA/"
plot_path = '/home/sabrina/Desktop/Schreibtisch/Arbeit_2025/windtunnel_software/TestPlots/'
out_dir = '/home/sabrina/Desktop/Schreibtisch/Arbeit_2025/'
# Header information
header_information = {
    "[main project name]": 'BFS',
    "[sub project name]": 'Boundary layer',
    "[wind tunnel facility]": 'WOTAN',
    "[model scale]": 200,
    "[Zref - model scale[mm]]": 200,
    "[wind direction [°]]": 0,
    "[Lref - model scale[m]]": 1/200,
    "[number of columns]": 3,
    "[directly measured components]": 'UW',
    "[confidence intervall (U,V,W)/Uref[-]]": 0.03
}
scale = header_information["[model scale]"]
scale_factor = 0.6568
data_nd = 0  # 0=dimensional, 1=non-dimensional
file_type = 'png'
mode = 1  # 1=vertical, 2=lateral, 3=convergence, 4=Reynolds, 5=longitudinal
checker = True
plot = True
scatter = True
save_data = True
new_VDI_ref = True
# Plot settings
z_axis_lim = (None, None)
dpi = 300
shown_z_value = (None, None)
# Error values
u_err = 0
v_err = 0
u_v_err = u_err
I_u_err = 0.1
I_v_err = 0.1
flux_err = 0.1
lux_err = 0.1
# Plotting flags
plot_convergence = True
plot_scatter_hist = True
plot_spectra = True
plot_wavelet = True
plot_turb_int_u = True
plot_turb_int_v = True
plot_fluxes = True
plot_lux = True
plot_wind_profile_alpha = True
plot_wind_profile_z0 = True
plot_combined_profiles = True
plot_reynolds = True
# VDI parameters
zref = 10
d0 = 0
Uref = 5
Kappa = 0.4
# Fitting options
fitAlpha = True
fitZ0 = False
alpha_value = 0.17
z_ref = 50
# Plot customization
component2Plot_u = "I_u"
component2Plot_v = "I_w"
Labels_turb = ["z0=5mm", "z0=0.1m", "z0=0.5m", "z0=2m", "Windkanalmessungen"]
xLabel_turb_u = r"$I_{u} $ [-]"
yLabel_turb = r"$Z_{fs} [m] $ [-]"
xLabel_turb_v = r"$I_{w} $ [-]"
show_Legend = True
LegendSize = 0.001
figSize = (5, 10)
figSize_lux = [13, 9]
# Spectra settings
comps2Plot = ["u"]
xLabel_spectra = r"$f\cdot z\cdot U^{-1}$ [-]"
yLabel_spectra = r"$f\cdot S_{uu}\cdot (\sigma_u\sigma_u)^{-1}$ [-]"
Labels_spectra = ["Windkanalmessung u", "Windkanalmessung u", "max u", "min u"]
showLegend_spectra = True
xAchse_spectra = (10**-3., 2.*10.**1.)
yAchse_spectra = (10**-4, 10.**0.)
refRangeVis = "maxMinLines"
# Flux settings
invertxAchse = True
xLabel_flux = r' - u' + '\'' + 'w' + '\' $\cdot U_{ref}^{-2}$ $[-]$'
yLabel_flux = r' $Z_{fs}$ [m]'
xAchse_flux = (0.0, 0.0070)
yAchse_flux = (0, 70)
showLegend_flux = False
# Lux settings
xLabel_lux = r"$L_{ux} $ [m]"
yLabel_lux = r"$Z_{fs} $ [m]"
xAchse_lux = [10, 1000]
yAchse_lux = [0, 300]
showLegend_lux = True
Labels_lux = ["Windkanalmessungen", "z0=10m (Theorie)", "z0=1m (Theorie)", "z0=0.1m (Theorie)",
              "z0=0.01 (Theorie)", "Feldmessung (niedrige)", "Feldmessung (hohe)"]
# Wind profile settings
components_wind = ["u"]
err_wind = 0.04
Labels_wind = ["u"]
xLabel_wind = r"$U/U_{ref} $ [-]"
yLabel_wind = r"$Z_{fs} $ [m]"
xAchse_wind = [0.0, 1.0]
yAchse_wind = [0, 40]
showLegend_wind = False

def main():
    wt.check_directory(plot_path)
    wt.check_directory(txt_path)
    if mode == 2:
        outdata_path = '/path/to/your/outdata_lat/'
    else:
        outdata_path = plot_path
    time_series = {}
    time_series.fromkeys(namelist)
    time_series_eq = {}
    time_series_eq.fromkeys(namelist)
    for name in namelist:
        files = wt.get_files(path, name)
        time_series[name] = {}
        time_series[name].fromkeys(files)
        time_series_eq[name] = {}
        time_series_eq[name].fromkeys(files)
    # Load and process data
    for name in namelist:
        for i, file in enumerate(files):
            ts = wt.Timeseries.from_file(path + file)
            ts.get_wind_comps(path + file)
            ts.get_wtref(wtref_path, name, index=i, vscale=scale_factor)
            print(file)
            print(ts.u)
            if data_nd == 0:
                print('Warning: Assuming that data is dimensional. If using non-dimensional input data, set variable data_nd to 1')
                ts.nondimensionalise()
                print("test")
                print(ts.u)
            ts.adapt_scale(scale, Lref=header_information["[Lref - model scale[m]]"])
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
            ts.save2file(file, out_dir=out_dir, header_information=header_information)
            time_series[name][file] = ts
            time_series_eq[name][file] = ts_eq
    if files == []:
        raise Exception('No Matching File Names Found. Please check namelist and/or path!')
    # Check positions
    if checker:
        for name in namelist:
            files = wt.get_files(path, name)
            if mode == 1 or mode == 3 or mode == 4:
                for i in range(np.size(files) - 2):
                    if time_series[name][files[i]].x != time_series[name][files[i + 1]].x:
                        raise Exception('Positions do not match! Check data file.')
                    if time_series[name][files[i]].y != time_series[name][files[i + 1]].y:
                        raise Exception('Positions do not match! Check data file.')
            if mode == 2:
                for i in range(np.size(files) - 2):
                    if time_series[name][files[i]].x != time_series[name][files[i + 1]].x:
                        raise Exception('Positions do not match! Check data file.')
                    if time_series[name][files[i]].z != time_series[name][files[i + 1]].z:
                        raise Exception('Positions do not match! Check data file.')
            if mode == 5:
                for i in range(np.size(files) - 2):
                    if time_series[name][files[i]].y != time_series[name][files[i + 1]].y:
                        raise Exception('Positions do not match! Check data file.')
                    if time_series[name][files[i]].z != time_series[name][files[i + 1]].z:
                        raise Exception('Positions do not match! Check data file.')
    # Calculate statistics
    wind_comps = {}
    wind_comps.fromkeys(namelist)
    wind_stats = {}
    wind_stats.fromkeys(namelist)
    turb_data = {}
    turb_data.fromkeys(namelist)
    lux_data = {}
    lux_data.fromkeys(namelist)
    spectra_data = {}
    spectra_data.fromkeys(namelist)
    for name in namelist:
        n_outliers_u = 0
        n_outliers_v = 0
        files = wt.get_files(path, name)
        wind_comps[name] = {}
        wind_comps[name].fromkeys(files)
        if mode != 3:
            wind_stats[name] = {}
            wind_stats[name].fromkeys(files)
            turb_data[name] = {}
            turb_data[name].fromkeys(files)
            lux_data[name] = {}
            lux_data[name].fromkeys(files)
            spectra_data[name] = {}
            spectra_data[name].fromkeys(files)
            for file in files:
                n_outliers_u = n_outliers_u + time_series[name][file].n_outliers_u
                n_outliers_v = n_outliers_v + time_series[name][file].n_outliers_v
                wind_comps[name][file] = time_series[name][file].wind_comp1, time_series[name][file].wind_comp2
    for name in namelist:
        files = wt.get_files(path, name)
        for file in files:
            wind_stats[name][file] = wt.calc_wind_stats_wght(time_series[name][file].t_transit,
                                                              time_series[name][file].u,
                                                              time_series[name][file].v)
            turb_data[name][file] = wt.calc_turb_data(time_series[name][file].u.dropna(),
                                                       time_series[name][file].v.dropna())
            dt = time_series[name][file].t_eq[1] - time_series[name][file].t_eq[0]
            if data_nd == 0:
                lux_data[name][file] = wt.calc_lux_data(dt, (time_series[name][file].u_eq.dropna().values *
                                                              time_series[name][file].wtref))
            if data_nd == 1:
                lux_data[name][file] = wt.calc_lux_data(dt, (time_series[name][file].u_eq.dropna().values))
    # Prepare arrays for plotting
    for name in namelist:
        x = []
        y = []
        heights = []
        mean_mag = []
        u_mean = []
        u_mean_wght = []
        u_std = []
        u_std_wght = []
        v_mean = []
        v_mean_wght = []
        v_std = []
        v_std_wght = []
        I_u = []
        I_v = []
        fluxes = []
        lux = []
        wdir = []
        wtref = []
        for file in files:
            mean_mag.append(time_series[name][file].mean_magnitude)
            u_mean.append(np.mean(time_series[name][file].u))
            u_std.append(np.std(time_series[name][file].u))
            wtref.append(time_series[name][file].wtref)
            if mode != 4:
                x.append((time_series[name][file].x))
                y.append((time_series[name][file].y))
                heights.append((time_series[name][file].z))
                u_mean_wght.append(time_series[name][file].weighted_component_mean[0])
                u_std_wght.append(np.sqrt(time_series[name][file].weighted_component_variance[0]))
                v_mean.append(np.mean(time_series[name][file].v))
                v_mean_wght.append(time_series[name][file].weighted_component_mean[1])
                v_std.append(np.std(time_series[name][file].v))
                v_std_wght.append(np.sqrt(time_series[name][file].weighted_component_variance[1]))
                wdir.append(time_series[name][file].mean_direction)
                I_u.append(turb_data[name][file][0])
                I_v.append(turb_data[name][file][1])
                fluxes.append(turb_data[name][file][2])
                lux.append(lux_data[name][file])
        if save_data:
            wt.check_directory(outdata_path)
            outfile = outdata_path + name + '.npz'
            np.savez(outfile, x=x, y=y, heights=heights, mean_mag=mean_mag, u_mean=u_mean, v_mean=v_mean,
                     I_u=I_u, I_v=I_v, fluxes=fluxes, lux=lux, wtref=wtref)
    return time_series, time_series_eq, files, namelist, wind_comps, spectra_data, heights, u_mean, mean_mag, I_u, I_v, fluxes, lux

def run_analysis(time_series, time_series_eq, files, namelist, wind_comps, spectra_data, heights, u_mean, mean_mag, I_u, I_v, fluxes, lux):
    for name in namelist:
        for file in files:
            ts = time_series[name][file]
            # Convergence plots
            if plot_convergence and (mode == 1 or mode == 2 or mode == 5) and scatter:
                u_convergence = wt.convergence_test(time_series_eq[name][file].u_eq)
                v_convergence = wt.convergence_test(time_series_eq[name][file].v_eq)
                yLabel = time_series[name][file].wind_comp1 + '_mean'
                fig1, ax1 = plt.subplots(1)
                wt.plot_convergence_test(u_convergence, ylabel=yLabel, Label=file[:-4], ax=ax1)
                plt.tight_layout()
                fig1.savefig(plot_path + 'convergence_' + time_series[name][file].wind_comp1 + '_mean' + '.' + file_type, dpi=1000, bbox_inches='tight')
                plt.show()
                fig2, ax2 = plt.subplots(1)
                wt.plot_convergence_test(v_convergence, ylabel=time_series[name][file].wind_comp2 + '_mean', ax=ax2)
                plt.tight_layout()
                fig2.savefig(plot_path + 'convergence_' + time_series[name][file].wind_comp2 + '_mean' + '.' + file_type, dpi=1000, bbox_inches='tight')
                plt.show()
            # Scatter and histograms
            if plot_scatter_hist and (mode == 1 or mode == 2 or mode == 5) and scatter:
                plt.figure(files.index(file) + 100)
                wt.plots.plot_scatter(time_series[name][file].u, time_series[name][file].v)
                plt.savefig(plot_path + 'scatter_' + file[:-4] + '.' + file_type)
                plt.show()
                plt.figure(files.index(file) + 200)
                wt.plots.plot_hist(time_series[name][file].u)
                plt.savefig(plot_path + 'hist_u_' + file[:-4] + '.' + file_type)
                plt.show()
                plt.figure(files.index(file) + 300)
                wt.plots.plot_hist(time_series[name][file].v)
                plt.savefig(plot_path + 'hist_v_' + file[:-4] + '.' + file_type)
                plt.show()
            # Spectra
            if plot_spectra and (mode == 1 or mode == 2 or mode == 5) and scatter:
                spectra_data[name][file] = wt.calc_spectra(
                    time_series_eq[name][file].u_eq.dropna() * time_series_eq[name][file].wtref,
                    time_series_eq[name][file].v_eq.dropna() * time_series_eq[name][file].wtref,
                    time_series_eq[name][file].t_eq[~np.isnan(time_series_eq[name][file].t_eq)],
                    time_series_eq[name][file].z)
                plt.figure(files.index(file) + 400)
                wt.plots.plot_spectra(spectra_data[name][file][0], spectra_data[name][file][1],
                                      spectra_data[name][file][2], spectra_data[name][file][4],
                                      spectra_data[name][file][5], wind_comps[name][file],
                                      time_series[name][file].z, ref_path=ref_path, xLabel=xLabel_spectra,
                                      yLabel=yLabel_spectra, Labels=Labels_spectra, showLegend=showLegend_spectra,
                                      xAchse=xAchse_spectra, yAchse=yAchse_spectra, refRangeVis=refRangeVis,
                                      comps2Plot=comps2Plot)
                plt.savefig(plot_path + 'spectra_' + file[:-4] + '.' + file_type)
                plt.show()
            # Wavelet
            if plot_wavelet:
                print('\n Start wavelet analysis for {}'.format(file))
                wavelet, scale = wt.calc_wavelet_transform(time_series[name][file].u_eq,
                                                            time_series[name][file].t_eq,
                                                            wavelet='morlet')
                y_val = time_series[name][file].z
                f_scale = y_val / (scale * np.mean(time_series[name][file].u_eq))
                plt.figure(55)
                plt.style.use('classic')
                wt.plots.plot_wavelet_transform(wavelet, scale, time_series[name][file].u_eq,
                                                time_series[name][file].t_eq, y_val)
                plt.savefig(plot_path + 'wavelet_analysis_' + file + '.' + file_type, bbox_inches='tight', dpi=300)
                print(' Saved plot to:' + plot_path + 'wavelet_analysis_' + file + '.' + file_type)
                plt.show()
                print(' Finished wavelet analysis for {}'.format(file))
    for name in namelist:
        # Turbulence intensity u
        if plot_turb_int_u and mode == 1:
            plt.style.use('default')
            fig = plt.figure(2, dpi=dpi)
            wt.plots.plot_turb_int(I_u, heights, yerr=I_u_err, component=component2Plot_u, ref_path=ref_path,
                                   new_ref=new_VDI_ref, cut=shown_z_value, Labels=Labels_turb, xLabel=xLabel_turb_u,
                                   yLabel=yLabel_turb)
            if (z_axis_lim[0] != None) and (z_axis_lim[1] != None):
                plt.ylim(z_axis_lim)
            plt.tight_layout()
            plt.savefig(plot_path + 'I_' + ts.wind_comp1 + '_' + name + '.' + file_type)
            plt.show()
        # Turbulence intensity v
        if plot_turb_int_v and mode == 1:
            plt.style.use('default')
            plt.figure(3, dpi=dpi)
            wt.plots.plot_turb_int(I_v, heights, yerr=I_v_err, component=component2Plot_v, ref_path=ref_path,
                                   new_ref=new_VDI_ref, cut=shown_z_value, Labels=Labels_turb, xLabel=xLabel_turb_v,
                                   yLabel=yLabel_turb)
            if (z_axis_lim[0] != None) and (z_axis_lim[1] != None):
                plt.ylim(z_axis_lim)
            plt.tight_layout()
            plt.savefig(plot_path + 'I_' + ts.wind_comp2 + '_' + name + '.' + file_type)
            plt.show()
        # Fluxes
        if plot_fluxes and mode == 1:
            if invertxAchse:
                fluxes = fluxes * np.ones_like(fluxes) * (-1)
            plt.style.use('default')
            plt.figure(4, dpi=dpi)
            wt.plots.plot_fluxes(fluxes, heights, yerr=flux_err, component='w', cut=shown_z_value,
                                 xLabel=xLabel_flux, yLabel=yLabel_flux, xAchse=xAchse_flux, yAchse=yAchse_flux,
                                 showLegend=showLegend_flux)
            plt.tight_layout()
            plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
            plt.show()
        # Lux
        if plot_lux and mode == 1:
            plt.style.use('default')
            plt.figure(6, dpi=dpi)
            wt.plots.plot_lux(lux, heights, err=lux_err, component='w', ref_path=ref_path, new_ref=new_VDI_ref,
                              xLabel=xLabel_lux, yLabel=yLabel_lux, Labels=Labels_lux, xAchse=xAchse_lux,
                              yAchse=yAchse_lux, showLegend=showLegend_lux, figSize=figSize_lux)
            plt.tight_layout()
            plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)
            plt.show()
        # Wind profile with alpha
        if plot_wind_profile_alpha and mode == 1 and plot:
            if fitAlpha:
                alpha, ref = calc_alpha(u_mean, heights)
                alpha = alpha_value
                x = np.arange(0, len(u_mean), 1)
                profile2 = [Uref * ((heights[i] - d0) / (z_ref - d0)) ** alpha for i in x] - np.ones_like(ref) * (0.26)
            plt.style.use('default')
            plt.figure(0, dpi=dpi)
            ret = wt.plots.plot_winddata(heights, u_mean=u_mean, mean_magnitude=mean_mag, components=components_wind,
                                         yerr=err_wind, xLabel=xLabel_wind, yLabel=yLabel_wind, Labels=Labels_wind,
                                         xAchse=xAchse_wind, yAchse=yAchse_wind, showLegend=showLegend_wind)
            plt.plot(profile2, heights, "k--", color="black")
            plt.legend([f"alpha: {alpha}", Labels_wind], loc="upper left")
            plt.show()
        # Wind profile with z0
        if plot_wind_profile_z0 and mode == 1 and plot:
            z_fit = np.linspace(10, 100, 100)
            if fitAlpha:
                alpha, ref = calc_alpha(u_mean, heights)
                alpha = alpha_value
                fit_alpha = [((z_fit[i] - 0) / (z_ref - 0)) ** alpha for i in range(100)]
                plt.style.use('default')
                plt.figure(0, dpi=dpi)
                ret = wt.plots.plot_winddata(heights, u_mean=u_mean, mean_magnitude=mean_mag,
                                             components=components_wind, yerr=err_wind, xLabel=xLabel_wind,
                                             yLabel=yLabel_wind, Labels=Labels_wind, xAchse=xAchse_wind,
                                             yAchse=yAchse_wind, showLegend=showLegend_wind)
                plt.plot(fit_alpha, z_fit, "k-", color="darkgrey", linewidth=0.5)
                plt.legend([f"Profilexponent α = {alpha:.3f} ± {err_wind}", Labels_wind[0]], loc="upper left")
                plt.show()
            elif fitZ0:
                z0 = calc_z0(u_mean, heights, d0=0., sfc_height=100.)
                z0 = z0[0]
                u_star = (Kappa) / np.log((zref - d0) / z0)
                fit_z0 = (u_star / Kappa) * np.log((z_fit - d0) / z0)
                plt.figure(0, dpi=dpi)
                ret = wt.plots.plot_winddata(heights, u_mean=u_mean, mean_magnitude=mean_mag,
                                             components=components_wind, yerr=err_wind, xLabel=xLabel_wind,
                                             yLabel=yLabel_wind, Labels=Labels_wind, xAchse=xAchse_wind,
                                             yAchse=yAchse_wind, showLegend=showLegend_wind)
                plt.plot(fit_z0, z_fit, "k-", color="darkgrey", linewidth=0.5)
                plt.legend([f"z0: {z0:.3f}", Labels_wind[0]], loc="upper left")
                plt.show()
        # Combined profiles
        if plot_combined_profiles and mode == 1 and plot:
            alpha = alpha_value
            z0 = 0.1410726917895896
            zref = 50
            d0 = 0
            kappa = 0.4
            z_fit = np.linspace(10, 70, 100)
            u_power = (z_fit / zref) ** alpha
            u_star = (kappa) / np.log((zref - d0) / z0)
            u_log = (u_star / kappa) * np.log((z_fit - d0) / z0)
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
            ax1.errorbar(u_mean, heights, xerr=0.04, c="darkred", fmt="o", label='Measurements', zorder=5,
                         ecolor="black", capsize=3.0)
            ax1.plot(u_power, z_fit, "k-", color="darkgrey", linewidth=0.5, label=f'Power Law (α={alpha:.3f})')
            ax1.set_xlabel('Velocity [-]')
            ax1.set_ylabel('Height (m)')
            ax1.set_title('Linear Scale')
            ax1.grid(True)
            ax1.legend()
            ax2.errorbar(u_mean, heights, xerr=0.04, fmt="o", c='darkred', label='Measurements', zorder=5,
                         ecolor="black", capsize=3.0)
            ax2.plot(u_log, z_fit, "k-", linewidth=0.5, label=f'Log Law (z₀={z0:.3f}m)', color="darkgrey")
            ax2.set_xlabel('Velocity (m/s)')
            ax2.set_ylabel('Height (m) - Log Scale')
            ax2.set_yscale('log')
            ax2.set_title('Logarithmic Scale')
            ax2.grid(True)
            ax2.legend()
            plt.tight_layout()
            plt.suptitle('Wind Profile Analysis', y=1.02, fontsize=14)
            plt.show()
        # Reynolds independence
        if plot_reynolds and mode == 4:
            print("Reynolds independence test")
            plt.style.use('default')
            plt.figure(0)
            wt.plots.plot_Re_independence(mean_mag, wtref, yerr=u_err, ymin=0, ymax=1)
            plt.tight_layout()
            plt.savefig(plot_path + 'Re_u' + name + '.' + file_type)
            plt.legend(loc="lower right")
            plt.show()
            break
    ### ADD NEW PLOTS HERE ###
    # Example: if plot_my_new_analysis:
    #     my_analysis_function(time_series, heights, param1=value1)
    print("Analysis complete")

if __name__ == "__main__":
    time_series, time_series_eq, files, namelist, wind_comps, spectra_data, heights, u_mean, mean_mag, I_u, I_v, fluxes, lux = main()
    run_analysis(time_series, time_series_eq, files, namelist, wind_comps, spectra_data, heights, u_mean, mean_mag, I_u, I_v, fluxes, lux)
