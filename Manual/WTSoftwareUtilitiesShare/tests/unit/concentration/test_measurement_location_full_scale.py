from __future__ import annotations

from pathlib import Path

import numpy as np

from windtunnel.concentration.utils import load_avg_file


def test_avg_parser_uses_measurement_height_for_full_scale_coordinates(tmp_path: Path):
    """
    Reported bug:
    If the parameter file contains z_measure=5mm (probe height), some outputs showed
    Z_fs=-7.2m when z_source=23mm and scale=400 (i.e. z_rel=-18mm).

    For reporting "measurement location", full scale height should be based on
    z_measure, not (z_measure - z_source).
    """

    scale = 400.0

    # Model-scale positions [mm]
    z_source = 23.0
    z_measure = 5.0
    z_rel = z_measure - z_source  # -18mm -> -7.2m in full scale if interpreted as absolute

    # Minimal avg file content with the same header patterns that load_avg_file parses.
    content = "\n".join(
        [
            "# General concentration measurement data:",
            "#",
            f"# geometric scale: 1:{scale}",
            "# Variables: "
            f"x (measurement relativ to source): 0 [mm], "
            f"y (measurement relativ to source): 0 [mm], "
            f"z (measurement relativ to source): {z_rel} [mm], "
            "x_source: 0 [mm], y_source: 0 [mm], z_source: 23 [mm], "
            "x_measure: 0 [mm], y_measure: 0 [mm], z_measure: 5 [mm], "
            "ambient temperature: 20.0 [°C], ambient pressure: 101325.000000 [Pa], "
            "mass flow rate 1.000000 [kg/s], reference length (model): 0.002500 [m], "
            "Tracer gas: X, mol. weight tracer: 0.029000 [mol/kg], gas factor: 1.000000, "
            "calibartion curve: 0.000000, wtref: 1.000000 [m/s], "
            "full scale wtref: 3.0000[m/s], full scale flow rate: 0.5000 [m^3/s]",
            '# "time [ms]" "c_star [-]" "net_concentration [ppmV]" "full_scale_concentration [ppmV]"',
            # data line (3+ columns required by parser)
            "0.1 0.2 0.3",
            "",
        ]
    )

    p = tmp_path / "_avg_dummy.txt.ts#0"
    p.write_text(content, encoding="utf-8")

    meta = load_avg_file(str(p))
    assert meta is not None

    # The key check: measurement height => 5mm * 400 / 1000 = 2m
    assert np.isclose(meta["z_fs"], z_measure / 1000.0 * scale)

    # And keep the relative-to-source value for debugging (optional)
    assert np.isclose(meta["z_rel_fs"], z_rel / 1000.0 * scale)


