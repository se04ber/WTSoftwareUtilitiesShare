from __future__ import annotations

import numpy as np

import windtunnel as wt


def test_pointconcentration_relative_coordinates_are_measure_minus_source():
    """
    Guardrail for a subtle but serious bug:
    PointConcentration.ambient_conditions() must compute x/y/z as (measure - source).
    """

    pc = wt.PointConcentration(
        time=np.array([0.0, 1.0]),
        wtref=np.array([1.0, 1.0]),
        slow_FID=np.array([0.0, 0.0]),
        fast_FID=np.array([0.0, 0.0]),
        open_rate=np.array([0.0, 0.0]),
    )

    pc.ambient_conditions(
        x_source=100.0,
        y_source=200.0,
        z_source=23.0,
        x_measure=110.0,
        y_measure=205.0,
        z_measure=5.0,
        pressure=101325.0,
        temperature=20.0,
        calibration_curve=0.0,
        mass_flow_controller=1.0,
        calibration_factor=0.0,
        config_name="cfg",
    )

    assert pc.x == 10.0
    assert pc.y == 5.0
    assert pc.z == -18.0
    np.testing.assert_allclose(pc.distance, np.sqrt(10.0**2 + 5.0**2 + (-18.0) ** 2))


