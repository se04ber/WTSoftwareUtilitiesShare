from __future__ import annotations

import numpy as np

import windtunnel as wt


def _make_pc():
    # Minimal PointConcentration object for flow-rate/c* computations.
    n = 10
    time = np.linspace(0, 1, n)
    wtref = np.ones(n)
    slow = np.zeros(n)
    fast = np.zeros(n)
    open_rate = np.full(n, 0.3)  # [%/10] in legacy assumptions
    pc = wt.PointConcentration(time, wtref, slow, fast, open_rate)

    pc.temperature = 23.0
    pc.pressure = 101768.0
    pc.standard_temp = 0.0
    pc.standard_pressure = 101325.0
    pc.full_scale_temp = 20.0
    pc.gas_factor = 1.0
    pc.mass_flow_controller = 0.300
    pc.calibration_curve = None
    pc.calibration_factor = None

    pc.convert_temperature()
    return pc


def test_calc_model_mass_flow_rate_using_max_flow_rate_no_calibration():
    pc = _make_pc()
    q = pc.calc_model_mass_flow_rate(usingMaxFlowRate="True", applyCalibration="False")

    expected = (
        (np.mean(pc.open_rate) * 10 * pc.gas_factor * pc.mass_flow_controller)
        * pc.temperature_K
        * pc.standard_pressure
        / (pc.pressure * pc.standard_temp_K)
    )
    np.testing.assert_allclose(q, expected)


def test_calc_model_mass_flow_rate_using_max_flow_rate_with_calibration():
    pc = _make_pc()
    pc.calibration_curve = 2.0
    pc.calibration_factor = 1.0
    q = pc.calc_model_mass_flow_rate(usingMaxFlowRate="True", applyCalibration="True")

    expected = (
        (
            np.mean(pc.open_rate)
            * 10
            * pc.gas_factor
            * pc.mass_flow_controller
            * pc.calibration_curve
            + pc.calibration_factor
        )
        * pc.temperature_K
        * pc.standard_pressure
        / (pc.pressure * pc.standard_temp_K)
    )
    np.testing.assert_allclose(q, expected)


def test_calc_model_mass_flow_rate_legacy_mode_with_calibration():
    pc = _make_pc()
    pc.calibration_curve = 2.0
    pc.calibration_factor = 1.0
    q = pc.calc_model_mass_flow_rate(usingMaxFlowRate="False", applyCalibration="True")

    expected = (
        pc.gas_factor
        * (np.mean(pc.open_rate) * 10 * pc.calibration_curve + pc.calibration_factor)
        * pc.temperature_K
        * pc.standard_pressure
        / (pc.pressure * pc.standard_temp_K)
    )
    np.testing.assert_allclose(q, expected)


