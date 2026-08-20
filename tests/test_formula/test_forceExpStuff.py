"""Test the relative change in polarization."""

import mrfmsim.formula.forceExpStuff as fes
import numpy as np
from mrfmsim.component import Sample
import pytest


@pytest.fixture
def sample_e():
    """Electron sample."""
    return Sample(spin="e", temperature=0.001, T1=1.0, T2=1.0, spin_density=10.0)


@pytest.fixture
def sample_h():
    """Nucleus sample."""
    return Sample(spin="1H", temperature=4.2, T1=10, T2=5e-6, spin_density=49.0)



def test_rel_dpol_sat_td_mw(sample_e):
    """Test rel_dpol_sat_td when offset is 0.

    When the offset is 0, the exp of negative result should be 0

    """
    ext_B_offset = np.zeros([3])

    rpol = np.exp(-fes.rel_dpol_sat_td_mw(
         ext_B_offset, 1, sample_e.Gamma, sample_e.T2, 1,1000
    ))

    assert np.array_equal(rpol, [0, 0, 0])


def test_rel_dpol_sat_td_mw_symmetry(sample_e):
    """Test rel_dpol_sat_td_mw is symmetric when the B_offset is positive and negative.

    """
    ext_B_offset_a = np.array([-2, -1, 0, 1, 2])
    ext_B_offset_b = np.array([2, 1, 0, -1, -2])

    rpol_a = np.exp(-fes.rel_dpol_sat_td_mw(
         ext_B_offset_a, 1, sample_e.Gamma, sample_e.T2, 1,1000
    ))

    rpol_b = np.exp(-fes.rel_dpol_sat_td_mw(
         ext_B_offset_b, 1, sample_e.Gamma, sample_e.T2, 1,1000
    ))

    assert np.array_equal(rpol_a, rpol_b)


def test_phase_one_pulse():
    """Test phase_one_pulse when offset is 0.
    
    """
    ext_B_offset = np.zeros([3])
    phase = fes.phase_one_pulse(ext_B_offset, 1)

    # R = np.exp(-r)
    assert np.array_equal(phase, [np.pi/2, np.pi/2, np.pi/2])

def test_phase_one_pulse_bigger_than_dB():
    """Test phase_one_pulse when offset is bigger than dB.
    
    """
    ext_B_offset = np.ones([3])
    phase = fes.phase_one_pulse(ext_B_offset, 0.5)

    # R = np.exp(-r)
    assert np.array_equal(phase, [0,0,0])

def test_phase_one_pulse_smaller_than_negdB():
    """Test phase_one_pulse when offset is smaller than -dB.
    
    """
    ext_B_offset = -np.ones([3])
    phase = fes.phase_one_pulse(ext_B_offset, 0.5)

    # R = np.exp(-r)
    assert np.array_equal(phase, [np.pi,np.pi,np.pi])

def test_positiveFeedback_efficiency_halfcycle_lowT1T2 ():
    """Test positiveFeedback_efficiency_halfcycle when T1 and T2 are low, so the force sequence shoule be square waves.
        The efficiency should be 2.
    
    """

    # R = np.exp(-r)
    assert np.isclose(np.abs(fes.positiveFeedback_efficiency_halfcycle(1000, 0.0001,1)),b=2,atol=1e-3)

def test_positiveFeedback_efficiency_halfcycle_highT1T2 ():
    """Test positiveFeedback_efficiency_halfcycle when T1 and T2 are high, so the force sequence shoule be fluctuating in low amplitude.
        The efficiency should be close to 0.
    
    """

    # R = np.exp(-r)
    assert np.isclose(np.abs(fes.positiveFeedback_efficiency_halfcycle(0.000001, 100,1)),b=0,atol=1e-3)



def test_positiveFeedback_efficiency_onetime_lowT1T2 ():
    """Test positiveFeedback_efficiency_onetime when T1 and T2 are low, so the force sequence shoule be a single square waves.
        The efficiency should be 0.
    
    """

    # R = np.exp(-r)
    assert np.isclose(np.abs(fes.positiveFeedback_efficiency_onetime(0.000001, 100,1,0)),b=0,atol=1e-3)

def test_positiveFeedback_efficiency_onetime_highT1T2 ():
    """Test positiveFeedback_efficiency_onetime when T1 and T2 are high, so the force sequence shoule be fluctuating in low amplitude.
        The efficiency should be close to 0.
    
    """

    # R = np.exp(-r)
    assert np.isclose(np.abs(fes.positiveFeedback_efficiency_onetime(10000, 1000000,1,0)),b=0,atol=1e-3)