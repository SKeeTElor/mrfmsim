import numpy as np

def rel_dpol_sat_td_mw(B_offset,B1,Gamma, T2, dB,f_c):
    """Relative change in polarization for time-dependent saturation.
    Ignoring cantilever's field change, considering introduce an external Bz field dB*cos(wm*t),wm>=w_cantilever, so the spin is saturated one or multiple times.

    The result is not a steady-state solution because it ignores T1 relaxation.
    """
    # ignore division error the Exp takes care of the inf, and nan
    # np.seterr(divide="ignore", invalid="ignore")

    atan_omega_i = dB*np.cos(0)+B_offset # T

    atan_omega_f = dB*np.cos(-np.pi)+B_offset # T

    div = np.divide(np.arctan(Gamma*T2*atan_omega_i)-np.arctan(Gamma*T2*atan_omega_f),dB*f_c*2*np.pi)
    r = div*B1**2*Gamma
    # R = np.exp(-r)
    return r

def phase_one_pulse(B_offset, dB):
    """Caluculating phase delay for one pulse saturation due to the resonant condition are met at different times within the cantilever oscillation cycle.
    """
    # Restrict the offset to the sensitive slice [-dB, dB].
    B_offset_phase = np.clip(B_offset, -dB, dB)
    phase_delay = np.arccos(B_offset_phase/dB)
    # R = np.exp(-r)
    return phase_delay

def positiveFeedback_efficiency_halfcycle(rel_dpol, T1,f_c):
    '''Considering use tip-spin interaction to do positive feedback on cantilever, calculating positive feedback efficiency'''

    r = rel_dpol
    R = np.exp(-r)
    w = f_c*2*np.pi
    R1 = np.exp(-np.pi/w/T1)

    posterm = (2-np.divide((1-R),(1-R*R1))*(1+R1)/(1+1/T1**2/w**2)-(1-R1)*np.divide(np.divide(1+R,1-R*R1),(1+r**2/np.pi**2)) )

    velterm = ((1-R)/(1-R*R1)*(1+R1)/(1+1/T1**2/w**2)/(T1*w)+(1-R1)/(1-R*R1)*(1+R)/(1+r**2/np.pi**2)*(r/np.pi))
    
    return np.sqrt(np.square(posterm)+np.square(velterm))*np.exp(1j*np.arctan(np.divide(velterm,posterm))) 

def positiveFeedback_efficiency_onetime(rel_dpol, T1,f_c,phase_delay):
    """Considering use tip-spin interaction to do positive feedback on cantilever, calculating positive feedback efficiency.
    Assuming saturation only happens once during the one cantilever oscillation cycle.

    Parameters:
    rel_dpol: float: relative change rate in polarization, unit:[dimensionless]
    T1: float: spin-lattice relaxation time, unit:[s]
    f_c: float: cantilever's resonance frequency, unit:[Hz]
    Returns:
    float: positive feedback efficiency, unit:[dimensionless]
    """

    r = rel_dpol
    R = np.exp(-r)
    w = f_c*2*np.pi
    R1 = np.exp(-2*np.pi/w/T1)
    
    return 1/(1+1/T1**2/w**2)*(1-R1)*np.divide((1-R),(1-R*R1))*np.sqrt(1+1/T1**2/w**2)*np.exp(-1j*np.arctan(1/(T1*w)))*np.exp(1j*phase_delay)

def force_td(Bzx, positivefeedback_eff, mz_eq, spin_density, grid_voxel,Q,k_c):
    """Calculate force account for the negative sign in the approximation.
    Parameter:
    Bzx: np.ndarray: x gradiant of Bz field, unit:[mT] 
    positivefeedback_eff: float: positive feedback efficiency, unit:[dimensionless],
    mz_eq: np.ndarray: equilibrium magnetization along z-axis, unit:[aN.nm/mT]
    spin_density: np.ndarray: spin density, unit:[spins/m^3]
    grid_voxel: float: volume of a single grid voxel, unit:[m^3]
    Q: float: cantilever's quality factor, unit:[dimensionless]
    k_c: float: cantilever's spring constant, unit:[aN/nm]

    Returns:
    amp_spin: np.array: amplitude of position signal induced by the force, unit:[nm]
    
    """
    return Bzx* positivefeedback_eff*mz_eq*spin_density*grid_voxel*Q/np.pi/k_c

