#suppressing a warning - UserWarning: Wswiglal-redir-stdio
import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
import lal
import lal as _lal

# Import lalsuite packages
import lalsimulation, copy
from lalsimulation import gwsignal

# Import astropy and gwpy 
import astropy.units as u

# pycbc has a number of very useful functions and wrappers 
from pycbc import waveform

from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.interpolate import interp1d

import sys

sys.path.insert(0, "../")
from src.regular_precession import *
from src.systems_lib import *
from src.mismatch_n_SNR import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams["figure.dpi"] = 200

import precession

error_handler = np.seterr(invalid="raise")

solar_mass = 4.92624076 * 1e-6              # [solar_mass] = sec
giga_parsec = 1.02927125 * 1e17             # [giga_parsec] = sec
year = 31557600                             # [year] = sec

rp_params = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 2,
    "phi_J": np.pi / 2,
    "mcz": 20 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.25,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 4.0,             # dimensionless
    "omega_tilde": 2.0,             # dimensionless
    "gamma_P": 0.0,
}

rp_params = redshifted_new_params(0.3, rp_params)

# Let's define some convenient astropy units
Msun          = u.solMass
dimensionless = u.dimensionless_unscaled
seconds       = u.s
Mpc           = u.Mpc
radians       = u.rad
Hz            = u.Hz

f_min = 20
f_cut = Regular_precession(rp_params).get_f_cut()
delta_f = 0.05
f_range = np.arange(f_min, f_cut, delta_f)

theta_1 = np.pi/4
theta_2 = np.pi/3
phi_1 = np.pi/6
phi_2 = np.pi/4

chi_1 = 1.0
chi_2 = 1.0


parameters_RP_to_IMRPhenomXPHM = {
    'chi1': chi_1,
    'chi2': chi_2,
    'chirp_mass': rp_params['mcz']/solar_mass,
    'eta': rp_params['eta'],
    'theta_S': rp_params['theta_S'],
    'phi_S': rp_params['phi_S'],
    'theta_LJ': Regular_precession(rp_params).get_theta_LJ(f_range[0]),
    'phi_LJ': Regular_precession(rp_params).get_phi_LJ(f_range[0]),
    'cos_i_JN': Regular_precession(rp_params).precession_angles()[0],
    'sin_i_JN': Regular_precession(rp_params).precession_angles()[1],
    'cos_o_XH': Regular_precession(rp_params).precession_angles()[2],
    'sin_o_XH': Regular_precession(rp_params).precession_angles()[3],
    'theta_1': theta_1,
    'theta_2': theta_2,
    'phi_1': phi_1,
    'phi_2': phi_2,
    'gamma_P': rp_params['gamma_P']
}

def get_m1_m2_IMRPhenomXPHM(params):
    """
    Get the component masses m1 and m2 from the parameters for the IMRPhenomXPHM waveform model.

    Parameters:
    params (dict): A dictionary containing the parameters of the system.

    Returns:
    tuple: A tuple containing the component masses m1 and m2.
    """
    # Extract the relevant parameters
    chirp_mass = params['chirp_mass']
    eta = params['eta']

    # Calculate m1 and m2 from the chirp mass and symmetric mass ratio
    q = (1 - 2 * eta - np.sqrt(1 - 4 * eta)) / (2 * eta)
    total_mass = chirp_mass / (eta**(3/5))
    
    m1 = total_mass / (1 + q)
    m2 = total_mass - m1

    return m1, m2

def get_Fp_Fc_IMRPhenomXPHM(params, psi):
    """
    Get the polarization sensitivity factors Fp and Fc for the IMRPhenomXPHM waveform model.
    
    Parameters:
    params (dict): A dictionary containing the parameters of the system.
    
    Returns:
    tuple: A tuple containing the polarization sensitivity factors Fp and Fc.
    """
    # Extract the relevant parameters
    theta_S = params['theta_S']
    phi_S = params['phi_S']
    cos_i_JN = params['cos_i_JN']
    sin_i_JN = params['sin_i_JN']
    cos_o_XH = params['cos_o_XH']
    sin_o_XH = params['sin_o_XH']
    
    #X_dot_l = np.cos(i_JN) * np.sin(Omega_XH) * np.sin(phi_S) - np.cos(i_JN) * np.cos(Omega_XH) * np.cos(theta_S) * np.cos(phi_S) + np.sin(i_JN) * np.sin(theta_S) * np.cos(phi_S)
    X_dot_l = cos_i_JN * sin_o_XH * np.sin(phi_S) - cos_i_JN * cos_o_XH * np.cos(theta_S) * np.cos(phi_S) + sin_i_JN * np.sin(theta_S) * np.cos(phi_S)
    #Y_dot_l = np.cos(Omega_XH) * np.sin(phi_S) + np.sin(Omega_XH) * np.cos(theta_S) * np.cos(phi_S)
    Y_dot_l = cos_o_XH * np.sin(phi_S) + sin_o_XH * np.cos(theta_S) * np.cos(phi_S)
    #X_dot_m = - np.cos(i_JN) * np.sin(Omega_XH) * np.cos(phi_S) - np.cos(i_JN) * np.cos(Omega_XH) * np.cos(theta_S) * np.sin(phi_S) + np.sin(i_JN) * np.sin(theta_S) * np.sin(phi_S)
    X_dot_m = - cos_i_JN * sin_o_XH * np.cos(phi_S) - cos_i_JN * cos_o_XH * np.cos(theta_S) * np.sin(phi_S) + sin_i_JN * np.sin(theta_S) * np.sin(phi_S)
    #Y_dot_m = - np.cos(Omega_XH) * np.cos(phi_S) + np.sin(Omega_XH) * np.cos(theta_S) * np.sin(phi_S)
    Y_dot_m = - cos_o_XH * np.cos(phi_S) + sin_o_XH * np.cos(theta_S) * np.sin(phi_S)
    # Z_dot_l = - \sin\iota_{JN} \sin\Omega_{XH} \sin\phi_S + \sin\iota_{JN} \cos\Omega_{XH} \cos\theta_S \cos\phi_S + \cos\iota_{JN} \sin\theta_S \cos\phi_S
    Z_dot_l = -sin_i_JN * sin_o_XH * np.sin(phi_S) + sin_i_JN * cos_o_XH * np.cos(theta_S) * np.cos(phi_S) + cos_i_JN * np.sin(theta_S) * np.cos(phi_S)
    # Z_dot_m = \sin\iota_{JN} \sin\Omega_{XH} \cos\phi_S + \sin\iota_{JN} \cos\Omega_{XH} \cos\theta_S \sin\phi_S + \cos\iota_{JN} \sin\theta_S \sin\phi_S
    Z_dot_m = sin_i_JN * sin_o_XH * np.cos(phi_S) + sin_i_JN * cos_o_XH * np.cos(theta_S) * np.sin(phi_S) + cos_i_JN * np.sin(theta_S) * np.sin(phi_S)
    
    # if psi != 0:
    #     X' = cos(psi) * X + sin(psi) * Y
    #     Y' = -sin(psi) * X + cos(psi) * Y
    #X_dot_l = np.cos(psi) * X_dot_l + np.sin(psi) * Y_dot_l
    #Y_dot_l = -np.sin(psi) * X_dot_l + np.cos(psi) * Y_dot_l
    #X_dot_m = np.cos(psi) * X_dot_m + np.sin(psi) * Y_dot_m
    #Y_dot_m = -np.sin(psi) * X_dot_m + np.cos(psi) * Y_dot_m
    
    P_dot_l = cos_i_JN * X_dot_l - sin_i_JN * Z_dot_l
    P_dot_m = cos_i_JN * X_dot_m - sin_i_JN * Z_dot_m
    
    Q_dot_l = Y_dot_l
    Q_dot_m = Y_dot_m

    Fp = 0.5 * (((P_dot_l)**2) - ((Q_dot_l)**2) - ((P_dot_m)**2) + ((Q_dot_m)**2))
    Fc = ((P_dot_l * Q_dot_l) - (P_dot_m * Q_dot_m))
    
    return Fp, Fc



def get_inclination_IMRPhenomXPHM(params):
    """
    Get the inclination angle from the parameters for the IMRPhenomXPHM waveform model.

    Parameters:
    params (dict): A dictionary containing the parameters of the system.

    Returns:
    float: The inclination angle in radians.
    """
    # Extract the relevant parameters
    theta_LJ = params['theta_LJ']
    phi_JL0 = get_phi_JL0_IMRPhenomXPHM(params)
    cos_i_JN = params['cos_i_JN']
    sin_i_JN = params['sin_i_JN']

    # Calculate the inclination angle using the provided parameters
    # inclination_angle = np.arccos(np.sin(theta_LJ) * np.cos(phi_JL0) * np.sin(i_JN) + np.cos(theta_LJ) * np.cos(i_JN))
    inclination_angle = np.arccos(np.sin(theta_LJ) * np.cos(phi_JL0) * sin_i_JN + np.cos(theta_LJ) * cos_i_JN)

    return inclination_angle

def get_phi_JL0_IMRPhenomXPHM(params):
    """
    Get the reference phase phi_JL0 from the parameters for the IMRPhenomXPHM waveform model.

    Parameters:
    params (dict): A dictionary containing the parameters of the system.

    Returns:
    float: The reference phase phi_JL0 in radians.
    """
    # Extract the relevant parameters
    theta1 = params['theta_1']
    phi1 = params['phi_1']
    theta2 = params['theta_2']
    phi2 = params['phi_2']
    chi1 = params['chi1']
    chi2 = params['chi2']
    
    m1, m2 = get_m1_m2_IMRPhenomXPHM(params)
    
    # \tan^{-1} \left( \frac{m_1^2 \chi_{1} \sin\theta_1 \sin\phi_1 + m_2^2 \chi_{2} \sin\theta_2 \sin\phi_2}{m_1^2 \chi_{1} \sin\theta_1 \cos\phi_1 + m_2^2 \chi_{2} \sin\theta_2 \cos\phi_2} \right)
    
    num = (m1**2 * chi1 * np.sin(theta1) * np.sin(phi1)) + (m2**2 * chi2 * np.sin(theta2) * np.sin(phi2))
    den = (m1**2 * chi1 * np.sin(theta1) * np.cos(phi1)) + (m2**2 * chi2 * np.sin(theta2) * np.cos(phi2))
    
    phi_JL0 = np.arctan2(num, den)
    
    # This choice should maybe match up with gamma_P, as phi_JL_0 = gamma_P - pi/2, but I am not sure if this is the correct way to do it. I will check this later.
    
    return phi_JL0

def get_phi_ref_IMRPhenomXPHM(params):
    """
    Get the reference phase from the parameters for the IMRPhenomXPHM waveform model.

    Parameters:
    params (dict): A dictionary containing the parameters of the system.

    Returns:
    float: The reference phase in radians.
    """
    # Extract the relevant parameters
    theta_LJ = params['theta_LJ']
    phi_JL0 = get_phi_JL0_IMRPhenomXPHM(params)
    cos_i_JN = params['cos_i_JN']
    sin_i_JN = params['sin_i_JN']
    gamma_P = params['gamma_P']
    
    # \tan^{-1} \left(\frac{ (\cos\gamma_P \sin\phi_{JL_0} - \sin\gamma_P \cos\theta_{LJ} \cos\phi_{JL_0}) \sin\iota_{JN} +  \sin\gamma_P \sin\theta_{LJ} \cos\iota_{JN}}{(\sin\gamma_P \sin\phi_{JL_0} + \cos\gamma_P \cos\theta_{LJ} \cos\phi_{JL_0} )\sin\iota_{JN} - \cos\gamma_P \cos\theta_{LJ} \cos\iota_{JN}}\right)
    
    # Calculate the reference phase using the provided parameters
    num = (np.cos(gamma_P) * np.sin(phi_JL0) - np.sin(gamma_P) * np.cos(theta_LJ) * np.cos(phi_JL0)) * sin_i_JN + (np.sin(gamma_P) * np.sin(theta_LJ) * cos_i_JN)
    den = (np.sin(gamma_P) * np.sin(phi_JL0) + np.cos(gamma_P) * np.cos(theta_LJ) * np.cos(phi_JL0)) * sin_i_JN - (np.cos(gamma_P) * np.sin(theta_LJ) * cos_i_JN)
    #num = np.cos(theta_LJ) * np.sin(phi_LJ) * sin_i_JN - np.sin(theta_LJ) *  cos_i_JN
    #den = np.cos(phi_LJ) * sin_i_JN

    reference_phase = np.arctan2(num, den)
    
    reference_phase = reference_phase % (2 * np.pi)  # Ensure the phase is in the range [0, 2*pi)
    reference_phase = np.array(reference_phase)  # Convert to numpy array if it's not already

    return reference_phase[0] # Return the first element if it's an array, otherwise return the value

def cartesian_spin_components_IMRPhenomXPHM(params):
    """
    Calculate the Cartesian components of the spin vectors for the two black holes.

    Parameters:
    params (dict): A dictionary containing the parameters of the system.

    Returns:
    tuple: A tuple containing the Cartesian components of the spin vectors for both black holes.
    """
    # Extract the relevant parameters
    theta_1 = params['theta_1']
    phi_1 = params['phi_1']
    theta_2 = params['theta_2']
    phi_2 = params['phi_2']

    chi1 = params['chi1']
    chi2 = params['chi2']

    # Calculate the Cartesian components of the spin vectors
    S1_x = chi1 * np.sin(theta_1) * np.cos(phi_1)
    S1_y = chi1 * np.sin(theta_1) * np.sin(phi_1)
    S1_z = chi1 * np.cos(theta_1)
    
    S2_x = chi2 * np.sin(theta_2) * np.cos(phi_2)
    S2_y = chi2 * np.sin(theta_2) * np.sin(phi_2)
    S2_z = chi2 * np.cos(theta_2)

    return (S1_x, S1_y, S1_z), (S2_x, S2_y, S2_z)

def r_at_20_Hz(mcz_in_solar, eta, f_ref = 20.0):
    """
    Calculate the orbital separation at a reference frequency (default is 20 Hz) for a binary system with given chirp mass and mass ratio.
    
    Parameters:
    mcz_in_solar (float): The chirp mass of the binary system in solar masses.
    eta (float): The symmetric mass ratio of the binary system.
    f_ref (float): The reference frequency in Hz (default is 20 Hz).
    
    Returns:
    float: The orbital separation at the reference frequency in units of the total mass.
    """
    
    total_mass = (eta**(-0.6))*mcz_in_solar*solar_mass
    r20 = ((np.pi*f_ref*total_mass)**(-2/3))
    return r20

def get_theta_Omega_tilde(theta1, theta2, deltaphi, q, chi1, chi2, r):
    """
    Calculate the precession parameters theta_tilde and Omega_tilde using the provided parameters.
    
    Parameters:
    theta1 (float): The angle between the spin of the first black hole and the orbital angular momentum.
    theta2 (float): The angle between the spin of the second black hole and the orbital angular momentum.
    deltaphi (float): The difference in the azimuthal angles of the spins.
    q (float): The mass ratio of the binary system (m2/m1).
    chi1 (float): The dimensionless spin magnitude of the first black hole.
    chi2 (float): The dimensionless spin magnitude of the second black hole.
    r (float): The orbital separation in units of the total mass.

    Returns:
    tuple: A tuple containing theta_tilde and Omega_tilde.
    """
    
    chieff =  precession.eval_chieff(theta1, theta2, q, chi1, chi2)

    #calculate asymptotic angular momentum kappa
    kappa = (chi1 * np.cos(theta1) + q**2 * chi2 * np.cos(theta2) )/(1+q)**2 + \
            (chi1**2 + q**4 *chi2**2 + 2*chi1*chi2*q**2 * (np.cos(theta1)*np.cos(theta2) + np.cos(deltaphi)*np.sin(theta1)*np.sin(theta2))) / (2*q*(1+q)**2*r**(1/2))
            
            
    bracket_theta = precession.eval_bracket_theta(kappa, r, chieff, q, chi1, chi2)
    bracket_omega = precession.eval_bracket_omega(kappa, r, chieff, q, chi1, chi2)
    
    new_bracket_theta = bracket_theta * ((r/6.0)**(1.0/2.0)) * 10
    new_bracket_omega = bracket_omega * ((r/6.0)**(5.0/2.0))/(1000*solar_mass)
    
    return new_bracket_theta, new_bracket_omega

def get_IMRPhenomXPHM_dict_from_RP_params(rp_params, chi_1, chi_2, theta1, theta2, phi1, phi2):
    
    deltaphi = phi2 - phi1

    q_rp = (1 - 2 * rp_params['eta'] - np.sqrt(1 - 4 * rp_params['eta'])) / (2 * rp_params['eta'])

    separation_at_20 = r_at_20_Hz(rp_params['mcz']/solar_mass, rp_params['eta'])
    
    total_mass = (rp_params['eta']**(-0.6))*rp_params['mcz']
    
    theta_tilde_new, omega_tilde_new = get_theta_Omega_tilde(theta1, theta2, deltaphi, q_rp, chi_1, chi_2, separation_at_20)
    
    theta_tilde_new = np.nan_to_num(theta_tilde_new, nan=0.0)
    omega_tilde_new = np.nan_to_num(omega_tilde_new, nan=0.0)

    rp_params['theta_tilde'] = theta_tilde_new
    rp_params['omega_tilde'] = omega_tilde_new

    theta_LJ_0 = Regular_precession(rp_params).get_theta_LJ(f=20)
    phi_LJ_0 = Regular_precession(rp_params).get_phi_LJ(f=20)

    cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = Regular_precession(rp_params).precession_angles()

    parameters_RP_to_IMRPhenomXPHM = {
        'chi1': chi_1,
        'chi2': chi_2,
        'chirp_mass': rp_params['mcz']/solar_mass,
        'eta': rp_params['eta'],
        'theta_S': rp_params['theta_S'],
        'phi_S': rp_params['phi_S'],
        'theta_LJ': theta_LJ_0,
        'phi_LJ': phi_LJ_0,
        'cos_i_JN': cos_i_JN,
        'sin_i_JN': sin_i_JN,
        'cos_o_XH': cos_o_XH,
        'sin_o_XH': sin_o_XH,
        'theta_1': theta1,
        'theta_2': theta2,
        'phi_1': phi1,
        'phi_2': phi2,
        'gamma_P': rp_params['gamma_P']
    }

    return parameters_RP_to_IMRPhenomXPHM

def get_hstrain_from_modes(rp_params, chi_1, chi_2, theta1, theta2, phi1, phi2, fmax, psi, return_freqs=False, mode_list=None):
    
    cos_theta_JN, _, _, _ = Regular_precession(rp_params).precession_angles()
    
    theta_JN = np.arccos(cos_theta_JN)
    
    params_IMRPhenomXPHM = get_IMRPhenomXPHM_dict_from_RP_params(rp_params, chi_1, chi_2, theta1, theta2, phi1, phi2)
    
    m1, m2 = get_m1_m2_IMRPhenomXPHM(params_IMRPhenomXPHM)
    Fp, Fc = get_Fp_Fc_IMRPhenomXPHM(params_IMRPhenomXPHM, psi)
    inclination = get_inclination_IMRPhenomXPHM(params_IMRPhenomXPHM)
    phi_ref = get_phi_ref_IMRPhenomXPHM(params_IMRPhenomXPHM)
    S1_vec, S2_vec = cartesian_spin_components_IMRPhenomXPHM(params_IMRPhenomXPHM)

    param_dict = {
        'mass1' : float(m1) * Msun,
        'mass2' : float(m2) * Msun,
        'spin1x' : float(S1_vec[0]) * dimensionless,
        'spin1y' : float(S1_vec[1]) * dimensionless,
        'spin1z' : float(S1_vec[2]) * dimensionless,
        'spin2x' : float(S2_vec[0]) * dimensionless,
        'spin2y' : float(S2_vec[1]) * dimensionless,
        'spin2z' : float(S2_vec[2]) * dimensionless,
        'deltaF' : 0.05 * Hz,
        'f22_start' : 20.0 * Hz,
        'f22_ref': 20.0 * Hz,
        'f_max': fmax * Hz,
        'phi_ref' : float(phi_ref) * radians,
        'distance' : float((rp_params['dist']/giga_parsec) * 1000) * Mpc,
        'inclination' : float(inclination) * radians,
        'condition' : 0
    }

    approximant = 'IMRPhenomXPHM'
    generator = gwsignal.core.waveform.LALCompactBinaryCoalescenceGenerator(approximant)
    
    hlm = gwsignal.core.waveform.GenerateFDModes(param_dict, generator)

    freqs = None
    hlm_by_lm = {}

    for key, val in hlm.items():
        if isinstance(key, str):
            if 'freq' in key.lower():          # 'frequency_array' or similar
                freqs = np.array(val)
            continue

        # Non-string key -> a mode
        try:
            lm = (int(key.l), int(key.m))      # SpinWeightedSphericalHarmonicMode object
        except AttributeError:
            try:
                lm = (int(key[0]), int(key[1]))  # tuple-like
            except Exception:
                print(f"Skipping unrecognised key: {key!r}")
                continue

        hlm_by_lm[lm] = np.array(val)

    if freqs is None:
        f_start = float(param_dict['f22_start'].value)
        f_max   = float(param_dict['f_max'].value)
        delta_f = float(param_dict['deltaF'].value)
        print("Warning: Frequency array not found in hlm output, constructing from parameters.")
        freqs   = np.arange(f_start, f_max + delta_f, delta_f)

    if mode_list is not None:
        requested = set(tuple(m) for m in mode_list)
        hlm_by_lm = {lm: v for lm, v in hlm_by_lm.items() if lm in requested}
        #print("Using modes:", list(hlm_by_lm.keys()))

    theta = theta_JN * u.rad if theta_JN is not None else float(param_dict['inclination'].value) * u.rad
    phi_val = 0 * u.rad  # Default phi = 0, can be modified if needed
    phi = float(phi_val.value) * u.rad
    
    hp = 0.
    hc = 0.

    for (l, m), hlm_arr in hlm_by_lm.items():
        mode_fn = gwsignal.core.gw.SpinWeightedSphericalHarmonicMode(-2, l, m)
        ylm = mode_fn(theta, phi)  # -2 for spin weight, l and m as mode indices

        hp += 0.5 * (ylm * hlm_arr + np.conj(ylm) * np.conj(hlm_arr)[::-1])
        hc += 1j*0.5 * (ylm * hlm_arr - np.conj(ylm) * np.conj(hlm_arr)[::-1])
        
    hp_array = np.array(hp)
    hc_array = np.array(hc)
    
    if return_freqs:
        return freqs, Fp * hp_array + Fc * hc_array
    
    else:
        return Fp * hp_array + Fc * hc_array


def get_hstrain_NP_align_from_modes(rp_params, chi_1, chi_2, theta1, theta2, phi1, phi2, fmax, psi, return_freqs=False, mode_list=None):
    
    cos_theta_JN, _, _, _ = Regular_precession(rp_params).precession_angles()
    
    theta_JN = np.arccos(cos_theta_JN)
    
    params_IMRPhenomXPHM = get_IMRPhenomXPHM_dict_from_RP_params(rp_params, chi_1, chi_2, theta1, theta2, phi1, phi2)
    
    m1, m2 = get_m1_m2_IMRPhenomXPHM(params_IMRPhenomXPHM)
    Fp, Fc = get_Fp_Fc_IMRPhenomXPHM(params_IMRPhenomXPHM, psi)
    inclination = get_inclination_IMRPhenomXPHM(params_IMRPhenomXPHM)
    phi_ref = get_phi_ref_IMRPhenomXPHM(params_IMRPhenomXPHM)
    S1_vec, S2_vec = cartesian_spin_components_IMRPhenomXPHM(params_IMRPhenomXPHM)

    param_dict = {
        'mass1' : float(m1) * Msun,
        'mass2' : float(m2) * Msun,
        'spin1x' : 0 * dimensionless,  
        'spin1y' : 0 * dimensionless,
        'spin1z' : float(S1_vec[2]) * dimensionless,
        'spin2x' : 0 * dimensionless,
        'spin2y' : 0 * dimensionless,
        'spin2z' : float(S2_vec[2]) * dimensionless,
        'deltaF' : 0.05 * Hz,
        'f22_start' : 20.0 * Hz,
        'f22_ref': 20.0 * Hz,
        'f_max': fmax * Hz,
        'phi_ref' : float(phi_ref) * radians,
        'distance' : float((rp_params['dist']/giga_parsec) * 1000) * Mpc,
        'inclination' : float(inclination) * radians,
        'condition' : 0
    }

    approximant = 'IMRPhenomXPHM'
    generator = gwsignal.core.waveform.LALCompactBinaryCoalescenceGenerator(approximant)
    
    hlm = gwsignal.core.waveform.GenerateFDModes(param_dict, generator)

    freqs = None
    hlm_by_lm = {}

    for key, val in hlm.items():
        if isinstance(key, str):
            if 'freq' in key.lower():          # 'frequency_array' or similar
                freqs = np.array(val)
            continue

        # Non-string key -> a mode
        try:
            lm = (int(key.l), int(key.m))      # SpinWeightedSphericalHarmonicMode object
        except AttributeError:
            try:
                lm = (int(key[0]), int(key[1]))  # tuple-like
            except Exception:
                print(f"Skipping unrecognised key: {key!r}")
                continue

        hlm_by_lm[lm] = np.array(val)

    if freqs is None:
        f_start = float(param_dict['f22_start'].value)
        f_max   = float(param_dict['f_max'].value)
        delta_f = float(param_dict['deltaF'].value)
        print("Warning: Frequency array not found in hlm output, constructing from parameters.")
        freqs   = np.arange(f_start, f_max + delta_f, delta_f)

    if mode_list is not None:
        requested = set(tuple(m) for m in mode_list)
        hlm_by_lm = {lm: v for lm, v in hlm_by_lm.items() if lm in requested}
        #print("Using modes:", list(hlm_by_lm.keys()))

    theta = theta_JN * u.rad if theta_JN is not None else float(param_dict['inclination'].value) * u.rad
    phi_val = 0 * u.rad  # Default phi = 0, can be modified if needed
    phi = float(phi_val.value) * u.rad
    
    hp = 0.
    hc = 0.

    for (l, m), hlm_arr in hlm_by_lm.items():
        mode_fn = gwsignal.core.gw.SpinWeightedSphericalHarmonicMode(-2, l, m)
        ylm = mode_fn(theta, phi)  # -2 for spin weight, l and m as mode indices

        hp += 0.5 * (ylm * hlm_arr + np.conj(ylm) * np.conj(hlm_arr)[::-1])
        hc += 1j*0.5 * (ylm * hlm_arr - np.conj(ylm) * np.conj(hlm_arr)[::-1])
        
    hp_array = np.array(hp)
    hc_array = np.array(hc)
    
    if return_freqs:
        return freqs, Fp * hp_array + Fc * hc_array
    
    else:
        return Fp * hp_array + Fc * hc_array


def get_RP_Waveform(rp_params, chi_1, chi_2, theta1, theta2, phi1, phi2, fmax, psi = 0, return_freqs=True):
    
    """
    Generate the precessing waveform using the Regular Precession model and return the strain.
    Parameters:
    rp_params (dict): A dictionary containing the parameters for the Regular Precession model.
    chi_1 (float): The dimensionless spin magnitude of the first black hole.
    chi_2 (float): The dimensionless spin magnitude of the second black hole.
    theta1 (float): The angle between the spin of the first black hole and the orbital angular momentum.
    theta2 (float): The angle between the spin of the second black hole and the orbital angular momentum.
    phi1 (float): The azimuthal angle of the spin of the first black hole.
    phi2 (float): The azimuthal angle of the spin of the second black hole.
    fmax (float): The maximum frequency up to which the waveform should be generated.
    psi (float): The polarization angle (default is 0).
    return_freqs (bool): Whether to return the frequency array along with the waveform (default is True).
    
    Returns:
    If return_freqs is True, returns a tuple (f_range, waveform_RP) where f_range is the array of frequencies and waveform_RP is the corresponding strain values.
    If return_freqs is False, returns only the waveform_RP array.
    """
    
    
    f_cut = Regular_precession(rp_params).get_f_cut()
    f_range = np.arange(20, f_cut, 0.05)

    separation_at_20 = r_at_20_Hz(rp_params['mcz']/solar_mass, rp_params['eta'])
    q_rp = (1 - 2 * rp_params['eta'] - np.sqrt(1 - 4 * rp_params['eta'])) / (2 * rp_params['eta'])
    total_mass = Regular_precession(rp_params).get_total_mass()
    theta_tilde, omega_tilde = get_theta_Omega_tilde(theta1, theta2, np.abs(phi2 - phi1), q_rp, chi_1, chi_2, separation_at_20)

    theta_tilde = np.nan_to_num(theta_tilde, nan=0.0)
    omega_tilde = np.nan_to_num(omega_tilde, nan=0.0)

    rp_params['theta_tilde'] = theta_tilde
    rp_params['omega_tilde'] = omega_tilde

    waveform_RP = Regular_precession(rp_params).precessing_strain(f_range)
    
    if return_freqs:
        return f_range, waveform_RP
    
    return waveform_RP


def to_pycbc_fs(h_arr, f_arr, delta_f):
    """Wrap into a pycbc FrequencySeries zero-padded from f=0."""
    f_start_idx = int(round(f_arr[0] / delta_f))
    padded = np.zeros(f_start_idx + len(h_arr), dtype=complex)
    padded[f_start_idx:] = h_arr
    return FrequencySeries(padded, delta_f=delta_f)

def mismatch_IMRPhenomXPHM(h1_arr, h2_arr, f1, f2, f_cut=None, delta_f=0.05):
    """
    Calculate the mismatch between the IMRPhenomXPHM waveform and the Regular Precession waveform.

    Parameters:
    h1_arr (numpy array): The complex strain values of the first waveform.
    h2_arr (numpy array): The complex strain values of the second waveform.
    f1 (numpy array): The frequency grid for the first waveform.
    f2 (numpy array): The frequency grid for the second waveform.
    f_cut (float, optional): The cutoff frequency. If None, no cutoff is applied.
    delta_f (float): The frequency resolution.

    Returns:
    float: The mismatch value between 0 and 1, where 0 indicates perfect agreement and 1 indicates complete disagreement.
    """
    
    h1_arr = np.array(h1_arr)
    h2_arr = np.array(h2_arr)

    f_low  = 20.0
    f_high = min(f1[-1], f2[-1])
    if f_cut is not None:
        f_high = min(f_high, f_cut)

    # Slice h1 to [f_low, f_high]
    mask1  = (f1 >= f_low) & (f1 <= f_high)
    f1_cut = f1[mask1]
    h1_cut = h1_arr[mask1]

    # Interpolate h2 onto f1's grid
    re_interp = interp1d(f2, np.real(h2_arr), bounds_error=False, fill_value=0.0)
    im_interp = interp1d(f2, np.imag(h2_arr), bounds_error=False, fill_value=0.0)
    h2_cut = re_interp(f1_cut) + 1j * im_interp(f1_cut)

    # Wrap as pycbc FrequencySeries (zero-padded from f=0)
    h1_fs = to_pycbc_fs(h1_cut, f1_cut, delta_f)
    h2_fs = to_pycbc_fs(h2_cut, f1_cut, delta_f)
    Sn_fs = to_pycbc_fs(Sn(f1_cut, delta_f=delta_f), f1_cut, delta_f)

    match, _, _ = optimized_match(
        h1_fs, h2_fs, Sn_fs,
        low_frequency_cutoff=f_low,
        high_frequency_cutoff=f_high,
        return_phase=True
    )
    
    return 1 - match