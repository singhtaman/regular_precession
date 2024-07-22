
import numpy as np
from scipy.integrate import simps

from pycbc.filter import optimized_match
from pycbc.types import FrequencySeries
import precession

error_handler = np.seterr(invalid="raise")

solar_mass = 4.92624076 * 1e-6              # [solar_mass] = sec
giga_parsec = 1.02927125 * 1e17             # [giga_parsec] = sec
year = 31557600                             # [year] = sec

mass2sec = 4.92624076e-6

#Regular precessing class
from regular_precession import *

#Importing the required functions for params
from systems_lib import *

# For cosmology: converting redshift to luminosity distance
from astropy.cosmology import FlatLambdaCDM


def Sn(f, delta_f=0.25, frequencySeries=True):
    """ ALIGO noise curve from arXiv:0903.0338
    _______________________________________________________________________________________
    Parameters used:
    f : array_like : frequency array
    delta_f : float : frequency step
    frequencySeries : bool : return FrequencySeries object
    _______________________________________________________________________________________
    Returns:
    Sn_val : Noise curve values in the form of
        FrequencySeries : if frequencySeries is True
        array_like : if frequencySeries is False
    _______________________________________________________________________________________
    Noise curve array for ALIGO
    """
    Sn_val = np.zeros_like(f)
    fs = 20
    for i in range(len(f)):
        if f[i] < fs:
            Sn_val[i] = np.inf
        else:
            S0 = 1E-49
            f0 = 215
            Sn_temp = np.power(f[i]/f0, -4.14) - 5 * np.power(f[i]/f0, -2) + 111 * ((1 - np.power(f[i]/f0, 2) + 0.5 * np.power(f[i]/f0, 4)) / (1 + 0.5 * np.power(f[i]/f0, 2)))
            Sn_val[i] = Sn_temp * S0
    if frequencySeries:
        return FrequencySeries(Sn_val, delta_f=delta_f)
    return Sn_val


def opt_mismatch_gammaP(rp_params, np_params, size_of_gammaP_arr):
    """ Optimal mismatch for each gammaP
    _______________________________________________________________________________________
    Parameters used:
    rp_params : dict : Regular precession parameters
    np_params : dict : Non-precessing parameters
    size_of_gammaP_arr : int : size of gammaP array
    _______________________________________________________________________________________
    Returns:
    gammaP_arr : array_like : gammaP array
    mismatch_arr : array_like : mismatch array
    ind_arr : array_like : index array
    phase_m_arr : array_like : phase array
    _______________________________________________________________________________________
    Mismatches between RP and NP signals for each gammaP
    """
    f_cut = Regular_precession(rp_params).get_f_cut()
    f_range = np.arange(20, f_cut, 0.25)
    
    delta_f = 0.25

    gammaP_arr = np.linspace(0, 2*np.pi, size_of_gammaP_arr)
    mismatch_arr = [] 
    ind_arr = []
    phase_m_arr = []
    
    psd = Sn(f_range, delta_f)
    
    for i in range(size_of_gammaP_arr):
        
        rp_params['gamma_P'] = gammaP_arr[i]
        
        precession_initial = Regular_precession(rp_params)
        NP_initial = Regular_precession(np_params)
        
        precessing_signal = precession_initial.precessing_strain(f_range, delta_f)
        non_precessing_signal = NP_initial.precessing_strain(f_range, delta_f)
        
        match, ind, phase_m = optimized_match(precessing_signal, non_precessing_signal, psd, return_phase=True)
        
        mismatch_arr.append(1 - match)
        ind_arr.append(ind)
        phase_m_arr.append(phase_m)
    
    return gammaP_arr, mismatch_arr, ind_arr, phase_m_arr



def opt_mismatch_extremas_gammaP(rp_params, np_params, size_of_gammaP_arr, return_tc_phic):
    """ Optimal mismatch extremas for array of gammaP
    _______________________________________________________________________________________
    Parameters used:
    rp_params : dict : Regular precession parameters
    np_params : dict : Non-precessing parameters
    size_of_gammaP_arr : int : size of gammaP array
    return_tc_phic : bool : return tc and phic
    _______________________________________________________________________________________
    Returns:
    min_mismatch : float : min mismatch
    max_mismatch : float : max mismatch
    min_gammaP : float : min gammaP
    max_gammaP : float : max gammaP
    and if return_tc_phic is True then - 
        min_ind : int : min index
        min_phase : float : min phase
    _______________________________________________________________________________________
    Extremas of mismatch values for array of gamma_P along with tc and phic if needed 
    """
    
    gammaP_arr, mismatch_arr, ind_arr, phase_m_arr = opt_mismatch_gammaP(rp_params, np_params, size_of_gammaP_arr)
    min_mismatch = np.min(mismatch_arr)
    max_mismatch = np.max(mismatch_arr)
    min_gammaP = gammaP_arr[np.argmin(mismatch_arr)]
    max_gammaP = gammaP_arr[np.argmax(mismatch_arr)]

    if return_tc_phic:
        min_ind = np.argmin(mismatch_arr)
        min_ind_a = ind_arr[min_ind]
        min_phase_a = phase_m_arr[min_ind]

        return min_mismatch, max_mismatch, min_gammaP, max_gammaP, min_ind_a, min_phase_a

    return min_mismatch, max_mismatch, min_gammaP, max_gammaP


def ideal_param(rp_params, min_gammaP):
    """ Ideal parameters for a given gamma_P_min
    _______________________________________________________________________________________
    Parameters used:
    rp_params : dict : Regular precession parameters
    min_gammaP : float : min gammaP
    _______________________________________________________________________________________
    Returns:
    rp_params : dict : Regular precession parameters (ideal ones)
    _______________________________________________________________________________________
    Sets gammaP to min_gammaP
    """
    rp_params["gamma_P"] = min_gammaP
    return rp_params



def get_SNRs(rp_params, np_params, redshift = 100)->dict:
    """ Get SNRs for RP signal and NP signal
    _______________________________________________________________________________________
    Parameters used:
    rp_params : dict : Regular precession parameters
    np_params : dict : Non-precessing parameters
    _______________________________________________________________________________________
    Returns:
    SNR_dict : dict : SNR values for RP signal, NP signal and cross terms
    _______________________________________________________________________________________
    SNR values for RP signal and NP signal along with cross terms -- Cross terms don't yield SNR
    """
    f_cut = Regular_precession(rp_params).get_f_cut()
    f_range = np.arange(20, f_cut, 0.25)

    psd = Sn(f_range)
    
    if redshift == 100: # mask for no redshift, only luminosity distance -- for testing
        rp_signal = Regular_precession(rp_params).precessing_strain(f_range)
        np_signal = Regular_precession(np_params).precessing_strain(f_range)
    else:
        rp_params["dist"] = redshift_to_luminosity_distance(redshift)*giga_parsec
        np_params["dist"] = redshift_to_luminosity_distance(redshift)*giga_parsec
        rp_signal = Regular_precession(rp_params).precessing_strain(f_range)
        np_signal = Regular_precession(np_params).precessing_strain(f_range)

    integrand_RP = np.conj(rp_signal)*rp_signal/psd
    integrated_inner_product_RP = simps(integrand_RP, f_range)
    SNR_rp = np.sqrt(4*integrated_inner_product_RP.real) # SNR for RP signal

    integrand_NP = np.conj(np_signal)*np_signal/psd
    integrated_inner_product_NP = simps(integrand_NP, f_range)
    SNR_np = np.sqrt(4*integrated_inner_product_NP.real) # SNR for NP signal

    integrand_cross1 = np.conj(rp_signal)*np_signal/psd
    integrated_inner_product_cross1 = simps(integrand_cross1, f_range)
    SNR_cross1 = np.sqrt(4*abs(integrated_inner_product_cross1.real)) # Cross term inner product -- doesn't yield SNR

    integrand_cross2 = np.conj(np_signal)*rp_signal/psd
    integrated_inner_product_cross2 = simps(integrand_cross2, f_range)
    SNR_cross2 = np.sqrt(4*abs(integrated_inner_product_cross2.real)) # Another cross term inner product -- doesn't yield SNR

    return {"SNR_RPRP": SNR_rp, 
            "SNR_NPNP": SNR_np, 
            "SNR_RPNP": SNR_cross1,
            "SNR_NPRP": SNR_cross2
            } 



def lindblom_inequality(rp_params, np_params, size_of_gammaP_arr, rp_SNR_term = True, np_SNR_term = False):
    """ Lindblom inequality
    ____________________________________________________________________________________________
    Parameters used:
    rp_params : dict : Regular precession parameters
    np_params : dict : Non-precessing parameters
    size_of_gammaP_arr : int : size of gammaP array
    cross_SNR_term : bool : cross SNR term
    rp_SNR_term : bool : RP SNR term
    np_SNR_term : bool : NP SNR term
    ____________________________________________________________________________________________
    Returns:
    lindblom : float : lindblom inequality value
    ____________________________________________________________________________________________
    """
    mismatch_min = opt_mismatch_extremas_gammaP(rp_params, np_params, size_of_gammaP_arr, return_tc_phic = False)
    SNR_dict = get_SNRs(rp_params, np_params)
    
    if rp_SNR_term:
        SNR_rp = SNR_dict["SNR_RPRP"]
        lindblom = mismatch_min - 1/(2*(SNR_rp**2))
        return lindblom
    
    elif np_SNR_term:
        SNR_np = SNR_dict["SNR_NPNP"]
        lindblom = mismatch_min - 1/(2*(SNR_np**2))
        return lindblom
    
    else:
        print("Please select a term to calculate the lindblom inequality")
        print("copy this in the function call : 'rp_SNR_term = true, np_SNR_term = False' for RP SNR after rp, np and size of gammaP array Parameters used")
        return None




def same_total_mass_diff_chirp(q, mcz, rp_params, np_params):
    """ Varying q while keeping total mass constant, returns new param dict and new chirp mass along with eta
    _______________________________________________________________________________________
    Parameters used:
    q : float : mass ratio
    mcz : float : chirp mass
    rp_params : dict : Regular precession parameters
    np_params : dict : Non-precessing parameters
    _______________________________________________________________________________________
    Returns:
    rp_params : dict : Regular precession parameters
    np_params : dict : Non-precessing parameters
    mcz_new : float : new chirp mass
    eta : float : eta
    _______________________________________________________________________________________
    """
    eta = q/(1+q)**2
    default_mcz = mcz
    total_mass_q1 = mcz/(0.25**0.6)
    mcz_new = total_mass_q1 * eta**0.6
    rp_params['eta'] = np_params['eta'] = eta
    rp_params['mcz'] = np_params['mcz'] = mcz_new*solar_mass
    return rp_params, np_params, mcz_new, eta

def random_angles_generator():
    """ Random angles generator
    _______________________________________________________________________________________
    Parameters used:
    None
    _______________________________________________________________________________________
    Returns:
    theta_J: float : thetaJ
    phi_J: float : phiJ
    theta_S: float : thetaS
    phi_S: float : phiS
    _______________________________________________________________________________________
    Generates random sky locations and binary orientations
    """
    n = 4
    get_random_numbers = np.random.random(n)
    
    # Random angles for J
    cos_theta_J = 2*get_random_numbers[0] - 1
    thetaJ = np.arccos(cos_theta_J)
    phiJ = 2*np.pi*get_random_numbers[1]
    
    # Random angles for S
    cos_theta_S = 2*get_random_numbers[2] - 1
    thetaS = np.arccos(cos_theta_S)
    phiS = 2*np.pi*get_random_numbers[3]
    
    return thetaJ, phiJ, thetaS, phiS



def random_J_angles_generator():
    """ Random J angles generator
    _______________________________________________________________________________________
    Parameters used:
    None
    _______________________________________________________________________________________
    Returns:
    thetaJ : float : thetaJ
    phiJ : float : phiJ
    _______________________________________________________________________________________
    Generates random binary orientations
    """
    n = 2
    get_random_numbers = np.random.random(n)
    
    # Random angles for J
    
    cos_theta_J = 2*get_random_numbers[0] - 1
    thetaJ = np.arccos(cos_theta_J)
    phiJ = 2*np.pi*get_random_numbers[1]
    
    return thetaJ, phiJ

def random_N_angles_generator():
    """ Random N angles generator
    _______________________________________________________________________________________
    Parameters used:
    None
    _______________________________________________________________________________________
    Returns:
    thetaS : float : thetaS
    phiS : float : phiS
    _______________________________________________________________________________________
    Generates random sky location angles
    """
    n = 2
    get_random_numbers = np.random.random(n)
    
    # Random angles for N
    
    cos_theta_S = 2*get_random_numbers[0] - 1
    thetaS = np.arccos(cos_theta_S)
    phiS = 2*np.pi*get_random_numbers[1]
    
    return thetaS, phiS


def get_param_distributions(q, s, N, r):
    """Gets precession parameters for a population of binaries
    _______________________________________________________________________________________
    Parameters used:
    q : float : mass ratio
    s : float : spins - equal spin case
    N : int : number of binaries in the population
    r : float : binary separation we want our distribution values at
    _______________________________________________________________________________________
    Returns:
    bracket_theta_arr : array_like : bracket theta array
    bracket_omega_arr : array_like : bracket omega array
    delta_theta_arr : array_like : delta theta array
    delta_omega_arr : array_like : delta omega array
    little_omega_arr : array_like : little omega array
    _______________________________________________________________________________________
    Calculates distributions of all 5 precession parameters for a population of binaries
    """
    vec = precession.isotropic_angles(N)

    # N binaries with random theta1, theta2 and deltaphi
    theta1_arr, theta2_arr, deltaphi_arr = vec[0], vec[1], vec[2]

    bracket_theta_arr = []
    bracket_omega_arr = []
    delta_theta_arr = []
    delta_omega_arr = []
    little_omega_arr = []


    for it in range(N) :
        # theta1, theta2 and deltaphi from the arrays
        theta1 = theta1_arr[it]
        theta2 = theta2_arr[it]
        deltaphi = deltaphi_arr[it]

        #assign spin s to both binaries
        chi1 = chi2 = s

        #calculate chi_eff
        chieff =  precession.eval_chieff(theta1, theta2, q, chi1, chi2)

        #calculate asymptotic angular momentum kappa
        kappa = (chi1 * np.cos(theta1) + q**2 * chi2 * np.cos(theta2) )/(1+q)**2 + \
                (chi1**2 + q**4 *chi2**2 + 2*chi1*chi2*q**2 * (np.cos(theta1)*np.cos(theta2) + np.cos(deltaphi)*np.sin(theta1)*np.sin(theta2))) / (2*q*(1+q)**2*r**(1/2))

        #based on above values calculate bracket omega, bracket theta, delta theta, delta omega, nut_freq
        
        bracket_theta = precession.eval_bracket_theta(kappa, r, chieff, q, chi1, chi2)
        bracket_omega = precession.eval_bracket_omega(kappa, r, chieff, q, chi1, chi2)
        delta_theta = precession.eval_delta_theta(kappa, r, chieff, q, chi1, chi2)
        delta_omega = precession.eval_delta_omega(kappa, r, chieff, q, chi1, chi2)
        little_omega = precession.eval_nutation_freq(kappa, r, chieff, q, chi1, chi2)


        #appending all the values to corresponding arrays
        bracket_theta_arr.append(bracket_theta)
        bracket_omega_arr.append(bracket_omega)
        delta_theta_arr.append(delta_theta)
        delta_omega_arr.append(delta_omega)
        little_omega_arr.append(little_omega)

    return bracket_theta_arr, bracket_omega_arr, delta_theta_arr, delta_omega_arr, little_omega_arr


def random_precession_params_iso(chi1 = 1, chi2 = 1, q = 1, r = 6):
    """ Random precession params:
    _______________________________________________________________________________________
    Parameters used:
    chi1 : float : chi1
    chi2 : float : chi2
    q : float : mass ratio
    r : float : binary separation
    _______________________________________________________________________________________
    Returns:
    theta_tilde : float : dimensionless theta
    omega_tilde : float : dimensionless omega
    bracket_theta : float : bracket theta
    bracket_omega : float : bracket omega
    _______________________________________________________________________________________
    draws theta_1, theta_2, deltaphi randomly and calculates bracket theta, bracket omega using the PRECESSION code
    returns dimensionless theta, omega and precessional averaged theta and omega
    """
    n = 3
    get_random_numbers_p = np.random.random(n)

    cos_theta_1 = 2*get_random_numbers_p[0] - 1
    theta_1 = np.arccos(cos_theta_1)
    cos_theta_2 = 2*get_random_numbers_p[1] - 1
    theta_2 = np.arccos(cos_theta_2)

    deltaphi = 2*np.pi*get_random_numbers_p[2]

    chieff = precession.eval_chieff(theta_1, theta_2, q, chi1, chi2)

    kappa = (chi1 * np.cos(theta_1) + q**2 * chi2 * np.cos(theta_2) )/(1+q)**2 + \
            (chi1**2 + q**4 *chi2**2 + 2*chi1*chi2*q**2 * (np.cos(theta_1)*np.cos(theta_2) + np.cos(deltaphi)*np.sin(theta_1)*np.sin(theta_2))) / (2*q*(1+q)**2*r**(1/2))
    
    bracket_theta = precession.eval_bracket_theta(kappa, r, chieff, q, chi1, chi2)
    bracket_omega = precession.eval_bracket_omega(kappa, r, chieff, q, chi1, chi2)
    
    theta_tilde = bracket_theta*10
    omega_tilde = bracket_omega/(1000*mass2sec)

    return theta_tilde, omega_tilde, bracket_theta, bracket_omega


def vary_params(rp_params, np_params, set_new_line_of_sight = False, set_new_precession_params = False, do_you_want_iso = True):
    """ Get random RP params: Varies J angles first, then N angles and then precession params
    _______________________________________________________________________________________
    Parameters used:
    rp_params : dict : Regular precession parameters
    np_params : dict : Non-precessing parameters
    set_new_line_of_sight : bool : set new line of sight
    set_new_precession_params : bool : set new precession params
    _______________________________________________________________________________________
    Returns:
    rp_params : dict : Regular precession parameters
    np_params : dict : Non-precessing parameters
    _______________________________________________________________________________________
    """
    theta_J, phi_J = random_J_angles_generator()
    if set_new_line_of_sight:
        theta_S, phi_S = random_N_angles_generator()
    else:
        theta_S = rp_params["theta_S"]
        phi_S = rp_params["phi_S"]
    if set_new_precession_params:
        if do_you_want_iso:
            theta_tilde, omega_tilde, bracket_theta, bracket_omega = random_precession_params_iso()
        
    else:
        theta_tilde = rp_params["theta_tilde"]
        omega_tilde = rp_params["omega_tilde"]
    
    rp_params["theta_J"] = theta_J
    rp_params["phi_J"] = phi_J
    np_params["theta_J"] = theta_J
    np_params["phi_J"] = phi_J
    
    rp_params["theta_S"] = theta_S
    rp_params["phi_S"] = phi_S
    np_params["theta_S"] = theta_S
    np_params["phi_S"] = phi_S
    
    rp_params["theta_tilde"] = theta_tilde
    rp_params["omega_tilde"] = omega_tilde
    
    return rp_params, np_params



def redshift_to_luminosity_distance(z)->float:
    """
    Convert redshift to luminosity distance in standard Lambda CDM cosmology.
    _______________________________________________________________________________________    
    Parameters used:
    z : float : redshift
    _______________________________________________________________________________________
    Returns:    
    luminosity_distance : float : luminosity distance in Gpc
    """
    # Define the cosmology
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315) # Planck 2018 results
    luminosity_distance = cosmo.luminosity_distance(z).to('Gpc').value
    return luminosity_distance



def redshifted_new_params(z, rp_params)->dict:
    """
    Get redshifted params
    _______________________________________________________________________________________    
    Parameters used:
    z : float : redshift
    rp_params : dict : Regular precession parameters
    _______________________________________________________________________________________
    Returns:    
    rp_params : dict : Regular precession parameters
    """
    rp_params["dist"] = redshift_to_luminosity_distance(z)*giga_parsec
    old_chirp = rp_params["mcz"]
    #new chirp mass
    rp_params["mcz"] = old_chirp*(1+z)
    return rp_params

