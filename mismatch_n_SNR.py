
import numpy as np
from scipy.integrate import simps

from pycbc.filter import optimized_match
from pycbc.types import FrequencySeries

error_handler = np.seterr(invalid="raise")

solar_mass = 4.92624076 * 1e-6              # [solar_mass] = sec
giga_parsec = 1.02927125 * 1e17             # [giga_parsec] = sec
year = 31557600                             # [year] = sec

# Importing the regular precession class
from regular_precession import *

# Importing system parameters
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


def ideal_params(rp_params, min_gammaP):
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
    
    rp_params = ideal_params(rp_params, min_gammaP)

    if return_tc_phic:
        min_ind = np.argmin(mismatch_arr)
        min_ind_a = ind_arr[min_ind]
        min_phase_a = phase_m_arr[min_ind]

        return min_mismatch, max_mismatch, min_gammaP, max_gammaP, min_ind_a, min_phase_a

    return min_mismatch, max_mismatch, min_gammaP, max_gammaP



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


# Some utilities

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