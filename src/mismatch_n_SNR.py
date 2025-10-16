"""
Mismatch and SNR calculations for regular precession waveforms.

This module provides functions to compute mismatches, SNRs, and related quantities for precessing and non-precessing systems.
"""

import numpy as np
from scipy.integrate import simps
import sys

from pycbc.filter import optimized_match
from pycbc.types import FrequencySeries

error_handler = np.seterr(invalid="raise")

# Importing the regular precession class
sys.path.insert(0, "../")
from src.regular_precession import *

# Importing system parameters
from src.systems_lib import *

# For cosmology: converting redshift to luminosity distance
from astropy.cosmology import FlatLambdaCDM


def Sn(f, delta_f=0.25, frequencySeries=True):
    """
    Compute the ALIGO noise curve from arXiv:0903.0338.

    Parameters
    ----------
    f : array_like
        Frequency array.
    delta_f : float, optional
        Frequency step (default 0.25).
    frequencySeries : bool, optional
        If True, return a FrequencySeries object (default True).

    Returns
    -------
    Sn_val : FrequencySeries or array_like
        Noise curve values as FrequencySeries if frequencySeries is True, else as array.
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

def get_mismatch_rp(source, template, update_tc_phic=True):
    """
    Calculate the mismatch between a source and template waveform.

    Parameters
    ----------
    source : dict
        Source system parameters.
    template : dict
        Template system parameters.
    update_tc_phic : bool, optional
        Whether to update template's t_c and phi_c (default True).

    Returns
    -------
    dict
        Dictionary with keys 'mismatch', 'phase', and 'ind'.
    """
    
    # Get the f_cut
    f_cut = Regular_precession(source).get_f_cut()
    # Make the frequency range
    f_range = np.arange(20, f_cut, 0.25)
    delta_f = 0.25
    
    # Noise PSD
    psd = Sn(f_range, delta_f)
    
    # Initialize the systems
    Source_init = Regular_precession(source)
    Template_init = Regular_precession(template)
    
    # Get the signals
    Source_signal = Source_init.precessing_strain(f_range, delta_f)
    Template_signal = Template_init.precessing_strain(f_range, delta_f)
    
    # Calculate the match
    match, ind, phase = optimized_match(Source_signal, Template_signal, psd, return_phase=True)
    
    # Mismatch
    mismatch = 1 - match
    
    if update_tc_phic:
        template['phi_c'] = phase
        template['t_c'] = template['t_c'] - ind * Template_signal.delta_t
    
    # return the mismatch
    return {'mismatch': mismatch, 'phase': phase, 'ind': ind}


def opt_mismatch_gammaP(rp_params, np_params, size_of_gammaP_arr):
    """
    Compute mismatch as a function of gamma_P.

    Parameters
    ----------
    rp_params : dict
        Regular precession parameters.
    np_params : dict
        Non-precessing parameters.
    size_of_gammaP_arr : int
        Number of gamma_P values to scan.

    Returns
    -------
    gammaP_arr : array_like
        Array of gamma_P values.
    mismatch_arr : array_like
        Array of mismatch values.
    ind_arr : array_like
        Array of indices of best match.
    phase_m_arr : array_like
        Array of phase values at best match.
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
    """
    Return ideal parameters for a given minimum gamma_P.

    Parameters
    ----------
    rp_params : dict
        Regular precession parameters.
    min_gammaP : float
        Value of gamma_P for which to return ideal params.

    Returns
    -------
    dict
        Updated rp_params.
    """
    rp_params["gamma_P"] = min_gammaP
    return rp_params


def opt_mismatch_extremas_gammaP(rp_params, np_params, size_of_gammaP_arr, return_tc_phic):
    """
    Find the minimum and maximum mismatch for an array of gamma_P values.

    Parameters
    ----------
    rp_params : dict
        Regular precession parameters.
    np_params : dict
        Non-precessing parameters.
    size_of_gammaP_arr : int
        Size of gamma_P array.
    return_tc_phic : bool
        Whether to return t_c and phi_c.

    Returns
    -------
    min_mismatch : float
        Minimum mismatch.
    max_mismatch : float
        Maximum mismatch.
    min_gammaP : float
        Gamma_P value at minimum mismatch.
    max_gammaP : float
        Gamma_P value at maximum mismatch.
    min_ind : int, optional
        Index of minimum mismatch (if return_tc_phic is True).
    min_phase : float, optional
        Phase at minimum mismatch (if return_tc_phic is True).
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



def get_SNRs(rp_params, np_params, redshift=100) -> dict:
    """
    Get SNRs for regular precession (RP) and non-precessing (NP) signals, and cross terms.

    Parameters
    ----------
    rp_params : dict
        Regular precession parameters.
    np_params : dict
        Non-precessing parameters.
    redshift : float, optional
        Redshift (default 100, serves as a placeholder for no redshift correction).

    Returns
    -------
    SNR_dict : dict
        Dictionary with SNR values for RP, NP, and cross terms.
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
    """
    Compute the Lindblom inequality for waveform distinguishability.

    Parameters
    ----------
    rp_params : dict
        Regular precession parameters.
    np_params : dict
        Non-precessing parameters.
    size_of_gammaP_arr : int
        Size of gamma_P array.
    rp_SNR_term : bool, optional
        Use RP SNR term (default True).
    np_SNR_term : bool, optional
        Use NP SNR term (default False).

    Returns
    -------
    lindblom : float or None
        Value of the Lindblom inequality, or None if no term selected.
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
    """
    Vary q while keeping total mass constant. Returns new parameter dicts and chirp mass.

    Parameters
    ----------
    q : float
        Mass ratio.
    mcz : float
        Chirp mass.
    rp_params : dict
        Regular precession parameters.
    np_params : dict
        Non-precessing parameters.

    Returns
    -------
    rp_params : dict
        Updated regular precession parameters.
    np_params : dict
        Updated non-precessing parameters.
    mcz_new : float
        New chirp mass.
    eta : float
        Symmetric mass ratio.
    """
    eta = q/(1+q)**2
    default_mcz = mcz
    total_mass_q1 = mcz/(0.25**0.6)
    mcz_new = total_mass_q1 * eta**0.6
    rp_params['eta'] = np_params['eta'] = eta
    rp_params['mcz'] = np_params['mcz'] = mcz_new*solar_mass
    return rp_params, np_params, mcz_new, eta


def redshift_to_luminosity_distance(z) -> float:
    """
    Convert redshift to luminosity distance in standard Lambda CDM cosmology.

    Parameters
    ----------
    z : float
        Redshift.

    Returns
    -------
    luminosity_distance : float
        Luminosity distance in Gpc.
    """
    # Define the cosmology
    cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315) # Planck 2018 results
    luminosity_distance = cosmo.luminosity_distance(z).to('Gpc').value
    return luminosity_distance



def redshifted_new_params(z, rp_params) -> dict:
    """
    Get redshifted parameters for a system.

    Parameters
    ----------
    z : float
        Redshift.
    rp_params : dict
        Regular precession parameters.

    Returns
    -------
    rp_params : dict
        Updated regular precession parameters with redshifted mass and distance.
    """
    rp_params["dist"] = redshift_to_luminosity_distance(z)*giga_parsec
    old_chirp = rp_params["mcz"]
    #new chirp mass
    rp_params["mcz"] = old_chirp*(1+z)
    return rp_params