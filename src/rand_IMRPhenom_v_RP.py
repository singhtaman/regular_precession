import sys
sys.path.append("../src")

from IMRPhenomXPHM_comparison import *

import warnings
from scipy.integrate import IntegrationWarning

warnings.filterwarnings("ignore", category=IntegrationWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", message=".*Elliptic intergal.*")

import os

z_d = 0.3

chi_1_d = 1.0
chi_2_d = 1.0

q_d = 1.0

rp_params = rp_params.copy()

rp_params["dist"] = redshift_to_luminosity_distance(z_d) * giga_parsec
rp_params["mcz"] = 20 * solar_mass

rp_params["eta"] = q_d/(1+q_d)**2

rp_params = redshifted_new_params(z_d, rp_params)

N = int(os.getenv("NumberofBinaries", 1000))

# Generate random orientations for the binaries
theta_J_arr = np.zeros(N)
phi_J_arr = np.zeros(N)

theta_S_arr = np.zeros(N)
phi_S_arr = np.zeros(N)

theta_1_arr = np.zeros(N)
phi_1_arr = np.zeros(N)

theta_2_arr = np.zeros(N)
phi_2_arr = np.zeros(N)

for i in range(N):
    number_of_random_vars = 8
    
    get_random_vars = np.random.random(number_of_random_vars)
    
    # Random J orientation
    cos_theta_J = 2 * get_random_vars[0] - 1
    theta_J = np.arccos(cos_theta_J)
    phi_J = 2 * np.pi * get_random_vars[1]
    
    # Random N orientation
    cos_theta_S = 2 * get_random_vars[2] - 1
    theta_S = np.arccos(cos_theta_S)
    phi_S = 2 * np.pi * get_random_vars[3]
    
    # Random spin orientations
    cos_theta_1 = 2 * get_random_vars[4] - 1
    theta_1 = np.arccos(cos_theta_1)
    phi_1 = 2 * np.pi * get_random_vars[5]

    cos_theta_2 = 2 * get_random_vars[6] - 1
    theta_2 = np.arccos(cos_theta_2)
    phi_2 = 2 * np.pi * get_random_vars[7]
    
    if theta_1 > np.pi or theta_2 > np.pi:
        raise ValueError("Theta values should be between 0 and pi")
    if phi_1 > 2*np.pi or phi_2 > 2*np.pi:
        raise ValueError("Phi values should be between 0 and 2*pi")
    
    theta_J_arr[i] = theta_J
    phi_J_arr[i] = phi_J
    
    theta_S_arr[i] = theta_S
    phi_S_arr[i] = phi_S
    
    theta_1_arr[i] = theta_1
    phi_1_arr[i] = phi_1
    
    theta_2_arr[i] = theta_2
    phi_2_arr[i] = phi_2
    
    
# Initialize empty lists for IMR mismatches
Mismatch_IMR_22_full_arr = []
Mismatch_IMR_22_cut_arr = []

Mismatch_IMR_2m_full_arr = []
Mismatch_IMR_2m_cut_arr = []

Mismatch_IMR_3m_full_arr = []
Mismatch_IMR_3m_cut_arr = []

Mismatch_IMR_all_full_arr = []
Mismatch_IMR_all_cut_arr = []

# Initialize empty lists for RP mismatches
Mismatch_RP_arr = []

fmax = 1024

params_package_arr = []

for theta_J, phi_J, theta_S, phi_S, theta_1, theta_2, phi_1, phi_2 in zip(theta_J_arr, phi_J_arr, theta_S_arr, phi_S_arr, theta_1_arr, theta_2_arr, phi_1_arr, phi_2_arr):
    # Update the parameters for the current binary
    rp_params["theta_J"] = theta_J
    rp_params["phi_J"] = phi_J
    
    rp_params["theta_S"] = theta_S
    rp_params["phi_S"] = phi_S
    
    rp_params["gamma_P"] = get_phi_JL0_IMRPhenomXPHM(get_IMRPhenomXPHM_dict_from_RP_params(rp_params, chi_1, chi_2, theta_1, theta_2, phi_1, phi_2)) - np.pi/2
    
    params_package_element = {
        "rp_params": rp_params.copy(),
        "theta_1": theta_1,
        "theta_2": theta_2,
        "phi_1": phi_1,
        "phi_2": phi_2
    }
    
    params_package_arr.append(params_package_element)
    
import multiprocessing as mp

def process_each_binary(param_package):
    # Unpack the parameters for the current binary
    params = param_package["rp_params"]
    theta_1 = param_package["theta_1"]
    theta_2 = param_package["theta_2"]
    phi_1 = param_package["phi_1"]
    phi_2 = param_package["phi_2"]
    
    
    # Now given the updated parameters, we can compute the precessing RP strain
    freq_range_RP, strain_RP = get_RP_Waveform(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, return_freqs=True)
    
    #we should put a check here to make sure that IMRPhenomXPHM can generate a waveform for these parameters, and if not, we can skip this binary and move on to the next one. This is important because there are some regions of the parameter space where IMRPhenomXPHM fails to generate waveforms, and we don't want our code to break because of that.
    
    try:
        # Compute the IMRPhenomXPHM strain for the same parameters with just the (2, 2) mode
        freq_range_IMR_22, strain_IMR_22 = get_hstrain_from_modes(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, psi=0, return_freqs=True, mode_list=[(2, -2), (2, 2)])
    except:
        print(f"Failed to generate IMRPhenomXPHM waveform for binary with parameters: {params} and spin angles: {theta_1}, {theta_2}, {phi_1}, {phi_2}. Skipping this binary.")
        return None, None, None
    
    # Compute the IMRPhenomXPHM strain for the same parameters with just the (2, m) modes
    freq_range_IMR_2m, strain_IMR_2m = get_hstrain_from_modes(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, psi=0, return_freqs=True, mode_list=[(2, -2), (2, -1), (2, 0), (2, 1), (2, 2)])
    
    # Compute the IMRPhenomXPHM strain for the same parameters with just the (2, m) and (3, m) modes
    freq_range_IMR_3m, strain_IMR_3m = get_hstrain_from_modes(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, psi=0, return_freqs=True, mode_list=[(2, -2), (2, -1), (2, 0), (2, 1), (2, 2), (3, -3), (3, -2), (3, -1), (3, 0), (3, 1), (3, 2), (3, 3)])
    
    # Compute the IMRPhenomXPHM strain for the same parameters with all modes
    freq_range_IMR_all, strain_IMR_all = get_hstrain_from_modes(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, psi=0, return_freqs=True, mode_list=None)
    
    # We also need non-precessing waveforms to compute the mismatch with the non-precessing limit of the RP model, so we can set the precession parameters to zero and compute the non-precessing RP waveform
    rp_params_non_prec = params.copy()
    chi_1_non_prec = 0.0
    chi_2_non_prec = 0.0
    
    freq_range_RP_non_prec, strain_RP_non_prec = get_RP_Waveform(rp_params_non_prec, chi_1_non_prec, chi_2_non_prec, theta_1, theta_2, phi_1, phi_2, fmax, return_freqs=True)
    
    # Lets also compute the non-precessing IMRPhenomXPHM waveform for the same parameters (reminder: The non-precessing limit has aligned spins - not necessarily zero spins - so we can set the spin angles to zero but keep the spin magnitudes the same)
    
    freq_range_IMR_22_non_prec, strain_IMR_22_non_prec = get_hstrain_NP_align_from_modes(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, psi=0, return_freqs=True, mode_list=[(2, -2), (2, 2)])
    
    freq_range_IMR_2m_non_prec, strain_IMR_2m_non_prec = get_hstrain_NP_align_from_modes(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, psi=0, return_freqs=True, mode_list=[(2, -2), (2, -1), (2, 0), (2, 1), (2, 2)])
    
    freq_range_IMR_3m_non_prec, strain_IMR_3m_non_prec = get_hstrain_NP_align_from_modes(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, psi=0, return_freqs=True, mode_list=[(2, -2), (2, -1), (2, 0), (2, 1), (2, 2), (3, -3), (3, -2), (3, -1), (3, 0), (3, 1), (3, 2), (3, 3)])
    
    freq_range_IMR_all_non_prec, strain_IMR_all_non_prec = get_hstrain_NP_align_from_modes(params, chi_1_d, chi_2_d, theta_1, theta_2, phi_1, phi_2, fmax, psi=0, return_freqs=True, mode_list=None)
    
    # Now we can compute the mismatches between the RP and NP RP waveforms, and the IMRPhenomXPHM waveforms with different mode content and their non-precessing limits
    
    mismatch_RP_NP_RegPre = mismatch_IMRPhenomXPHM(strain_RP, strain_RP_non_prec, freq_range_RP, freq_range_RP_non_prec, f_cut = freq_range_RP[-1])
    
    mismatch_RP_NP_IMR_22_full = mismatch_IMRPhenomXPHM(strain_IMR_22, strain_IMR_22_non_prec, freq_range_IMR_22, freq_range_IMR_22_non_prec, f_cut = freq_range_IMR_22[-1])
    
    mismatch_RP_NP_IMR_22_cut = mismatch_IMRPhenomXPHM(strain_IMR_22, strain_IMR_22_non_prec, freq_range_IMR_22, freq_range_IMR_22_non_prec, f_cut = freq_range_RP[-1])
    
    mismatch_RP_NP_IMR_2m_full = mismatch_IMRPhenomXPHM(strain_IMR_2m, strain_IMR_2m_non_prec, freq_range_IMR_2m, freq_range_IMR_2m_non_prec, f_cut = freq_range_IMR_2m[-1])
    
    mismatch_RP_NP_IMR_2m_cut = mismatch_IMRPhenomXPHM(strain_IMR_2m, strain_IMR_2m_non_prec, freq_range_IMR_2m, freq_range_IMR_2m_non_prec, f_cut = freq_range_RP[-1])
    
    mismatch_RP_NP_IMR_3m_full = mismatch_IMRPhenomXPHM(strain_IMR_3m, strain_IMR_3m_non_prec, freq_range_IMR_3m, freq_range_IMR_3m_non_prec, f_cut = freq_range_IMR_3m[-1])
    
    mismatch_RP_NP_IMR_3m_cut = mismatch_IMRPhenomXPHM(strain_IMR_3m, strain_IMR_3m_non_prec, freq_range_IMR_3m, freq_range_IMR_3m_non_prec, f_cut = freq_range_RP[-1])
    
    mismatch_RP_NP_IMR_all_full = mismatch_IMRPhenomXPHM(strain_IMR_all, strain_IMR_all_non_prec, freq_range_IMR_all, freq_range_IMR_all_non_prec, f_cut = freq_range_IMR_all[-1])
    
    mismatch_RP_NP_IMR_all_cut = mismatch_IMRPhenomXPHM(strain_IMR_all, strain_IMR_all_non_prec, freq_range_IMR_all, freq_range_IMR_all_non_prec, f_cut = freq_range_RP[-1])
    
    mismatch_full_dict = {
        "2, 2": mismatch_RP_NP_IMR_22_full,
        "2, m": mismatch_RP_NP_IMR_2m_full,
        "3, m": mismatch_RP_NP_IMR_3m_full,
        "all": mismatch_RP_NP_IMR_all_full
    }
    
    mismatch_cut_dict = {
        "2, 2": mismatch_RP_NP_IMR_22_cut,
        "2, m": mismatch_RP_NP_IMR_2m_cut,
        "3, m": mismatch_RP_NP_IMR_3m_cut,
        "all": mismatch_RP_NP_IMR_all_cut
    }
    
    return mismatch_RP_NP_RegPre, mismatch_full_dict, mismatch_cut_dict


# Initialize empty lists
Mismatch_compare_RP_arr = []
Mismatch_compare_IMR_full_arr = []
Mismatch_compare_IMR_cut_arr = []

if __name__ == "__main__":
    pool = mp.Pool(mp.cpu_count())
    
    results = pool.map(process_each_binary, params_package_arr)
    
    pool.close()

    Mismatch_compare_RP_arr = [result[0] for result in results]
    Mismatch_compare_IMR_full_arr = [result[1] for result in results]
    Mismatch_compare_IMR_cut_arr = [result[2] for result in results]
    
    
Dict_to_save_for_RP_IMR_comparison = {
    "Parameters": params_package_arr,
    "Mismatch RP vs NP RegPre": Mismatch_compare_RP_arr,
    "Mismatch RP vs NP IMR full": Mismatch_compare_IMR_full_arr,
    "Mismatch RP vs NP IMR cut": Mismatch_compare_IMR_cut_arr
}

import pickle

with open("saved_data/Comparison_RP_IMRPhenomXPHM_mcz20_q1_chi1_1_chi2_1_z03.pkl", "wb") as f:
    pickle.dump(Dict_to_save_for_RP_IMR_comparison, f)