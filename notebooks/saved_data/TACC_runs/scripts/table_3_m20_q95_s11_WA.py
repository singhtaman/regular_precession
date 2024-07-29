"""
This is a script to calculate the mismatch between two waveforms with random parameters.
"""
import sys

sys.path.insert(0, "../../../../")


from mismatch_n_SNR import *
from MC_skyloc_binorient_plus_prec import *
import numpy as np
import matplotlib.pyplot as plt


#default parameters
mcz_d = 20 * solar_mass
q_d = 0.95
z_d = 0.1

mcz_arr = np.array([10, 20, 40])*solar_mass
q_arr = np.array([0.1, 0.5, 1])
z_arr = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1])

chi1_d = 1.0
chi2_d = 1.0

chi1_arr = np.array([0.1, 0.2, 0.5, 1.0])
chi2_arr = np.array([0.1, 0.2, 0.5, 1.0])


rp_params = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 2,
    "phi_J": np.pi / 2,
    "mcz": mcz_d,
    "dist": redshift_to_luminosity_distance(z_d) * giga_parsec,
    "eta": q_d/(1+q_d)**2,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 4.0,  # dimensionless
    "omega_tilde": 2.0,  # dimensionless
    "gamma_P": 0.0
}

np_params = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 2,
    "phi_J": np.pi / 2,
    "mcz": mcz_d,
    "dist": redshift_to_luminosity_distance(z_d) * giga_parsec,
    "eta": q_d/(1+q_d)**2,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 0.0,
    "omega_tilde": 0.0,
    "gamma_P": 0.0
}


binary_spin_dist = 'WA'

# Total number of binary orientations
N = 10000

theta_rand_J_arr, phi_rand_J_arr, theta_rand_S_arr, phi_rand_S_arr, theta_1_rand_S_arr, theta_2_rand_S_arr, delta_phi_rand_S_arr = zip(*[random_angles_generator_along_w_prec(binary_spin_dist) for _ in range(N)])

rp_param_JnS_arr = []
np_param_JnS_arr = []

f_ref = 20.0

eta_d = q_d/(1+q_d)**2

import multiprocessing as mp

total_mass = (eta_d**(-0.6))*mcz_d

r20 = ((np.pi*f_ref*total_mass)**(-2/3))

for theta_J, phi_J, theta_S, phi_S, theta_1, theta_2, delta_phi in zip(theta_rand_J_arr, phi_rand_J_arr, theta_rand_S_arr, phi_rand_S_arr, theta_1_rand_S_arr, theta_2_rand_S_arr, delta_phi_rand_S_arr):
    rp_params["theta_J"] = theta_J
    rp_params["phi_J"] = phi_J
    
    rp_params["theta_S"] = theta_S
    rp_params["phi_S"] = phi_S
    
    r20 = ((np.pi*f_ref*total_mass)**(-2/3))
    
    if binary_spin_dist == 'ISO':
        theta_tilde_new, omega_tilde_new = get_prec_params_t1_t2_dphi(theta_1, theta_2, delta_phi, q_d, r20, chi1=chi1_d, chi2=chi2_d)
        
    elif binary_spin_dist == 'WA' or binary_spin_dist == 'SA':
        theta_tilde_new, omega_tilde_new = prec_av_inf_to_f20_rp_params(theta_1, theta_2, delta_phi, q_d, chi1_d, chi2_d, r20=r20, binary_spin_orientation=binary_spin_dist)
    
    rp_params["theta_tilde"] = theta_tilde_new
    rp_params["omega_tilde"] = omega_tilde_new
    
    np_params = rp_params.copy()
    np_params["theta_tilde"] = 0.0
    np_params["omega_tilde"] = 0.0
    
    rp_param_JnS_arr.append(rp_params.copy())
    np_param_JnS_arr.append(np_params.copy())



def process_rp_params(rp_params):
    
    rp_params['eta'] = q_arr[2]/(1+q_arr[2])**2
        
    SNR_dict_z_arr = []
    mismatch_dict_z_arr = []
    
    for z in z_arr:
        rp_params['mcz'] = mcz_arr[1]
        
        rp_params = redshifted_new_params(z, rp_params)
        
        np_params = rp_params.copy()
        np_params["theta_tilde"] = 0.0
        np_params["omega_tilde"] = 0.0
        
        min_mismatch, max_mismatch, min_gammaP, max_gammaP, min_ind, min_phase = opt_mismatch_extremas_gammaP(rp_params, np_params, size_of_gammaP_arr = 101, return_tc_phic = True)
        
        mismatch_record = {
            "min_mismatch": min_mismatch,
            "max_mismatch": max_mismatch,
            "min_gammaP": min_gammaP,
            "max_gammaP": max_gammaP,
            "min_ind": min_ind,
            "min_phase": min_phase,
        }
        
        rp_params["gamma_P"] = min_gammaP
        
        SNR_dict = get_SNRs(rp_params, np_params)
        
        SNR_dict_z_arr.append(SNR_dict)
        mismatch_dict_z_arr.append(mismatch_record)

    return SNR_dict_z_arr, mismatch_dict_z_arr

# Initialize empty lists
SNR_dict_JnS_arr = []
Mismatch_dict_JnS_arr = []

if __name__ == "__main__":
    pool = mp.Pool(mp.cpu_count())

    results = pool.map(process_rp_params, rp_param_JnS_arr)
    
    pool.close()

    SNR_dict_JnS_arr = [result[0] for result in results]
    Mismatch_dict_JnS_arr = [result[1] for result in results]
    
"""
# Parallelize processing for each rp_params
with mp.Pool() as pool:
    results = pool.map(process_rp_params, rp_param_JnS_arr)
    for SNR_dict_z_arr, mismatch_dict_z_arr in results:
        SNR_dict_JnS_arr.append(SNR_dict_z_arr)
        Mismatch_dict_JnS_arr.append(mismatch_dict_z_arr)
"""
    
DICT_Par_Mis_SNR = {
    "rp_param_JnS_arr": rp_param_JnS_arr,
    "np_param_JnS_arr": np_param_JnS_arr,
    "Mismatch_dict_JnS_arr": Mismatch_dict_JnS_arr,
    "SNR_dict_JnS_arr": SNR_dict_JnS_arr
}

# Need to save the data in a file

import pickle

#mass info

with open("../table_3/WA/DICT__param_mismatch_SNR_for_m20_q95_s11_WA_parallel_new.pkl", "wb") as f: # _new so it doesn't overwrite the previously saved file
    pickle.dump(DICT_Par_Mis_SNR, f)

print("Done") # while running an interactive shell, a print statement to say that its done