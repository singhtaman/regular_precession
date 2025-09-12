import numpy as np

import precession
import sys

error_handler = np.seterr(invalid="raise")

solar_mass = 4.92624076 * 1e-6              # [solar_mass] = sec
giga_parsec = 1.02927125 * 1e17             # [giga_parsec] = sec
year = 31557600                             # [year] = sec

# Importing the regular precession class
sys.path.insert(0, "../")
from src.regular_precession import *

# Importing system parameters
from src.systems_lib import *

# Importing the mismatch_SNR script functions
from src.mismatch_n_SNR import *


def random_angles_generator():
    """
    Generate random angles for J and S vectors.

    Returns
    -------
    thetaJ : float
        Polar angle for J.
    phiJ : float
        Azimuthal angle for J.
    thetaS : float
        Polar angle for S.
    phiS : float
        Azimuthal angle for S.
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


def get_param_distributions(q, s, N, r):
    """
    Get precession parameters for a population of binaries.

    Parameters
    ----------
    q : float
        Mass ratio.
    s : float
        Spin (equal spin case).
    N : int
        Number of binaries in the population.
    r : float
        Binary separation.

    Returns
    -------
    bracket_theta_arr : array_like
        Array of bracket theta values.
    bracket_omega_arr : array_like
        Array of bracket omega values.
    delta_theta_arr : array_like
        Array of delta theta values.
    delta_omega_arr : array_like
        Array of delta omega values.
    little_omega_arr : array_like
        Array of little omega values.
    """
   
    vec = precession.isotropic_angles(N)

    # N binaries with random theta1, theta2 and deltaphi
    theta1_arr, theta2_arr, deltaphi_arr = vec[0], vec[1], vec[2]

    #initializing all the arrays
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


def random_precession_params_iso(chi1=1, chi2=1, q=1, r=6):
    """
    Generate random precession parameters for isotropic spins.

    Parameters
    ----------
    chi1 : float, optional
        Spin of first object (default 1).
    chi2 : float, optional
        Spin of second object (default 1).
    q : float, optional
        Mass ratio (default 1).
    r : float, optional
        Binary separation (default 6).

    Returns
    -------
    theta_tilde : float
        Dimensionless theta.
    omega_tilde : float
        Dimensionless omega.
    bracket_theta : float
        Bracket theta.
    bracket_omega : float
        Bracket omega.
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
    omega_tilde = bracket_omega/(1000*solar_mass)

    return theta_tilde, omega_tilde, bracket_theta, bracket_omega


def random_angles_generator_along_w_prec(arg_for_distribution):
    """
    Generate random angles for J, S, and optionally spin angles for different distributions.

    Parameters
    ----------
    arg_for_distribution : str or None
        Distribution type ('none', 'ISO', 'WA', 'SA').

    Returns
    -------
    tuple
        Angles for J, S, and possibly spin angles depending on distribution.
    """
    n = 7
    get_random_numbers = np.random.random(n)
    
    # Random angles for J
    cos_theta_J = 2*get_random_numbers[0] - 1
    thetaJ = np.arccos(cos_theta_J)
    phiJ = 2*np.pi*get_random_numbers[1]
    
    # Random angles for S
    cos_theta_S = 2*get_random_numbers[2] - 1
    thetaS = np.arccos(cos_theta_S)
    phiS = 2*np.pi*get_random_numbers[3]
    
    if arg_for_distribution == 'none' or arg_for_distribution is None:
        
        return thetaJ, phiJ, thetaS, phiS
    
    elif arg_for_distribution == 'ISO':
        cos_theta_1 = 2*get_random_numbers[4] - 1
        theta_1 = np.arccos(cos_theta_1)
        
        cos_theta_2 = 2*get_random_numbers[5] - 1
        theta_2 = np.arccos(cos_theta_2)
        
        deltaphi = 2*np.pi*get_random_numbers[6]
        
        return thetaJ, phiJ, thetaS, phiS, theta_1, theta_2, deltaphi
    
    elif arg_for_distribution == 'WA':
        cos_theta_1_inf = (1 - np.cos(np.pi/6))*get_random_numbers[4] + np.cos(np.pi/6)
        theta_1_inf = np.arccos(cos_theta_1_inf)
        
        cos_theta_2_inf = (1 - np.cos(np.pi/6))*get_random_numbers[5] + np.cos(np.pi/6)
        theta_2_inf = np.arccos(cos_theta_2_inf)
        
        deltaphi = 2*np.pi*get_random_numbers[6]
        
        return thetaJ, phiJ, thetaS, phiS, theta_1_inf, theta_2_inf, deltaphi
    
    elif arg_for_distribution == 'SA':
        cos_theta_1_inf = (1 - np.cos(np.pi/18))*get_random_numbers[4] + np.cos(np.pi/18)
        theta_1_inf = np.arccos(cos_theta_1_inf)
        
        cos_theta_2_inf = (1 - np.cos(np.pi/18))*get_random_numbers[5] + np.cos(np.pi/18)
        theta_2_inf = np.arccos(cos_theta_2_inf)
        
        deltaphi = 2*np.pi*get_random_numbers[6]
        
        return thetaJ, phiJ, thetaS, phiS, theta_1_inf, theta_2_inf, deltaphi
    
    else:
        print("Please enter a valid argument for the distribution")
        return None

def prec_av_inf_to_f20_rp_params(theta_1_inf, theta_2_inf, delta_phi, q, chi1, chi2, r20, binary_spin_orientation):
    """
    Precessional average inspiral from infinity to f = 20 Hz.

    Parameters
    ----------
    theta_1_inf : float
        Theta1 at infinity.
    theta_2_inf : float
        Theta2 at infinity.
    delta_phi : float
        Delta phi.
    q : float
        Mass ratio.
    chi1 : float
        Spin of first object.
    chi2 : float
        Spin of second object.
    r20 : float
        Binary separation at f = 20 Hz.
    binary_spin_orientation : str
        Spin orientation ('WA', 'SA', etc.).

    Returns
    -------
    theta_tilde : float
        Dimensionless precession amplitude.
    omega_tilde : float
        Dimensionless precession frequency.
    """
    if binary_spin_orientation == 'ISO':
        print('Please call the function with the correct spin orientation: WA or SA')
        print('If you need to calculate precession parameters for isotropic spin orientation, please call the function: get_prec_params_t1_t2_dphi with the correct arguments')
        return warnings.warn('Spin orientation is isotropic, dont need to calculate the precessional averaged values, wrong function call')
    
    elif binary_spin_orientation == 'none' or binary_spin_orientation is None:
        return warnings.warn('Spin orientation is not specified, please specify the spin orientation, so that precessional averaged values can be calculated, if needed.')
    
    else:
        u20 = precession.eval_u(r20, q)
        u_arr = np.linspace(0., u20, 101)
        
        inspiral_prec_av = precession.inspiral_precav(theta1= theta_1_inf, theta2=theta_2_inf, deltaphi=delta_phi, u=u_arr[100], chi1=chi1, chi2=chi2, q=q, requested_outputs=['theta1','theta2','deltaphi'])
        
        theta_1_r20 = inspiral_prec_av['theta1']
        theta_2_r20 = inspiral_prec_av['theta2']
        delta_phi_r20 = inspiral_prec_av['deltaphi']
        
        kappa_r20_arr = precession.eval_kappa(theta1=theta_1_r20, theta2=theta_2_r20, deltaphi=delta_phi_r20, r=r20, q=q, chi1=chi1, chi2=chi2)
        kappa_r20 = kappa_r20_arr.flatten()
        
        chieff_r20_arr = precession.eval_chieff(theta_1_r20, theta_2_r20, q, chi1, chi2)
        chieff_r20 = chieff_r20_arr.flatten()
        
        bracket_theta = precession.eval_bracket_theta(kappa_r20, r20, chieff_r20, q, chi1, chi2)
        bracket_omega = precession.eval_bracket_omega(kappa_r20, r20, chieff_r20, q, chi1, chi2)
        
        new_bracket_theta_r20 = bracket_theta * ((r20/6.0)**(1.0/2.0))
        new_bracket_omega_r20 = bracket_omega * ((r20/6.0)**(5.0/2.0))
        
        theta_tilde = new_bracket_theta_r20*10
        omega_tilde = new_bracket_omega_r20/(1000*solar_mass)
        
        return theta_tilde, omega_tilde



def get_prec_params_t1_t2_dphi(theta_1, theta_2, delta_phi, q, r, chi1, chi2):
    """
    Get precession parameters for given spin angles and binary separation.

    Parameters
    ----------
    theta_1 : float
        Spin angle 1.
    theta_2 : float
        Spin angle 2.
    delta_phi : float
        Delta phi.
    q : float
        Mass ratio.
    r : float
        Binary separation.
    chi1 : float
        Spin of first object.
    chi2 : float
        Spin of second object.

    Returns
    -------
    theta_tilde : float
        Dimensionless precession amplitude.
    omega_tilde : float
        Dimensionless precession frequency.
    bracket_theta : float
        Precession amplitude.
    bracket_omega : float
        Precession frequency.
    """
    chieff = precession.eval_chieff(theta_1, theta_2, q, chi1, chi2)
    kappa = (chi1 * np.cos(theta_1) + q**2 * chi2 * np.cos(theta_2) )/(1+q)**2 + \
            (chi1**2 + q**4 *chi2**2 + 2*chi1*chi2*q**2 * (np.cos(theta_1)*np.cos(theta_2) + np.cos(delta_phi)*np.sin(theta_1)*np.sin(theta_2))) / (2*q*(1+q)**2*r**(1/2))
    
    bracket_theta = precession.eval_bracket_theta(kappa, r, chieff, q, chi1, chi2)
    bracket_omega = precession.eval_bracket_omega(kappa, r, chieff, q, chi1, chi2)
    
    if r == 6.0:
        theta_tilde = bracket_theta*10
        omega_tilde = bracket_omega/(1000*solar_mass)
        
    else:
        new_bracket_theta = bracket_theta * ((r/6.0)**(1.0/2.0))
        new_bracket_omega = bracket_omega * ((r/6.0)**(5.0/2.0))
        
        theta_tilde = new_bracket_theta*10
        omega_tilde = new_bracket_omega/(1000*solar_mass)
        
    
    return theta_tilde, omega_tilde