""" Regular Precession: When L moves in a cone around J with an opening angle theta_tilde that changes on a radiation reaction timescale, 
                        with frequency omega_tilde (also changing on the same timescale) and a phase gamma_P.
                        
Model presented in following paper : *Add arXiv number here or DOI identifier*
"""

import numpy as np
from scipy.integrate import odeint
from pycbc.types import FrequencySeries

error_handler = np.seterr(invalid="raise")

#suppressing a warning - UserWarning: Wswiglal-redir-stdio
import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
import lal
import lal as _lal

# Regular precessing class

class Regular_precession():

    def __init__(self, params) -> None:
        
        self.params = params
        # non-precessiing parameters
        self.theta_S = params['theta_S']            # Sky inclination -- polar angle for the line of sight vector in detector frame
        self.phi_S = params['phi_S']                # Sky azimuthal angle -- azimuthal angle for the line of sight vector in detector frame
        self.theta_J = params['theta_J']            # Source binary plane inclination -- polar angle for the total angular momentum vector in detector frame
        self.phi_J = params['phi_J']                # Source binary plane azimuthal angle -- azimuthal angle for the total angular momentum vector in detector frame
        self.mcz = params['mcz']                    # chirp mass [s] -- chirp mass in seconds : M_c = (m1*m2)**(3/5) / (m1 + m2)**(1/5)
        self.dist = params['dist']                  # distance to the source -- distance to the source (usually in Gpc)
        self.eta = params['eta']                    # symmetric mass ratio -- eta = m1*m2 / (m1 + m2)**2
        self.t_c = params['t_c']                    # coalescence time
        self.phi_c = params['phi_c']                # coalescence phase

        # regular precession parameters
        self.theta_tilde = params['theta_tilde']    # Dimensionless precession amplitude -- stands for 10 times the opening angle of the cone at binary separation r = 6M        
        self.omega_tilde = params['omega_tilde']    # Dimensionless precession frequency -- stands for 1000 times the precession frequency at binary separation r = 6M for a solar mass binary 
        self.gamma_P = params['gamma_P']            # Phase of the precession -- stands for the phase of the precession when the binary enters the detector band

        # some converters/constants
        
        self.SOLMASS2SEC = 4.92624076 * 1e-6        # solar mass -> seconds
        self.GIGAPC2SEC = 1.02927125 * 1e17         # gigaparsec -> seconds
        self.FMIN = 20                              # lower frequency of the detector sensitivity band [Hz]

    def get_total_mass(self):
        """ For the total mass of the binary system from chirp mass [seconds]
        ___________________________________________________________________________________________
        Parameters used:
        mcz: chirp mass [s]
        eta: symmetric mass ratio
        ___________________________________________________________________________________________
        Returns:
        total mass of the binary system [s]
        ___________________________________________________________________________________________
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        return self.mcz/(self.eta**(3/5))           # total mass of the binary system [s]: M = M_c / eta**(3/5)

    def get_f_cut(self):
        """ Cut-off frequency: frequency where the binary coalesces
        From equation 13
        ___________________________________________________________________________________________
        Parameters used:
        total mass of the binary system [s] -- get_total_mass() function : requires chirp mass and symmetric mass ratio
        ___________________________________________________________________________________________
        Returns:
        f_cut: cut-off frequency [Hz]
        ___________________________________________________________________________________________
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        return 1/(6**(3/2)*np.pi*self.get_total_mass()) # Equation 13 -- cut-off frequency [Hz]: f_cut = 1/(r_{ISCO}**(3/2) * pi * M)
    

    def get_theta_LJ(self, f):
        """ Angle between L and J at a given frequency [rad] -- in the L-J plane
        From equation 18a
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency where the angle is to be calculated [Hz]
        theta_tilde: Dimensionless precession amplitude -- stands for 10 times the opening angle of the cone at binary separation r = 6M
        ___________________________________________________________________________________________
        Returns:
        theta_LJ: Angle between L and J at a given frequency [rad]
        ___________________________________________________________________________________________
        get_f_cut() function is used to get the cut-off frequency
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        return 0.1*self.theta_tilde*(f/self.get_f_cut())**(1/3) # Equation 18a  -- Angle between L and J at a given frequency [rad]: theta_LJ = 0.1 * theta_tilde * (f/f_cut)**(1/3)
    
    def get_phi_LJ(self, f):
        """ Angle between projection of L in the x-y plane and x axis in source frame at a given frequency [rad] -- in the plane perpendicular to J
        From equation 18b and 19
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency where the angle is to be calculated [Hz]
        omega_tilde: Dimensionless precession frequency -- stands for 1000 times the precession frequency at binary separation r = 6M for a solar mass binary
        gamma_P: Phase of the precession -- stands for the phase of the precession when the binary enters the detector band
        ___________________________________________________________________________________________
        Returns:
        phi_LJ: Angle between projection of L in the x-y plane and x axis in source frame at a given frequency [rad]
        ___________________________________________________________________________________________
        get_f_cut() function is used to get the cut-off frequency
        get_total_mass() function is used to get the total mass of the binary system
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ***Keep an eye on this function for non-precessing case, this might be an issue***
        """
        phi_LJ_amp = (5000 * self.omega_tilde) / (96 * (self.get_total_mass()/self.SOLMASS2SEC) * (np.pi**(8/3)) * (self.mcz**(5/3)) * (self.get_f_cut()**(5/3)))
        return phi_LJ_amp * (1/self.FMIN - 1/f) + self.gamma_P # Equation 19 (also uses equation 18b) -- Angle between projection of L in the x-y plane and x axis in source frame at a given frequency [rad]: phi_LJ = phi_LJ_amp * (1/f_min - 1/f) + gamma_P
        
        
    def amp_prefactor(self) -> float:
        """ Amplitude prefactor calculated using chirp mass and distance
        From equation 6
        ___________________________________________________________________________________________
        Parameters used:
        mcz: chirp mass [s]
        dist: distance to the source (usually in Gpc) in units of seconds
        ___________________________________________________________________________________________
        Returns:
        amp_prefactor: amplitude prefactor
        ___________________________________________________________________________________________
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        """
        amp_prefactor = np.sqrt(5/96) * (np.pi**(-2/3)) * (self.mcz**(5/6)) / self.dist # from equation 6 -- A = sqrt(5/96) * (pi**(-2/3)) * (M_c**(5/6)) / D_L
        return amp_prefactor

    def precession_angles(self):
        """ Angles important for the precession model -- angles between the total angular momentum and the line of sight, and the angles between certain axes in different frames
        
        i_JN: angle between the total angular momentum and the line of sight
        o_XH: angle between the x-axis of the source frame and H in the detector frame -- H = V X N where V is the unit vector in the direction of zenith and N is the line of sight vector
        
        From equation A4, A6a, and A6b
        ___________________________________________________________________________________________
        Parameters used:
        theta_J: binary plane inclination -- polar angle for the total angular momentum vector in detector frame
        phi_J: binary plane azimuthal angle -- azimuthal angle for the total angular momentum vector in detector frame
        theta_S: source inclination -- polar angle for the line of sight vector in detector frame
        phi_S: source azimuthal angle -- azimuthal angle for the line of sight vector in detector frame
        ___________________________________________________________________________________________
        Returns:
        cos_i_JN: cosine of the angle between the total angular momentum and the line of sight
        sin_i_JN: sine of the angle between the total angular momentum and the line of sight
        cos_o_XH: cosine of the angle between the x-axis of the source frame and H in the detector frame
        sin_o_XH: sine of the angle between the x-axis of the source frame and H in the detector frame
        ___________________________________________________________________________________________
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        i_JN : Inclination angle
        o_XH : Longitude of ascending node 
        """
        cos_i_JN = np.sin(self.theta_J) * np.sin(self.theta_S) * np.cos(self.phi_J - self.phi_S) + np.cos(self.theta_J) * np.cos(self.theta_S) # Equation A4 -- cos(i_JN) = sin(theta_J) * sin(theta_S) * cos(phi_J - phi_S) + cos(theta_J) * cos(theta_S)
        sin_i_JN = np.sqrt(1 - cos_i_JN ** 2)
        
        if sin_i_JN == 0:
            cos_o_XH = 1
            sin_o_XH = 0
            """tan_o_XH = (np.sin(theta_J) * np.sin(phi_J - phi_S)) / (np.cos(theta_J) * np.sin(theta_S) * np.cos(phi_J - phi_S) + np.sin(theta_J) * np.cos(theta_S))
            cos_o_XH = 1 / np.sqrt(1 + tan_o_XH ** 2)
            sin_o_XH = np.sqrt(1 - cos_o_XH ** 2)"""
        else:
            cos_o_XH = (np.cos(self.theta_S) * np.sin(self.theta_J) * np.cos(self.phi_J - self.phi_S) - np.sin(self.theta_S) * np.cos(self.theta_J)) / (sin_i_JN) # Equation A6b cos(Omega_{XH}) = (cos(theta_S) * sin(theta_J) * cos(phi_J - phi_S) - sin(theta_S) * cos(theta_J)) / sin(i_JN)
            sin_o_XH = (np.sin(self.theta_J) * np.sin(self.phi_J - self.phi_S)) / (sin_i_JN) # Equation A6a sin(Omega_{XH}) = (sin(theta_J) * sin(phi_J - phi_S)) / sin(i_JN)
        return cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH
    
    def LdotN(self, f):
        """ Cosine of the angle between L and N
        From equation A10 (also A8)
        ___________________________________________________________________________________________
        Parameters used:
        cos_i_JN: cosine of the angle between the total angular momentum and the line of sight : Requires theta_J, phi_J, theta_S, phi_S
        sin_i_JN: sine of the angle between the total angular momentum and the line of sight : Requires theta_J, phi_J, theta_S, phi_S
        ___________________________________________________________________________________________
        Returns:
        LdotN: Dot product of L and N: Dot product of the angular momentum vector and the line of sight vector
        ___________________________________________________________________________________________
        get_theta_LJ(f) function is used to get the angle between L and J at a given frequency
        get_phi_LJ(f) function is used to get the angle between projection of L in the x-y plane and x axis in source frame at a given frequency
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = self.precession_angles()
        LdotN = np.sin(self.get_theta_LJ(f)) * sin_i_JN * np.sin(self.get_phi_LJ(f)) + np.cos(self.get_theta_LJ(f)) * cos_i_JN  
        return LdotN

    def beam_pattern_amplitude_and_phase(self, f):
        """
        Beam pattern functions for + and x polarizations in the form described in equations 2a and 2b.
        From equation 3, 4a, and 4b
        ___________________________________________________________________________________________
        Parameters used:
        theta_S: Sky inclination -- polar angle for the line of sight vector in detector frame
        phi_S: Sky azimuthal angle -- azimuthal angle for the line of sight vector in detector frame
        theta_J: source binary plane inclination -- polar angle for the total angular momentum vector in detector frame
        phi_J: source binary plane azimuthal angle -- azimuthal angle for the total angular momentum vector in detector frame
        ___________________________________________________________________________________________
        Returns:
        C_amp: amplitude of the beam pattern functions for + and x polarizations
        sin_2pa: sine of 2 times the polarization angle + alpha -- for x polarization
        cos_2pa: cosine of 2 times the polarization angle + alpha -- for + polarization
        ___________________________________________________________________________________________
        get_theta_LJ(f) function is used to get the angle between L and J at a given frequency
        get_phi_LJ(f) function is used to get the angle between projection of L in the x-y plane and x axis in source frame at a given frequency
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = self.precession_angles()
        
        # for C
        C_amp = np.sqrt((0.25 * ((1 + (np.cos(self.theta_S))**2)**2) * ((np.cos(2 * self.phi_S))**2)) + ((np.cos(self.theta_S))**2 * (np.sin(2 * self.phi_S))**2)) # Equation 4a -- C = sqrt(0.25 * (1 + cos(theta_S)**2)**2 * (cos(2 * phi_S)**2) + (cos(theta_S)**2 * sin(2 * phi_S)**2))

        # define alpha based on equation 4b
        sin_alpha = np.cos(self.theta_S) * np.sin(2 * self.phi_S) / C_amp
        cos_alpha = (1 + np.cos(self.theta_S)**2) * np.cos(2 * self.phi_S) / (2 * C_amp)
        
        # define tan_psi from equation 3
        num_psi = np.sin(self.get_theta_LJ(f)) * (np.cos(self.get_phi_LJ(f)) * sin_o_XH + np.sin(self.get_phi_LJ(f)) * cos_i_JN * cos_o_XH ) - np.cos(self.get_theta_LJ(f)) * sin_i_JN * cos_o_XH
        den_psi = np.sin(self.get_theta_LJ(f)) * (np.cos(self.get_phi_LJ(f)) * cos_o_XH - np.sin(self.get_phi_LJ(f)) * cos_i_JN * sin_o_XH ) + np.cos(self.get_theta_LJ(f)) * sin_i_JN * sin_o_XH
        
        if abs(cos_i_JN) == 1:
            o_XH = np.arctan2(sin_o_XH, cos_o_XH)
            tan_psi = np.tan(o_XH + cos_i_JN*self.get_phi_LJ(f)) # making sure it works for the face-on case
        else:
            tan_psi = num_psi / den_psi
        
        psi = np.arctan(tan_psi)
        
        #define  2 * psi + alpha
        if (2*psi%np.pi).any() == 0:
            alpha = np.arctan2(sin_alpha, cos_alpha)
            sin_2pa = np.sin(2*psi + alpha)
            cos_2pa = np.cos(2*psi + alpha)
        else:
            sin_2pa = (2 * cos_alpha * tan_psi + sin_alpha * (1 - (tan_psi)**2)) / (1 + (tan_psi)**2) # Combining equations 3 and 4b -- sin(2 * psi + alpha) = (2 * cos(alpha) * tan(psi) + sin(alpha) * (1 - tan(psi)**2)) / (1 + tan(psi)**2)
            cos_2pa = (cos_alpha * (1 - (tan_psi)**2) - 2 * sin_alpha * tan_psi) / (1 + (tan_psi)**2) # Combining equations 3 and 4b -- cos(2 * psi + alpha) = (cos(alpha) * (1 - tan(psi)**2) - 2 * sin(alpha) * tan(psi)) / (1 + tan(psi)**2)
        
        return C_amp, sin_2pa, cos_2pa

    def amplitude(self, f) -> np.array:
        """ GW amplitude -- as in equation 10 -- based on Apostolatos+ 1994
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency at which the amplitude is to be calculated [Hz]
        C_amp: amplitude of the beam pattern functions for + and x polarizations : Requires theta_S, phi_S
        sin_2pa: sine of 2 times the polarization angle + alpha -- for x polarization : Requires theta_S, phi_S, theta_J, phi_J
        cos_2pa: cosine of 2 times the polarization angle + alpha -- for + polarization : Requires theta_S, phi_S, theta_J, phi_J
        amp_prefactor: amplitude prefactor calculated using chirp mass and distance : Requires mcz, dist
        ___________________________________________________________________________________________
        Returns:
        amp: amplitude of the GW signal
        ___________________________________________________________________________________________
        LdotN(f) function returns Dot product of L and N: Dot product of the angular momentum vector and the line of sight vector at a given frequency
        beam_pattern_amplitude_and_phase(f) function returns the beam pattern functions amplitude and phase for + and x polarizations
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        LdotN = self.LdotN(f)
        C_amp, sin_2pa, cos_2pa = self.beam_pattern_amplitude_and_phase(f)

        amp = self.amp_prefactor() * C_amp * f**(-7/6) * np.sqrt(4 * LdotN**2 * sin_2pa**2 + cos_2pa**2 * (1+LdotN**2)**2)
        return amp

    def phase_phi_P(self, f):
        """ Polarization phase of the GW signal -- as in equation 11 -- based on Apostolatos+ 1994
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency at which the phase is to be calculated [Hz]
        sin_2pa: sine of 2 times the polarization angle + alpha -- for x polarization : Requires theta_S, phi_S, theta_J, phi_J
        cos_2pa: cosine of 2 times the polarization angle + alpha -- for + polarization : Requires theta_S, phi_S, theta_J, phi_J
        ___________________________________________________________________________________________
        Returns:
        phi_p: polarization phase of the GW signal
        ___________________________________________________________________________________________
        LdotN(f) function returns Dot product of L and N: Dot product of the angular momentum vector and the line of sight vector at a given frequency
        beam_pattern_amplitude_and_phase(f) function returns the beam pattern functions amplitude and phase for + and x polarizations
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        LdotN = self.LdotN(f)
        C_amp, sin_2pa, cos_2pa = self.beam_pattern_amplitude_and_phase(f)

        phi_p_temp = np.arctan2(2 * LdotN * sin_2pa, (1 + LdotN**2) * cos_2pa)
        phi_p = np.unwrap(phi_p_temp, discont=np.pi)
        return phi_p
    
    def f_dot(self, f):
        """ df/dt from Cutler & Flanagan 1994
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency at which the derivative is to be calculated [Hz]
        mcz: chirp mass [s]
        ___________________________________________________________________________________________
        Returns:
        f_dot: df/dt at a given frequency
        ___________________________________________________________________________________________
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        prefactor = (96/5) * np.pi**(8/3) * self.mcz**(5/3) * f**(11/3) # Leaving out the higher order terms
        return prefactor #* (1 - (743/336 + (11/4) * self.eta) * (np.pi * self.get_total_mass() * f)**(2/3) + 4 * np.pi * (np.pi * self.get_total_mass() * f)) #### Higher order terms

    def integrand_delta_phi(self, y, f):
        """Integrand for delta phi p (equations in Apostolatos 1994, and appendix A - Eq A18)
        For the correction to the precessional phase
        ___________________________________________________________________________________________
        Parameters used:
        y: variable for the integral : default is 0
        f: frequency at which the integrand is to be calculated [Hz]
        cos_i_JN: cosine of the angle between the total angular momentum and the line of sight : Requires theta_J, phi_J, theta_S, phi_S
        sin_i_JN: sine of the angle between the total angular momentum and the line of sight : Requires theta_J, phi_J, theta_S, phi_S
        omega_tilde: Dimensionless precession frequency -- stands for 1000 times the precession frequency at binary separation r = 6M for a solar mass binary
        ___________________________________________________________________________________________
        Returns:
        integrand_delta_phi: integrand for the delta phi p
        ___________________________________________________________________________________________
        LdotN(f) function returns Dot product of L and N: Dot product of the angular momentum vector and the line of sight vector at a given frequency
        get_theta_LJ(f) function returns the angle between L and J at a given frequency
        get_phi_LJ(f) function returns the angle between projection of L in the x-y plane and x axis in source frame at a given frequency
        f_dot(f) function returns the frequency derivative of the GW signal at a given frequency
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        LdotN = self.LdotN(f)
        cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = self.precession_angles()
        f_dot = self.f_dot(f)
        
        Omega_LJ = 1000 * self.omega_tilde * (f / self.get_f_cut())**(5/3) / (self.get_total_mass()/self.SOLMASS2SEC) # Equation 18b -- Omega_LJ = 1000 * omega_tilde * (f / f_cut)**(5/3) / (M / M_solar)  
        if abs(cos_i_JN) == 1:
            integrand_delta_phi = - Omega_LJ * np.cos(self.get_theta_LJ(f)) / f_dot # for face-on case
            
        else:
            integrand_delta_phi = (LdotN / (1 - LdotN**2)) * Omega_LJ * np.sin(self.get_theta_LJ(f)) * ( np.cos(self.get_theta_LJ(f)) * sin_i_JN * np.sin(self.get_phi_LJ(f)) - np.sin(self.get_theta_LJ(f)) * cos_i_JN ) / f_dot
        
        return integrand_delta_phi

    def phase_delta_phi(self, f):
        """ integrate the delta_phi integrand from 0 to f to get the phase correction
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency at which the phase is to be calculated [Hz]
        ___________________________________________________________________________________________
        Returns:
        integral: integral of the integrand_delta_phi
        ___________________________________________________________________________________________
        integrand_delta_phi(f) function returns the integrand for the delta phi p
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        integral = odeint(self.integrand_delta_phi, 0, f)
        return np.squeeze(integral)

    def Psi(self, f):
        """ GW phase upto 2 PN order
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency at which the phase is to be calculated [Hz]
        t_c: coalescence time
        phi_c: coalescence phase
        mcz: chirp mass [s]
        eta: symmetric mass ratio
        total mass of the binary system [s] -- get_total_mass() function : Requires chirp mass and symmetric mass ratio
        ___________________________________________________________________________________________
        Returns:
        Psi: GW phase of the GW signal
        ___________________________________________________________________________________________
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        x = (np.pi*self.get_total_mass()*f)**(2/3)
        Psi = (2 * np.pi * f * self.t_c) - self.phi_c - np.pi/4 + ((3/4)*(8*np.pi*self.mcz*f)**(-5/3)) * (1 + (20/9)*(743/336 + (11/4)*self.eta)*x - 16*np.pi*x**(3/2))
        return Psi

    
    def precessing_strain(self, f, delta_f=0.25, frequencySeries=True):
        """ GW signal with regular precession
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency at which the strain is to be calculated [Hz]
        delta_f: frequency resolution [Hz] : default is 0.25
        frequencySeries: boolean to return FrequencySeries object : default is True
        ___________________________________________________________________________________________
        Returns:
        precessing_strain: GW signal with regular precession
        ___________________________________________________________________________________________
        amplitude(f) function returns the amplitude of the GW signal at a given frequency
        Psi(f) function returns the GW phase of the GW signal at a given frequency
        phase_phi_P(f) function returns the polarization phase of the GW signal at a given frequency
        phase_delta_phi(f) function returns the phase correction to the precessional phase at a given frequency
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        precessing_strain = self.amplitude(f) * np.exp(1j*(self.Psi(f) - self.phase_phi_P(f) - 2 * self.phase_delta_phi(f)))
        if frequencySeries:
            return FrequencySeries(precessing_strain, delta_f, delta_f)
        return precessing_strain
    
    def cos_theta_L(self, f):
        """ For evolution of the orbital angular momentum vector in the detector frame
        From equation A8
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency at which the angle is to be calculated [Hz]
        cos_i_JN: cosine of the angle between the total angular momentum and the line of sight : Requires theta_J, phi_J, theta_S, phi_S
        sin_i_JN: sine of the angle between the total angular momentum and the line of sight : Requires theta_J, phi_J, theta_S, phi_S
        cos_o_XH: cosine of the angle between the x-axis of the source frame and H in the detector frame : Requires theta_J, phi_J, theta_S, phi_S
        sin_o_XH: sine of the angle between the x-axis of the source frame and H in the detector frame : Requires theta_J, phi_J, theta_S, phi_S
        ___________________________________________________________________________________________
        Returns:
        L_z: Cosine of the polar angle for the orbital angular momentum vector
        """
        cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = self.precession_angles()
        #from equation A8
        L_z = (np.sin(self.get_theta_LJ(f)) * (np.cos(self.get_phi_LJ(f)) * sin_o_XH + np.sin(self.get_phi_LJ(f)) * cos_i_JN * cos_o_XH) - sin_i_JN * cos_o_XH * np.cos(self.get_theta_LJ(f))) * np.sin(self.theta_S) + (np.sin(self.get_theta_LJ(f)) * np.sin(self.get_phi_LJ(f)) * sin_i_JN + np.cos(self.get_theta_LJ(f)) * cos_i_JN) * np.cos(self.theta_S)
        return L_z
    
    def phi_L(self, f):
        """ For evolution of the orbital angular momentum vector in the detector frame
        From equation A8
        ___________________________________________________________________________________________
        Parameters used:
        f: frequency at which the angle is to be calculated [Hz]
        cos_i_JN: cosine of the angle between the total angular momentum and the line of sight : Requires theta_J, phi_J, theta_S, phi_S
        sin_i_JN: sine of the angle between the total angular momentum and the line of sight : Requires theta_J, phi_J, theta_S, phi_S
        cos_o_XH: cosine of the angle between the x-axis of the source frame and H in the detector frame : Requires theta_J, phi_J, theta_S, phi_S
        sin_o_XH: sine of the angle between the x-axis of the source frame and H in the detector frame : Requires theta_J, phi_J, theta_S, phi_S
        ___________________________________________________________________________________________
        Returns:
        Phi_L: Phase of the orbital angular momentum vector
        ___________________________________________________________________________________________
        get_theta_LJ: Angle between L and J at a given frequency [rad] -- in the L-J plane
        get_phi_LJ: Phase of the L-J plane, this quantity is also frequency dependent
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        #from equation a8
        cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = self.precession_angles()
        L_H = np.sin(self.get_theta_LJ(f)) * (np.cos(self.get_phi_LJ(f)) * cos_o_XH - np.sin(self.get_phi_LJ(f)) * cos_i_JN * sin_o_XH) + sin_i_JN * sin_o_XH * np.cos(self.get_theta_LJ(f))
        L_V = np.sin(self.get_theta_LJ(f)) * (np.cos(self.get_phi_LJ(f)) * sin_o_XH + np.sin(self.get_phi_LJ(f)) * cos_i_JN * cos_o_XH) - sin_i_JN * cos_o_XH * np.cos(self.get_theta_LJ(f))
        L_N = np.sin(self.get_theta_LJ(f)) * np.sin(self.get_phi_LJ(f)) * sin_i_JN + np.cos(self.get_theta_LJ(f)) * cos_i_JN
        
        L_x = - np.sin(self.phi_S) * L_H - np.cos(self.theta_S) * np.cos(self.phi_S) * L_V + np.sin(self.theta_S) * np.cos(self.phi_S) * L_N
        L_y = np.cos(self.phi_S) * L_H - np.cos(self.theta_S) * np.sin(self.phi_S) * L_V + np.sin(self.theta_S) * np.sin(self.phi_S) * L_N
        Phi_L = np.arctan2(L_y, L_x)
        #Phi_L_ur = np.unwrap(Phi_L, discont = np.pi)
        return Phi_L#_ur
