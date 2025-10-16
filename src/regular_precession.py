"""
Regular Precession: When L moves in a cone around J with an opening angle theta_tilde that changes on a radiation reaction timescale,
with frequency omega_tilde (also changing on the same timescale) and a phase gamma_P.

Model presented in following paper : arXiv:2509.10628 [gr-qc]
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
        """
        Calculate the total mass of the binary system from the chirp mass and symmetric mass ratio [s], M = M_c / eta**(3/5).

        Parameters
        ----------
        mcz : float
            Chirp mass [s].
        eta : float
            Symmetric mass ratio [dimensionless].

        Returns
        -------
        total_mass : float
            Total mass of the binary system.
        """
        return self.mcz/(self.eta**(3/5))           # total mass of the binary system [s]: M = M_c / eta**(3/5)

    def get_f_cut(self):
        """
        Compute the cut-off frequency where the binary coalesces, equation 13 in the paper (https://arxiv.org/pdf/2509.10628), f_cut = 1/(6**(3/2)*pi*M), where M is the total mass of the binary.
        
        Parameters
        ----------
        get_total_mass : function call
            Computes the total mass of the binary system given the chirp mass and symmetric mass ratio, M = M_c / eta**(3/5).

        Returns
        -------
        f_cut : float
            Cut-off frequency.
        """
        return 1/(6**(3/2)*np.pi*self.get_total_mass()) # Equation 13 -- cut-off frequency [Hz]: f_cut = 1/(r_{ISCO}**(3/2) * pi * M)
    

    def get_theta_LJ(self, f):
        """
        Angle between L and J at a given frequency [rad] in the L-J plane, equation 18a in the paper (https://arxiv.org/pdf/2509.10628), theta_LJ = 0.1/(4*eta) * theta_tilde * (f/f_cut)**(1/3), where f_cut is the cut-off frequency, theta_tilde is the dimensionless precession amplitude.

        Parameters
        ----------
        f : float or array_like
            Frequency where the angle is to be calculated [Hz].
        eta : float
            Symmetric mass ratio.
        theta_tilde : float
            Dimensionless precession amplitude parameter.
        get_f_cut : function call
            Computes the cut-off frequency, f_cut = 1/(6**(3/2)*pi*M), where M is the total mass of the binary.

        Returns
        -------
        theta_LJ : float or array_like
            Angle between L and J at a given frequency.
        """
        return (0.1/(4*self.eta))*self.theta_tilde*(f/self.get_f_cut())**(1/3) # Equation 18a  -- Angle between L and J at a given frequency [rad]: theta_LJ = 0.1/(4*eta) * theta_tilde * (f/f_cut)**(1/3)
    
    def get_phi_LJ(self, f):
        """
        Angle between projection of L in the x-y plane and x axis in source frame at a given frequency [rad], equation 19 in the paper (https://arxiv.org/pdf/2509.10628), phi_LJ = gamma_P + \int_{f_min}^f (Omega_LJ (df'/dt)**(-1)) df', where Omega_LJ is the precession frequency.

        Parameters
        ----------
        f : float or array_like
            Frequency where the angle is to be calculated [Hz].
        omega_tilde : float
            Dimensionless precession frequency parameter [Dimensionless].
        mcz : float
            Chirp mass [s].
        gamma_P : float
            Precession phase at reference frequency [rad].
        get_total_mass : function call
            Computes the total mass of the binary system, M = M_c / eta**(3/5), given the chirp mass M_c and symmetric mass ratio eta [s].
        get_f_cut : function call
            Computes the cut-off frequency, f_cut = 1/(6**(3/2)*pi*M), where M is the total mass of the binary [Hz].

        Returns
        -------
        phi_LJ : float or array_like
            Angle between projection of L in the x-y plane and x axis in source frame at a given frequency.
        """
        phi_LJ_amp = (5000 * self.omega_tilde) / (96 * (self.get_total_mass()/self.SOLMASS2SEC) * (np.pi**(8/3)) * (self.mcz**(5/3)) * (self.get_f_cut()**(5/3)))
        return phi_LJ_amp * (1/self.FMIN - 1/f) + self.gamma_P # Equation 19 (also uses equation 18b) -- Angle between projection of L in the x-y plane and x axis in source frame at a given frequency [rad]: phi_LJ = phi_LJ_amp * (1/f_min - 1/f) + gamma_P
        
        
    def amp_prefactor(self) -> float:
        """
        Amplitude prefactor calculated using chirp mass and distance, equation 6 in the paper (https://arxiv.org/pdf/2509.10628), A = sqrt(5/96) * (pi**(-2/3)) * (M_c**(5/6)) / D_L, where M_c is the chirp mass and D_L is the luminosity distance.

        Parameters
        ----------
        mcz : float
            Chirp mass [s].
        dist : float
            Luminosity distance [s].

        Returns
        -------
        amp_prefactor : float
            Amplitude prefactor.
        """
        amp_prefactor = np.sqrt(5/96) * (np.pi**(-2/3)) * (self.mcz**(5/6)) / self.dist # from equation 6 -- A = sqrt(5/96) * (pi**(-2/3)) * (M_c**(5/6)) / D_L
        return amp_prefactor

    def precession_angles(self):
        """
        Compute angles important for the precession model, specifically the angle between J and N, and the angle between the x-axis of the source frame and H in the sky frame. Equations A4, A6a, A6b in the paper (https://arxiv.org/abs/2509.10628).

        Parameters
        ----------
        theta_J : float
            Binary plane inclination - polar angle for J in detector frame [rad].
        phi_J : float
            Binary plane azimuthal angle for J in detector frame [rad].
        theta_S : float
            Sky inclination - polar angle for line of sight in detector frame [rad].
        phi_S : float
            Sky azimuthal angle for line of sight in detector frame [rad].

        Returns
        -------
        cos_i_JN : float
            Cosine of the angle between the total angular momentum and the line of sight.
        sin_i_JN : float
            Sine of the angle between the total angular momentum and the line of sight.
        cos_o_XH : float
            Cosine of the angle between the x-axis of the source frame and H in the detector frame.
        sin_o_XH : float
            Sine of the angle between the x-axis of the source frame and H in the detector frame.
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
        """
        Cosine of the angle between L and N, equation A10 in the paper (https://arxiv.org/pdf/2509.10628).

        Parameters
        ----------
        f : float or array_like
            Frequency at which to calculate the dot product [Hz].
        get_theta_LJ : function call
            Computes the angle between L and J, equation 18a.
        get_phi_LJ : function call
            Computes the azimuthal angle of L in source frame, equation 18b.
        precession_angles : function call
            Computes orientation angles, equations A4, A6a, A6b.

        Returns
        -------
        LdotN : float or array_like
            Dot product of L and N (cosine of angle between L and N).
        """
        cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = self.precession_angles()
        LdotN = np.sin(self.get_theta_LJ(f)) * sin_i_JN * np.sin(self.get_phi_LJ(f)) + np.cos(self.get_theta_LJ(f)) * cos_i_JN  
        return LdotN

    def beam_pattern_amplitude_and_phase(self, f):
        """
        Beam pattern functions for + and x polarizations, equations 3, 4a, 4b in the paper (https://arxiv.org/pdf/2509.10628).

        Parameters
        ----------
        f : float or array_like
            Frequency at which to calculate the beam pattern [Hz].
        theta_S : float
            Sky inclination angle [rad].
        phi_S : float
            Sky azimuthal angle [rad].
        get_theta_LJ : function call
            Computes the angle between L and J, equation 18a.
        get_phi_LJ : function call
            Computes the azimuthal angle of L in source frame, equation 18b.
        precession_angles : function call
            Computes orientation angles, equations A4, A6a, A6b.

        Returns
        -------
        C_amp : float or array_like
            Amplitude of the beam pattern functions for + and x polarizations.
        sin_2pa : float or array_like
            Sine of 2 times the polarization angle + alpha (for x polarization).
        cos_2pa : float or array_like
            Cosine of 2 times the polarization angle + alpha (for + polarization).
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
        """
        GW amplitude, equation 10 in the paper (https://arxiv.org/pdf/2509.10628) - also in Apostolatos+ 1994 as equation 7a.
        
        Parameters
        ----------
        f : float or array_like
            Frequency at which the amplitude is to be calculated [Hz].
        LdotN : function call
            Computes the dot product between L and N, equation A10.
        beam_pattern_amplitude_and_phase : function call
            Computes beam pattern functions, equations 3, 4a, 4b.
        amp_prefactor : function call
            Computes amplitude prefactor, equation 6.

        Returns
        -------
        amp : float or array_like
            Amplitude of the GW signal.
        """
        LdotN = self.LdotN(f)
        C_amp, sin_2pa, cos_2pa = self.beam_pattern_amplitude_and_phase(f)

        amp = self.amp_prefactor() * C_amp * f**(-7/6) * np.sqrt(4 * LdotN**2 * sin_2pa**2 + cos_2pa**2 * (1+LdotN**2)**2)
        return amp

    def phase_phi_P(self, f):
        """
        Polarization phase of the GW signal, equation 11 in the paper (https://arxiv.org/pdf/2509.10628).
        
        Parameters
        ----------
        f : float or array_like
            Frequency at which the phase is to be calculated [Hz].
        LdotN : function call
            Computes the dot product between L and N, equation A10.
        beam_pattern_amplitude_and_phase : function call
            Computes beam pattern functions, equations 3, 4a, 4b.

        Returns
        -------
        phi_p : float or array_like
            Polarization phase of the GW signal.
        """
        LdotN = self.LdotN(f)
        C_amp, sin_2pa, cos_2pa = self.beam_pattern_amplitude_and_phase(f)

        phi_p_temp = np.arctan2(2 * LdotN * sin_2pa, (1 + LdotN**2) * cos_2pa)
        phi_p = np.unwrap(phi_p_temp, discont=np.pi)
        return phi_p
    
    def f_dot(self, f):
        """
        df/dt from Cutler & Flanagan 1994, equation 20 in the paper (https://arxiv.org/pdf/2509.10628).

        Parameters
        ----------
        f : float or array_like
            Frequency at which the derivative is to be calculated [Hz].
        mcz : float
            Chirp mass [s].

        Returns
        -------
        f_dot : float or array_like
            df/dt at a given frequency [Hz/s].
        """
        prefactor = (96/5) * np.pi**(8/3) * self.mcz**(5/3) * f**(11/3) # Leaving out the higher order terms
        return prefactor #* (1 - (743/336 + (11/4) * self.eta) * (np.pi * self.get_total_mass() * f)**(2/3) + 4 * np.pi * (np.pi * self.get_total_mass() * f)) #### Higher order terms

    def integrand_delta_phi(self, y, f):
        """
        Integrand for delta phi p, equation 12 and A19 in the paper (https://arxiv.org/pdf/2509.10628).

        Parameters
        ----------
        y : float
            Variable for the integral (typically 0 for initial condition).
        f : float or array_like
            Frequency at which the integrand is to be calculated [Hz].
        omega_tilde : float
            Dimensionless precession frequency parameter.
        LdotN : function call
            Computes the dot product between L and N, equation A10.
        precession_angles : function call
            Computes orientation angles, equations A11, A12, A13.
        f_dot : function call
            Computes frequency derivative, equation 20.
        get_theta_LJ : function call
            Computes the angle between L and J, equation A14.
        get_phi_LJ : function call
            Computes the azimuthal angle of L in source frame, equation A15.
        get_f_cut : function call
            Computes the cut-off frequency, equation A16.
        get_total_mass : function call
            Computes the total mass of the binary, equation A17.

        Returns
        -------
        integrand_delta_phi : float or array_like
            Integrand for the delta phi p.
        """
        LdotN = self.LdotN(f)
        cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = self.precession_angles()
        f_dot = self.f_dot(f)
        
        Omega_LJ = 1000 * self.omega_tilde * (f / self.get_f_cut())**(5/3) / (self.get_total_mass()/self.SOLMASS2SEC) # Equation 18b -- Omega_LJ = 1000 * omega_tilde * (f / f_cut)**(5/3) / (M / M_solar)  
        if abs(cos_i_JN) == 1:
            integrand_delta_phi = - Omega_LJ * np.cos(self.get_theta_LJ(f)) / f_dot # for face-on case
            
        else:
            integrand_delta_phi = (LdotN / (1 - LdotN**2)) * ((Omega_LJ * np.sin(self.get_theta_LJ(f)) * ( np.cos(self.get_theta_LJ(f)) * sin_i_JN * np.sin(self.get_phi_LJ(f)) - np.sin(self.get_theta_LJ(f)) * cos_i_JN ) / f_dot) - ( (self.get_theta_LJ(f)/(3*f)) * np.cos(self.get_phi_LJ(f)) * sin_i_JN)) # added correction term with theta_LJ/3f cos phi_LJ sin i_JN term
        
        return integrand_delta_phi

    def phase_delta_phi(self, f):
        """
        Integrates the delta_phi integrand from 0 to f to get the phase correction to precession phase, integrating equation 12/A19 in the paper (https://arxiv.org/pdf/2509.10628).

        Parameters
        ----------
        f : float or array_like
            Frequency at which the phase is to be calculated [Hz].
        integrand_delta_phi : function call
            Computes the integrand for phase correction, equation 12/A19.

        Returns
        -------
        integral : float or array_like
            Integral of the integrand_delta_phi.
        """
        integral = odeint(self.integrand_delta_phi, 0, f)
        return np.squeeze(integral)

    def Psi(self, f):
        """
        GW phase up to 2 PN order, equation 7 in https://arxiv.org/pdf/2509.10628.
        
        Parameters
        ----------
        f : float or array_like
            Frequency at which the phase is to be calculated [Hz].
        mcz : float
            Chirp mass [s].
        eta : float
            Symmetric mass ratio.
        t_c : float
            Coalescence time [s].
        phi_c : float
            Coalescence phase [rad].
        get_total_mass : function call
            Computes the total mass of the binary.

        Returns
        -------
        Psi : float or array_like
            GW phase of the GW signal.
        """
        x = (np.pi*self.get_total_mass()*f)**(2/3)
        Psi = (2 * np.pi * f * self.t_c) - self.phi_c - np.pi/4 + ((3/4)*(8*np.pi*self.mcz*f)**(-5/3)) * (1 + (20/9)*(743/336 + (11/4)*self.eta)*x - 16*np.pi*x**(3/2))
        return Psi

    
    def precessing_strain(self, f, delta_f=0.25, frequencySeries=True):
        """
        GW signal with regular precession, equation 9 in https://arxiv.org/pdf/2509.10628.

        Parameters
        ----------
        f : float or array_like
            Frequency at which the strain is to be calculated [Hz].
        delta_f : float, optional
            Frequency resolution [Hz]. Default is 0.25.
        frequencySeries : bool, optional
            Whether to return a FrequencySeries object. Default is True.
        amplitude : function call
            Computes the GW amplitude, equation 10.
        Psi : function call
            Computes the GW phase, equation 7.
        phase_phi_P : function call
            Computes the polarization phase, equation 11.
        phase_delta_phi : function call
            Computes the precession phase correction, equation A19.

        Returns
        -------
        precessing_strain : array or FrequencySeries
            GW signal with regular precession.
        """
        precessing_strain = self.amplitude(f) * np.exp(1j*(self.Psi(f) - self.phase_phi_P(f) - 2 * self.phase_delta_phi(f)))
        if frequencySeries:
            return FrequencySeries(precessing_strain, delta_f, delta_f)
        return precessing_strain
    
    def cos_theta_L(self, f):
        """
        Evolution of the orbital angular momentum vector in the detector frame (cosine of polar angle).

        Parameters
        ----------
        f : float or array_like
            Frequency at which the angle is to be calculated [Hz].
        theta_S : float
            Sky inclination angle [rad].
        get_theta_LJ : function call
            Computes the angle between L and J.
        get_phi_LJ : function call
            Computes the azimuthal angle of L in source frame.
        precession_angles : function call
            Computes orientation angles, equations A4, A6a, A6b.

        Returns
        -------
        L_z : float or array_like
            Cosine of the polar angle for the orbital angular momentum vector.
        """
        cos_i_JN, sin_i_JN, cos_o_XH, sin_o_XH = self.precession_angles()
        #from equation A8
        L_z = (np.sin(self.get_theta_LJ(f)) * (np.cos(self.get_phi_LJ(f)) * sin_o_XH + np.sin(self.get_phi_LJ(f)) * cos_i_JN * cos_o_XH) - sin_i_JN * cos_o_XH * np.cos(self.get_theta_LJ(f))) * np.sin(self.theta_S) + (np.sin(self.get_theta_LJ(f)) * np.sin(self.get_phi_LJ(f)) * sin_i_JN + np.cos(self.get_theta_LJ(f)) * cos_i_JN) * np.cos(self.theta_S)
        return L_z
    
    def phi_L(self, f):
        """
        Evolution of the orbital angular momentum vector in the detector frame (azimuthal angle).

        Parameters
        ----------
        f : float or array_like
            Frequency at which the angle is to be calculated [Hz].
        theta_S : float
            Sky inclination angle [rad].
        phi_S : float
            Sky azimuthal angle [rad].
        get_theta_LJ : function call
            Computes the angle between L and J, equation 18a.
        get_phi_LJ : function call
            Computes the azimuthal angle of L in source frame, equation 19.
        precession_angles : function call
            Computes orientation angles, equations A4, A6a, A6b.

        Returns
        -------
        Phi_L : float or array_like
            Azimuthal angle of the orbital angular momentum vector.
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
