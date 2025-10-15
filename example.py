#!/usr/bin/env python3
"""
Simple example demonstrating basic usage of the regular_precession package.

This script creates a regular precessing binary system and demonstrates
basic waveform generation capabilities.
"""

import numpy as np
import sys
import os

# Add the src directory to the path for importing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

solar_mass = 4.92624076 * 1e-6          # [solar_mass] = sec
giga_parsec = 1.02927125 * 1e17         # [giga_parsec] = sec 
# You can also find this in systems_lib.py - impoert it from there if needed

try:
    from regular_precession import Regular_precession
    print("✓ Successfully imported regular_precession module")
except ImportError as e:
    print(f"✗ Failed to import regular_precession: {e}")
    print("Please ensure all dependencies are installed.")
    sys.exit(1)

def main():
    """Run a basic example of regular precession waveform generation."""
    
    print("\n" + "="*50)
    print("Regular Precession Example")
    print("="*50)
    
    # Define parameters for a typical regular precessing binary
    params = {
        "theta_S": np.pi / 4,           # Sky location (N) - polar angle
        "phi_S": 0.0,                   # Sky location (N) - azimuthal angle
        "theta_J": np.pi / 4,           # Binary orientation (J) - polar angle
        "phi_J": 0.0,                   # Binary orientation (J) - azimuthal angle
        "mcz": 10 * solar_mass,         # Chirp mass in seconds
        "dist": 1.5 * giga_parsec,      # Distance in seconds
        "eta": 0.25,                    # Symmetric mass ratio
        "t_c": 0.0,                     # Time of coalescence
        "phi_c": 0.0,                   # Phase of coalescence
        "theta_tilde": 4.0,             # Precession amplitude - dimensionless
        "omega_tilde": 2.0,             # Precession frequency - dimensionless
        "gamma_P": 0.0,                 # Precessional phase at reference
    }
    
    print(f"Creating regular precessing binary with:")
    print(f"  Chirp mass: {params['mcz']:.1f} M☉")
    print(f"  Distance: {params['dist']:.1f} Gpc")
    print(f"  Symmetric mass ratio: {params['eta']:.2f}")
    print(f"  Precession amplitude: {params['theta_tilde']:.2f}")
    print(f"  Precession frequency: {params['omega_tilde']:.2f}")
    
    try:
        # Create regular precession instance
        rp = Regular_precession(params)
        print("\n✓ Successfully created Regular_precession instance")
        
        # Need to define a frequency range for waveform generation
        f_min = 20.0  # Minimum frequency [Hz]
        f_cut = rp.get_f_cut()  # Cutoff frequency [Hz]
        delta_f = 0.25  # Frequency step [Hz]
        f_range = np.arange(f_min, f_cut, delta_f)
        
        # Get total mass
        total_mass = rp.get_total_mass()
        print(f"✓ Total mass: {total_mass:.1f} M☉")
        print(f"\nAmplitude and phase as frequency arrays:")
        print(f"  Amplitudes: {rp.amplitude(f_range)}")
        print(f"  Phases: {rp.Psi(f_range) - (rp.phase_phi_P(f_range) + 2*rp.phase_delta_phi(f_range))}")
        
    except Exception as e:
        print(f"✗ Error during waveform generation: {e}")
        return 1
    
    print(f"\n{'='*50}")
    print("Example completed successfully!")
    print("See notebooks/ directory for more detailed examples.")
    print("="*50)
    
    return 0

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)