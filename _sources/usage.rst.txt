Usage
=====

Basic Usage
-----------

The main scientific code is in Python modules in the repository root. Example usage:

.. code-block:: python

   import numpy as np
   from src.regular_precession import Regular_precession
   
   # Physical constants (defined in systems_lib.py)
   solar_mass = 4.92624076 * 1e-6      # [solar_mass] = sec
   giga_parsec = 1.02927125 * 1e17     # [giga_parsec] = sec
   
   # Define parameters for a regular precessing binary
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
   
   # Create regular precession instance
   rp = Regular_precession(params)
   
   # Generate waveform data
   f_min = 20.0  # Minimum frequency [Hz]
   f_cut = rp.get_f_cut()  # Cutoff frequency [Hz]
   delta_f = 0.25  # Frequency step [Hz]
   f_range = np.arange(f_min, f_cut, delta_f)
   
   # Get amplitudes and phases
   amplitudes = rp.amplitude(f_range)
   phases = rp.Psi(f_range) - (rp.phase_phi_P(f_range) + 2*rp.phase_delta_phi(f_range))

Jupyter Notebooks
-----------------

See the `notebooks/` directory for example analyses and figures.

For detailed API documentation of all modules, see the :doc:`api` page.
