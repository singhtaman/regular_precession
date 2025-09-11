Usage
=====

Basic Usage
-----------

The main scientific code is in Python modules in the repository root. Example usage:

.. code-block:: python

   from regular_precession import Regular_precession
   params = {
       'theta_S': 0.5, 
       'phi_S': 0.0, 
       'theta_J': 1.0, 
       'phi_J': 0.5,
       'mcz': 10, 
       'dist': 1.5, 
       'eta': 0.25, 
       't_c': 0, 
       'phi_c': 0,
       'theta_tilde': 4.0, 
       'omega_tilde': 2.0, 
       'gamma_P': 0.0
   }
   rp = Regular_precession(params)
   amp = rp.amplitude(30)

Jupyter Notebooks
-----------------

See the `notebooks/` directory for example analyses and figures.

API Reference
-------------

.. automodule:: regular_precession
   :members:
   :undoc-members:
   :show-inheritance:
