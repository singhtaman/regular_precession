# regular_precession

[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-2509.10628-b31b1b.svg)](https://arxiv.org/abs/2509.10628)

### Regular precession waveforms and their distinguishability from non-precessing signals

This repository contains codes supporting the work on regular precession waveforms (arXiv: 2509.10628).

## Contents

- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Repository Structure](#repository-structure)
- [Pre-requisites](#pre-requisites)
- [Documentation](#documentation)
- [Citation](#citation)
- [Acknowledgement](#acknowledgement)
- [License](#license)


## Description

Stellar binary black holes coalescing in LIGO band can have misaligned spins depending on formation channels of the binary. With misaligned spins, the orbital angular momentum ($\vec{L}$) precesses and nutates about the total angular momentum ($\vec{J}$) leaving behind modulations in the GW signal. Regular precession is defined when $\vec{L}$ moves around in a cone around $\vec{J}$ (with no nutation). In this special case of a vast family of generically precessing sources, the direction $\hat{J}$ is nearly constant on the radiation reaction timescale.

We introduce 2 dimensionless parameters in the PN approximation to the GW waveforms. These define the opening angle of the cone that $\vec{L}$ makes with $\vec{J}$ and how fast $\vec{L}$ moves around $\vec{J}$ in the regular precession case. Then we investigate how these parameters (along with a nuisance parameter: precessional phase) effect the GW amplitude and phase. We also calculate mismatches between a precessing signal and a non-precessing template bank to establish criterion for distinguishablity of precession and probability of observing a precessing signal with current aLIGO sensitivities. We explore different spin alignment scenarios and spin magnitudes to find corresponding distributions of precession parameters and Monte Carlo over sky locations & binary orientations at different cosmological redshifts to discuss the fractions of distinguishable precessing sources for current sensitivities.

In this repository, you will find python scripts and jupyter notebooks that support the data and findings we report in our work.


## Installation

### Option 1: Install from source (recommended)
```bash
git clone https://github.com/singhtaman/regular_precession.git
cd regular_precession
pip install -e .
```

### Option 2: Direct pip install
```bash
pip install git+https://github.com/singhtaman/regular_precession.git
```

## Usage

### Basic Example

```python
from src.regular_precession import Regular_precession

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

# Create regular precession waveform
rp = Regular_precession(params)
# Use the class methods to generate waveforms and analyze precession effects
```

### Quick Test

After installation, run the provided example script to test your setup:

```bash
python example.py
```

## Examples

The `notebooks/` directory contains comprehensive Jupyter notebooks demonstrating various aspects of the regular precession model:

- **`Amplitude.ipynb`**: Explores amplitude modulations in regular precessing signals
- **`Phase-phi_p_and_delta_phi.ipynb`**: Analyzes phase evolution and precession effects
- **`Modulation_L.ipynb`**: Studies orbital angular momentum modulations
- **`Edge_on_systems.ipynb`**: Examines edge-on binary systems
- **`Validations_PN_approx.ipynb`**: Validates the post-Newtonian approximations used
- **Mismatch analysis notebooks**: Calculate distinguishability metrics and Monte Carlo studies

## Repository Structure

```
regular_precession/
├── LICENSE                       # MIT license file
├── README.md                     # This file
├── pyproject.toml               # Project configuration and dependencies
├── requirements.txt             # Alternative dependency specification
├── example.py                   # Simple example script for testing
├── src/                         # Source code directory
│   ├── regular_precession.py    # Main regular precession class
│   ├── mismatch_n_SNR.py       # Mismatch and SNR calculations
│   ├── systems_lib.py          # Binary system utilities
│   └── MC_skyloc_binorient_plus_prec.py  # Monte Carlo simulations
├── notebooks/                   # Jupyter notebooks with examples and analysis
│   ├── Amplitude.ipynb          # Amplitude modulation analysis
│   ├── Phase-phi_p_and_delta_phi.ipynb  # Phase evolution studies
│   ├── Modulation_L.ipynb       # Orbital angular momentum modulations
│   ├── Edge_on_systems.ipynb    # Edge-on binary analysis
│   ├── Validations_PN_approx.ipynb  # PN approximation validation
│   ├── optimized_mismatch_*.ipynb    # Mismatch contour studies
│   ├── table_*.ipynb           # Monte Carlo fraction studies
│   ├── figs/                   # Generated figures and plots
│   └── saved_data/             # Cached computation results and data
│       ├── *.pkl              # Pickled data files
│       └── TACC_runs/         # High-performance computing results
└── docs/                       # Documentation (Sphinx-generated)
    ├── *.html                  # HTML documentation files
    ├── *.rst                   # RestructuredText source files
    ├── _static/                # Static assets (CSS, JS, images)
    └── _modules/               # API documentation modules
```


## Pre-requisites

This project requires Python ≥ 3.9 and the following key dependencies:

- [`lalsuite`](https://pypi.org/project/lalsuite/) - LIGO Algorithm Library for gravitational wave data analysis
- [`PyCBC`](https://pycbc.org) - Python toolkit for gravitational wave astronomy
- Standard scientific Python packages: `numpy`, `scipy`, `matplotlib`, `jupyter`

### Installing LAL and PyCBC

For installation guidance, we recommend Davide Gerosa's excellent [installation notes](https://davidegerosa.com/installlal/).

#### Quick installation via conda (recommended):
```bash
conda install -c conda-forge lalsuite pycbc
```

#### Alternative pip installation:
```bash
pip install lalsuite pycbc
```

## Documentation

Full API documentation is available in the `docs/` directory. You can build the documentation locally or view the generated HTML files:

- **Online**: [singhtaman.github.io/regular_precession](https://singhtaman.github.io/regular_precession)
- **Local HTML**: Open `docs/index.html` in your browser for the main documentation
- **API Reference**: `docs/api.html` contains detailed API references
- **Usage Guide**: `docs/usage.html` provides usage examples

### Building Documentation Locally

To rebuild the documentation after making changes:

```bash
./rebuild_docs.sh
```

Or manually with Sphinx:

```bash
cd docs
sphinx-build -b html . _build/html
```

## Citation

If you use this code in your research, please cite:

```bibtex
@article{Singh2024RegularPrecession,
       author = {{Singh}, Tamanjyot and {Stoikos}, Evangelos and {Ali}, Saif and {Steinle}, Nathan and {Kesden}, Michael and {King}, Lindsay},
        title = "{Detecting regular precession using a new gravitational waveform model directly parameterized by both precession amplitude and frequency}",
      journal = {arXiv e-prints},
     keywords = {General Relativity and Quantum Cosmology, High Energy Astrophysical Phenomena},
         year = 2025,
        month = sep,
          eid = {arXiv:2509.10628},
        pages = {arXiv:2509.10628},
          doi = {10.48550/arXiv.2509.10628},
archivePrefix = {arXiv},
       eprint = {2509.10628},
 primaryClass = {gr-qc},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025arXiv250910628S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

``` 

## Acknowledgement

This work was supported by the National Science Foundation Grant No. PHY-2309320. 

We also acknowledge the Texas Advanced Computing Center ([TACC](http://www.tacc.utexas.edu)) at The University of Texas at Austin for providing HPC resources that have contributed to the research results reported within this repository ([Stanzione et al., 2020](https://doi.org/10.1145/3093338.3093385)).


## License

This project is licensed under the [MIT License](LICENSE).