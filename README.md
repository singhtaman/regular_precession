# regular_precession

### Regular precession waveforms and their distinguishability from non-precessing signals

This repository contains codes supporting the unpublished (soon to be submitted) work on regular precession waveforms.

## Contents

- [Description](#description)
- [Pre-requisites](#pre-requisites)
- [Acknowlegement](#acknowledgement)
- [License](#license)


## Description

Stellar binary black holes coalescing in LIGO band can have misaligned spins depending on formation channels of the binary. With misaligned spins, the orbital angular momentum ($\vec{L}$) precesses and nutates about the total angular momentum ($\vec{J}$) leaving behind modulations in the GW signal. Regular precession is defined when $\vec{L}$ moves around in a cone around $\vec{J}$ (with no nutation). In this special case of a vast family of generically precessing sources, the direction $\hat{J}$ is nearly constant on the radiation reaction timescale.

We introduce 2 dimensionless parameters in the PN approximation to the GW waveforms. These define the opening angle of the cone that $\vec{L}$ makes with $\vec{J}$ and how fast $\vec{L}$ moves around $\vec{J}$ in the regular precession case. Then we investigate how these parameters (along with a nuisance parameter: precessional phase) effect the GW amplitude and phase. We also calculate mismatches between a precessing signal and a non-precessing template bank to establish criterion for distinguishablity of precession and probability of observing a precessing signal with current aLIGO sensitivities. We explore different spin alignment scenarios and spin magnitudes to find corresponding distributions of precession parameters and Monte Carlo over sky locations & binary orientations at different cosmological redshifts to discuss the fractions of distinguishable precessing sources for current sensitivities.

In this repository, you will find python scripts and jupyter notebooks that support the data and findings we report in our (soon to be submitted) work.


## Pre-requisites
This project requires the installation of [`lal`](https://pypi.org/project/lalsuite/) and [`PyCBC`](https://pycbc.org). For a guide on how to install, I find Davide Gerosa's [notes](https://davidegerosa.com/installlal/) realy helpful. 

## Acknowledgement

This work was supported by the National Science Foundation Grant No. PHY-2309320. 

We also acknowledge the Texas Advanced Computing Center ([TACC](http://www.tacc.utexas.edu)) at The University of Texas at Austin for providing HPC resources that have contributed to the research results reported within this repository [Stanzione et al](https://doi.org/10.1145/3093338.3093385).


## License

This project is licensed under the [MIT License](LICENSE).