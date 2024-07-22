# regular_precession

### Regular precession waveforms and their distinguishability from non precessing signals

This repository contains codes supporting the unpublished (soon to be submitted) work on regular precession waveforms.

## Contents

- [Description](#description)
- [Pre-requisites](#prerequisites)
- [License](#license)
- [Acknowlegement](#acknowledgement)

## Description

Stellar binary black holes coalescing in LIGO band can have misaligned spins depending on formation channels of the binary. With misaligned spins, the orbital angular momentum (L) precesses and nutates about the total angular momentum (J) leaving behind modulations in the GW signal. Regular precession is defined when L moves around in a cone around J (no nutation). In this special case of a vast family of generically precessing sources, the direction of J is nearly constant on the radiation reaction timescale.

We introduce 2 dimensionless parameters in the PN approximation frequency domain waveforms. These define the opening angle of the cone and how fast L moves around J in the regualr precession case. Then we investigate how these parameters (along with a nuisance parameter: precessional phase) effect the GW amplitude and phase. We also calculate mismatches between a precessing signal and a non-precessing template bank to establish criterion for distinguishablity of precession and probability of observing a precessing signal with current aLIGO sensitivities.

In this repository, you will find `python` scripts and jupyter notebooks that support the data and findings we report in our (unpublished) work.


## Pre-requisites
This project requires the installation of [`PyCBC`](https://pycbc.org). For a guide on how to install, I find Davide Gerosa's [notes](https://davidegerosa.com/installlal/) realy helpful. 

## Acknowledgement

This work was supported by the National Science Foundation Grant No. PHY-2309320. We also acknowledge the Texas Advanced Computing Center ([TACC](http://www.tacc.utexas.edu)) at The University of Texas at Austin for providing HPC resources that have contributed to the research results reported within this paper [Stanzione et al](https://doi.org/10.1145/3093338.3093385).


## License

This project is licensed under the [MIT License](LICENSE).