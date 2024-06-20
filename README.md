# Computational Spectroscopy
Copyright 2024, Zeke A. Piskulich and Qiang Cui (Boston University).

Software is adapted from programs originally developed Copyright 2023 Dr. Ashley Borkowski, Dr. Hasini Senanayake, and Professor Ward Thompson (University of Kansas).

## Introduction

This code is a computational spectroscopic toolkit for calculating spectroscopic observables from molecular dynamics simulations.

Thus far, the currently implemented spectroscopies are as follows:

* One Dimensional Infrared Spectroscopy
   * Dipole Spectral Density
   * Frequency Distribution
* Two Dimentional Infrared Spectroscopy
* One Dimensional Raman Spectroscopy
* One Dimensional Sum-Frequency Generation Spectroscopy
* Frequency Frequency Time Correlation Function

This code works through the vibrational frequency mapping approach developed by Corcelli, Skinner, and others, and uses this approach to calculate spectroscopic observables from classical molecular dynamics simulations.

For a detailed description of the theory underlying this approach, please see the following publication:

[ Insert Publication Here ]

## Limitations

For currently identified limitations, please see the issue tracker.

## Requirements

1) GCC Fortran Compiler. I have tested the following codes with GGC 10.8
2) FFTW3 Library
3) HDF5 Library

Additionally, for installation, it uses a CMAKE build process, so you will need a version of cmake installed during this process.

## Installation

This package can be installed with cmake. The general approach for such a build is the following:

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../bin/
make
make install
```

## Usage 

### Input File Format

```
#noh
1710 
# ntimes, timestep, ncorr, nskip
2500 4.0d0 250 250
#wfreq
4.0
# flag_hist, nhist
.false. 100
#w1min w1max w2min w2max
3000.0 3800.0 3000.0 3800.0
#zmin, zmax, zhist
-35 35 50
#zcenter
  34.78215
# NTw
5
#Tw(1), Tw(2)...Tw(NTw)
0.0 200.0 400.0 600.0 800.0
```

Here, the meaning of the input arguments are described:

* noh: number of ohs (typically 2x number of waters)
* ntimes: number of times to calculate
* timestep: timestep of configurations in calculation
* ncorr: correlation length (in steps)
* nskip: steps to skip between origins
* wfreq: resolution of frequencies
* flag_hist: Should it do histogramming? (deprecated)
* nhist: number of histogram bins
* w1_min, w1max : frequency range for spectrum
* w3_min, w3_max : frequency range for 2d spectrum (second axis)
* zmin, zmax: min and max z distance for 2d spectral density
* zhist: number of z dimension histogram bins
* zcenter: center of box (for sfg)
* Ntw: number of waiting times
* Tw(1)...: waiting times in fs (list)

### Command Line Arguments

Command line arguments for the present software act effectively to override the inputs recieved from the input files for specific use cases (for instance, if you mostly want the spectra you calculate to use set parameters, but for the FFCF want longer TCFs).

The CLI argument options that are available are the following:

    -h, help
        -> Displays a useful help message.
    -ncorr [value]
        -> overrides ncorr value from input file
    -nskip [value]
        -> overrides nskip value from input file
    -in [filename]
        -> Reads filename as input file 
           rather than spectra.in
    -map [filename]
        -> Reads filename as empirical map file
           rather than empirical_map.in
    -tag [tagname]
        -> This CLI labels output files with tagname
           so sfg_tcf.dat would become 
           new_sfg_tcf.dat with the command -tag new_
    -avfreq [freq]
        -> Used for the FFCF calculation to     
           override the average frequency in
           in scenarios where many trajectories
           are being used in parallel.
    -zmin [value] (Not Yet Used)
        -> Set the minimum z position for molecules 
           to be included in the calculation. Based 
           on the Oxygen atom
    -zmax [value] (Not Yet Used)
        -> Set the maximum z position for molecules 
           to be included in the calculation. Based 
           on the Oxygen atom

