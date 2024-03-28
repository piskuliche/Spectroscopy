# Computational Spectroscopy
Copyright 2024, Zeke A. Piskulich (Boston University)

Software is adapted verison of programs copyright 2023 Ashley Borkowski, Ward Thompson (University of Kansas)

## Introduction

This code is a computational spectroscopic toolkit for calculating spectroscopic observables from molecular dynamics simulations.

Thus far, the currently implemented spectroscopies are as follows:

* One Dimensional Infrared Spectroscopy
* Two Dimentional Infrared Spectroscopy
* One Dimensional Raman Spectroscopy
* One Dimensional Sum-Frequency Generation Spectroscopy

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


