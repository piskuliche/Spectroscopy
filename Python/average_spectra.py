"""
Program: average_spectra.py

This program is aimed at averaging spectra for calculations with done with a trajectory launching spectral approach.

Copyright, Zeke Piskulich, 2023. 

"""

import numpy as np
import os


def Average_SFG(basename,  start_idx, stop_idx, outfile='sfg_spectra.dat'):
    """
    This function averages the SFG spectra across directories, and then saves it to a file.

    Args:
        basename (str): The basename of the file containing the spectra.
        start_idx (int): The first directory to read.
        stop_idx (int): The last directory to read.
        outfile (str): The name of the file to save the averaged spectra to.

    Returns:
        spectra (tuple): A tuple containing the frequencies, real and imaginary parts of the spectra.

    Outputs:
        A file containing the averaged spectra.
    """

    # Read in the first file to get the frequency array
    wfreq = np.genfromtxt("run_dir/0/%s" % basename, usecols=0, unpack=True)
    real_part = []
    imag_part = []
    for i in range(start_idx, stop_idx):
        real_tmp, imag_tmp = np.genfromtxt(
            "run_dir/%d/%s" % (i, basename), usecols=(1, 2), unpack=True)
        real_part.append(real_tmp)
        imag_part.append(imag_tmp)

    real_spec = np.average(real_part, axis=0)
    imag_spec = np.average(imag_part, axis=0)

    np.savetxt(outfile, np.c_[wfreq, real_spec, imag_spec])

    return (wfreq, real_spec, imag_spec)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fname', type=str, default='sfg_spectra.dat',
                        help='Name of the file containing the spectra')
    parser.add_argument('-start', type=int, default=0,
                        help='First file to read')
    parser.add_argument('-stop', type=int, default=100,
                        help='Final file to read')
    parser.add_argument('-spectra', type=int, default=0,
                        help='[0] IR [1] Raman [2] SFG [3] 2D-IR')
    args = parser.parse_args()

    basename = args.fname
    start_idx = args.start
    stop_idx = args.stop
    spec_type = args.spectra

    if (spec_type == 2):
        Average_SFG(basename, start_idx, stop_idx)
