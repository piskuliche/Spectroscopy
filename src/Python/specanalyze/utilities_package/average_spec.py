import numpy as np
import os
import pickle

def Average_Driver(basename, start_idx, stop_idx, outfile='ir_spectra.dat', spectype='1D', tag=''):
    """ This averages the spectra across directories, and then saves it to a file.

    Parameters:
    -----------
    basename : str
        The basename of the file containing the spectra.
    start_idx : int
        The first directory to read.
    stop_idx : int
        The last directory to read.
    outfile : str
        The name of the file to save the averaged spectra to.
    spectype : str
        The type of spectrum to average. Options are '1D' or '2D'.
    tag : str
        A tag to add to the beginning of the file name.

    Returns:
    --------
    spectra : tuple
        A tuple containing the frequencies, real and imaginary parts of the spectra.

    """
    if spectype == '1D':
        Average_1D_Spectra(f"{tag}{basename}", start_idx, stop_idx,
                            outfile=f"{tag}averaged_{basename}")
    elif spectype == '2D':
        Average_2D_Spectra(f"{tag}{basename}", start_idx, stop_idx,
                            outfile=f"{tag}averaged_{basename}")
    elif spectype == 'SFG':
        Average_SFG(f"{tag}{basename}", start_idx, stop_idx,
                            outfile=f"{tag}averaged_{basename}")
    else:
        raise ValueError("The valid options are '1D' or '2D'")


def Average_1D_Spectra(basename, start_idx, stop_idx, outfile='ir_spectra.dat'):
    """ This averages the spectra across directories, and then saves it to a file.


    Parameters:
    -----------
    basename : str
        The basename of the file containing the spectra.
    start_idx : int
        The first directory to read.
    stop_idx : int
        The last directory to read.
    outfile : str
        The name of the file to save the averaged spectra to.

    Returns:
    --------
    spectra : tuple
        A tuple containing the frequencies, real and imaginary parts of the spectra.

    """
    wfreq = np.genfromtxt(
        f"run_dir/{start_idx}/{basename}", usecols=0, unpack=True)
    spectral = []
    for i in range(start_idx, stop_idx):
        if i%100 == 0:
            print(i)
        temporary = np.genfromtxt(
            f"run_dir/{i}/{basename}", usecols=1, unpack=True)
        spectral.append(temporary)
    outpick = outfile.replace(".dat", ".pckl")
    pickle.dump(spectral, open(f"{outpick}", 'wb'))
    pickle.dump(wfreq, open(f"wfreq_{outpick}", 'wb'))

    averaged = np.average(spectral, axis=0)
    np.savetxt(outfile, np.c_[wfreq, averaged])
    return

def Average_2D_Spectra(basename, start_idx, stop_idx, outfile='ir_spectra.dat'):
    """ This averages the spectra across directories, and then saves it to a file.


    Parameters:
    -----------
    basename : str
        The basename of the file containing the spectra.
    start_idx : int
        The first directory to read.
    stop_idx : int
        The last directory to read.
    outfile : str
        The name of the file to save the averaged spectra to.

    Returns:
    --------
    spectra : tuple
        A tuple containing the frequencies, real and imaginary parts of the spectra.

    """
    wfreq1, wfreq2 = np.genfromtxt(
        f"run_dir/{start_idx}/{basename}", usecols=(0,1), unpack=True)
    spectral = []
    for i in range(start_idx, stop_idx):
        if i%100 == 0:
            print(i)
        temporary = np.genfromtxt(
            f"run_dir/{i}/{basename}", usecols=2, unpack=True)
        spectral.append(temporary)
    outpick = outfile.replace(".dat", ".pckl")

    averaged = np.average(spectral, axis=0)
    np.savetxt(outfile, np.c_[wfreq1,wfreq2, averaged])
    return


def Average_SFG(basename,  start_idx, stop_idx, outfile='sfg_spectrum.dat'):
    """
    This function averages the SFG spectra across directories, and then saves it to a file.

    Parameters:
    -----------
    basename : str
        The basename of the file containing the spectra.
    start_idx : int
        The first directory to read.
    stop_idx : int
        The last directory to read.
    outfile : str
        The name of the file to save the averaged spectra to.

    Returns:
    --------
    spectra : tuple
        A tuple containing the frequencies, real and imaginary parts of the spectra.


    """

    # Read in the first file to get the frequency array
    wfreq = np.genfromtxt(
        f"run_dir/{start_idx}/{basename}", usecols=0, unpack=True)
    real_part = []
    imag_part = []
    for i in range(start_idx, stop_idx):
        if i%100 == 0:
            print(i)
        real_tmp, imag_tmp = np.genfromtxt(
            "run_dir/%d/%s" % (i, basename), usecols=(1, 2), unpack=True)
        real_part.append(real_tmp)
        imag_part.append(imag_tmp)

    real_spec = np.average(real_part, axis=0)
    imag_spec = np.average(imag_part, axis=0)
    
    outpick = outfile.replace(".dat", ".pckl")
    pickle.dump(real_part, open(f"real_{outpick}", 'wb'))
    pickle.dump(imag_part, open(f"imag_{outpick}", 'wb'))
    pickle.dump(wfreq, open(f"wfreq_{outpick}", 'wb'))

    np.savetxt(outfile, np.c_[wfreq, real_spec, imag_spec])

    return (wfreq, real_spec, imag_spec)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fname', type=str, default='sfg_spectrum.dat',
                        help='Name of the file containing the spectra')
    parser.add_argument('-start', type=int, default=0,
                        help='First file to read')
    parser.add_argument('-stop', type=int, default=100,
                        help='Final file to read')
    parser.add_argument('-spectra', type=int, default=0,
                        help='[0] IR [1] Raman [2] SFG [3] 2D-IR')
    parser.add_argument('-tag', type=str, default="",
                        help='Tag to beginning of file')
    args = parser.parse_args()

    basename = args.fname
    start_idx = args.start
    stop_idx = args.stop
    spec_type = args.spectra
    tag=args.tag
    if spec_type == 0:
        for name in ["ir_spectrum.dat", "freq_dist.dat", "spec_dens.dat"]:
            Average_Spectra(f"{tag}{name}", start_idx, stop_idx,
                            outfile=f"{tag}averaged_{name}")
    elif spec_type == 1:
        for name in ["iso_spectrum.dat"]:
            Average_Spectra(f"{tag}{name}", start_idx, stop_idx,
                            outfile=f"{tag}averaged_{name}")
    elif spec_type == 2:
        Average_SFG(f"{tag}{basename}", start_idx, stop_idx)
    elif spec_type == 3:
        for name in ["2dir_imag_AAA.dat", "2dir_real_AAA.dat"]:
            waiting_times = ["00000","00200", "00400", "00600", "00800"]
            for wait in waiting_times:
                full_name = name.replace("AAA",wait)
                Average_2D_Spectra(f"{tag}{full_name}", start_idx, stop_idx,
                    outfile=f"{tag}averaged_{full_name}")
    elif spec_type == 4:
        w01=[]
        for i in range(start_idx,stop_idx):
            w01.append(np.genfromtxt(f"run_dir/{i}/{tag}w01_avg.dat"))
        avg = np.average(w01)
        print(f"average freq: {avg} cm^-1")
    elif spec_type == 5:
        name = "ffcf.dat"
        Average_Spectra(f"{tag}{name}", start_idx, stop_idx, outfile=f"{tag}averaged_{name}")
    else:
        raise ValueError("The valid options are [0-3]")

