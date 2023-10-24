import numpy as np


if __name__ == "__main__":

    import argparse
    import os

    parser = argparse.ArgumentParser(description='Average directories')
    parser.add_argument('-f', type=str, help='2dir_real_00000.dat')
    parser.add_argument('-nstart', type=int, help='Start directory')
    parser.add_argument('-nstop',  type=int, help='Stop directory')
    parser.add_argument('-nstep',  type=int, help='Step directory')
    args = parser.parse_args()

    filename = args.f
    nstart = args.nstart
    nstop = args.nstop
    nstep = args.nstep

    freq1, freq3 = np.genfromtxt("run_dir/" + str(nstart)+'/'+filename,
                                 dtype=float, usecols=(0, 1), autostrip=True, unpack=True)

    ir_spectra = []
    for i in range(nstart, nstop, nstep):
        ir_tmp = np.genfromtxt("run_dir/"+str(nstart)+'/'+filename,
                               dtype=float, usecols=(2), autostrip=True, unpack=True)
        ir_spectra.append(ir_tmp)
    ir = np.array(ir_spectra)
    ir_avg = np.average(ir, axis=0)
    os.system("mkdir -p average")
    np.savetxt("average/"+filename, np.transpose([freq1, freq3, ir_avg]), fmt=[
               '%13.8f', '%13.8f', '%13.8f'], delimiter=' ')
