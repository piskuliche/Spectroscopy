import h5py
import numpy as np


class AverageFreq:
    def __init__(self, start_idx=0, stop_idx=1000, mapfile="ir_base/empirical_map.in"):
        self.start_idx = 0
        self.stop_idx = 0
        self.mapfile = mapfile
        self.cmiper_au = 2.1947463 * 10 ** 5  # cm^-1 to atomic units

        # Initialize Storage
        self.freqs = None
        self.map_coeffs = None

    def read_freq_map(self):
        """ The frequency map is located on the second line of the mapfile.

        The frequency map is a quadratic map, and the line stores that coefficients
        of the quadratic function. The coefficients are stored in the following order:
        [a, b, c] where the quadratic function is f(x) = a + bx + cx^2.
        """
        coeffs = None
        with open(self.mapfile, "r") as f:
            f.readline()
            line = f.readline().strip()
            coeffs = [float(x) for x in line.split()]
        self.map_coeffs = coeffs

    def read_field_file(self, filename):
        """ Read the frequencies from the HDF5 file.

        The frequencies are stored in the HDF5 file under the dataset "frequencies".
        """
        fdata = []
        with h5py.File(filename, 'r') as data:
            fdata.extend(data[key][()] for key in data.keys() if "dot" in key)
        fdata = np.array(fdata)
        freqs = self.quad_map(fdata, *self.map_coeffs)*self.cmiperau
        return fdata.flatten()

    def generate_av_frequency(self):
        """ Generate the average frequency from the frequency map. """
        for i in range(self.start_idx, self.stop_idx):
            freqs = self.read_field_file(f"run_dir/{i}/field.h5")
            av_freqs = np.average(freqs)
            self.freqs.append(av_freqs)
        return

    @staticmethod
    def quad_map(x, a, b, c):
        return a + b * x + c * x ** 2.0


if __name__ == "__main__":

    av = AverageFreq(mapfile='empirical_map.in')
    av.read_freq_map()
    freqs = av.read_field_file('field.h5')
    print(np.average(freqs))
