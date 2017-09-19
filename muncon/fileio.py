from abc import ABCMeta, abstractmethod
import numpy as np
import os
from math import cos, sin, radians
import xml.etree.ElementTree as ET
from muncon.utils import Usnp


class AbstractFileHandler:
    __metaclass__ = ABCMeta

    def __init__(self, path):
        self.usnp = Usnp()
        self.path = path

    @abstractmethod
    def read(self, **kwargs):
        pass

    @abstractmethod
    def write(self, **kwargs):
        pass

    def get_usnp(self):
        return self.usnp


class MeasFile(AbstractFileHandler):
    """
    This class handles file operations on NIST MUF Measurement files.
    """
    def __init__(self, path):
        super(MeasFile, self).__init__(path)
        self.use_MC_mean = False
        self.saved_measpath = ""
        self.mc_usnps = []

    def read(self, **kwargs):
        meas_xml = ET.parse(self.path)
        meas_root = meas_xml.getroot()
        self.saved_measpath = os.path.normpath(meas_root.attrib['FileName'].replace("\\", "/"))
        for element in meas_root.findall(".//MonteCarloPerturbedSParams//SubItem[@Index='1']"):
            samplepath = self.localize_path(element.attrib['Text'])
            sample = MeasSnpFile(samplepath)
            sample.read()
            self.mc_usnps.append(sample.get_usnp())
        super(MeasFile, self).read()

    def localize_path(self, path):
        samplepath = os.path.normpath(path.replace("\\", "/"))
        saved_base = os.path.commonpath([self.saved_measpath, samplepath])
        meas_suffix = self.saved_measpath.replace(saved_base, "")
        local_prefix = self.path.replace(meas_suffix, "")
        local_samplepath = samplepath.replace(saved_base, local_prefix)
        return local_samplepath

    def write(self, **kwargs):
        super(MeasFile, self).write()


class SnpFile(AbstractFileHandler):
    """
    This class handles file operations on s-parameter files.
    """
    def __init__(self, path):
        super(SnpFile, self).__init__(path)
        self.column_header = ""
        self.format = ""
        self.f_multiplier = 1
        self.usnp.set_ports(int(path.split(".")[-1][1]))

    def read(self, **kwargs):
        super(SnpFile, self).read()
        with open(self.path, "r") as f:
            freqs = []
            sparams_real = []
            sparams_imag = []
            while True:
                line = f.readline()
                if not line:
                    break
                elif line[0] == "!":
                    self.usnp.add_comment(line)
                elif line[0] == "#":
                    self.column_header = line
                    self.parse_column_header()
                else:
                    lineparts = [float(i) for i in line.split()]
                    freqs.append(lineparts[0]*self.f_multiplier)
                    sparams = []
                    for i in range(1, int((len(lineparts)-1)/2+1)):
                        sparams.extend(self.get_ri_sparam([lineparts[2*i-1], lineparts[2*i]]))
                    sparams_real.append(sparams[0:7:2])
                    sparams_imag.append(sparams[1:8:2])
            self.usnp.set_freqs(np.asarray(freqs))
            self.usnp.set_sparams(np.vectorize(complex)(sparams_real, sparams_imag))
            cov_zeros = np.zeros((len(freqs), self.usnp.get_ports()**2, self.usnp.get_ports()**2))
            self.usnp.set_covariance(np.vectorize(complex)(cov_zeros, cov_zeros))

    def get_ri_sparam(self, sparam):
        if self.format == "MA":
            ri_sparam = [sparam[0]*cos(radians(sparam[1])), sparam[0]*sin(radians(sparam[1]))]
            return ri_sparam
        elif self.format == "DB":
            linear_mag = 10**(sparam[0]/20)
            ri_sparam = [linear_mag * cos(radians(sparam[1])), linear_mag * sin(radians(sparam[1]))]
            return ri_sparam
        else:  # Handles RI and no column header
            return sparam

    def parse_column_header(self):
        lineparts = self.column_header.split()
        f_units = lineparts[1]
        f_unit_dict = {
            'Hz': 1,
            'KHz': 1e3,
            'MHz': 1e6,
            'GHz': 1e9
        }
        self.f_multiplier = f_unit_dict[f_units]
        self.format = lineparts[3]
        z0 = float(lineparts[5])
        self.usnp.set_z0(z0)

        return

    def write(self, **kwargs):
        super(SnpFile, self).write()
        return


class MeasSnpFile(SnpFile):
    """
    This derived class sets some default properties which are missing from the MUF snp file as it can contain no header.
    """

    def __init__(self, path):
        super(MeasSnpFile, self).__init__(path)

    def read(self, **kwargs):
        self.format = "RI"
        self.f_multiplier = 1e9  # GHz
        self.usnp.set_z0(50)
        super(MeasSnpFile, self).read()
        self.usnp.set_covariance(None)  # Saves much unneeded space

    def write(self, **kwargs):
        super(MeasSnpFile, self).write()


class DsdFile(AbstractFileHandler):
    """
    This class handles file operations on Keysight DSD files.
    """
    def __init__(self, path):
        super(DsdFile, self).__init__(path)

    def read(self, **kwargs):
        super(DsdFile, self).read()

    def write(self, **kwargs):
        super(DsdFile, self).write()
