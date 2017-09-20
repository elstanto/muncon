from abc import ABCMeta, abstractmethod
import numpy as np
import os
from math import cos, sin, radians
import xml.etree.ElementTree as ET
from muncon.utils import Usnp, swap_s12_s21, pad_to_2port


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

    def set_usnp(self, usnp):
        self.usnp = usnp


class MeasFile(AbstractFileHandler):
    """
    This class handles file operations on NIST MUF Measurement files.
    """
    def __init__(self, path):
        super(MeasFile, self).__init__(path)
        self.saved_measpath = ""
        self.estimate = None
        self.mc_usnps = []
        self.usnp = None

    def read(self, **kwargs):
        meas_xml = ET.parse(self.path)
        meas_root = meas_xml.getroot()
        self.saved_measpath = os.path.normpath(meas_root.attrib['FileName'].replace("\\", "/"))
        estimate_path = meas_root.find(".//MeasSParams//SubItem[@Index='1']")
        estimate_path = self.localize_path(estimate_path.attrib['Text'])
        estimate = MeasSnpFile(estimate_path)
        estimate.read()
        self.estimate = estimate.get_usnp()
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

    def build_covariance(self, use_mc_mean):

        pass

    def write(self, **kwargs):
        super(MeasFile, self).write()


class SnpFile(AbstractFileHandler):
    """
    This class handles file operations on s-parameter files.
    """
    f_unit_dict = {
        'Hz': 1,
        'KHz': 1e3,
        'MHz': 1e6,
        'GHz': 1e9
    }

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
        self.f_multiplier = SnpFile.f_unit_dict[f_units]
        self.format = lineparts[3]
        z0 = float(lineparts[5])
        self.usnp.set_z0(z0)

    def write(self, **kwargs):
        super(SnpFile, self).write()
        usnp = self.get_usnp()
        sparams = usnp.get_sparams()
        self.f_multiplier = 1  # Hz
        self.format = "RI"
        self.build_column_header()
        if usnp.get_ports() == 1:
            sparams = pad_to_2port(sparams)
        with open(self.path, 'w') as f:
            for header_line in self.get_usnp().get_comments():
                f.write(header_line)
            f.write(self.column_header)
            for ifreq in range(1, len(usnp.get_freqs())):
                f.write(str(usnp.get_freqs()[ifreq] / self.f_multiplier) + " ")
                f.write(" ".join(" ".join([str(np.real(e)), str(np.imag(e))]) for e in sparams[ifreq]))
                f.write("\n")

    def build_column_header(self):
        units = list(SnpFile.f_unit_dict.keys())[list(SnpFile.f_unit_dict.values()).index(self.f_multiplier)]
        self.column_header = " ".join(["#", units, "S", self.format, "R", str(int(self.get_usnp().get_z0()[0]))]) + "\n"


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
        title = "MMS4 Definition Standard Data\n"
        filename_header = "Filename: " + self.path + "\n"
        header = "This file contains the standard definition S-parameters and their uncertainty covariance matrix.\n" \
                 "A comment line begins with the exclamation mark '!'.\n" \
                 "Each frequency corresponds to a line with 74 columns, ordered as below:\n" \
                 "freq(GHz) Nports S11r S11i S12r S12i S21r S21i S22r S22i (continues below)\n" \
                 "V11,11rr V11,11ri V11,11ir V11,11ii V11,12rr V11,12ri V11,12ir V11,12ii V11,21rr V11,21ri V11,21ir " \
                 "V11,21ii V11,22rr V11,22ri V11,22ir V11,22ii (cont.)\n" \
                 "V12,11rr V12,11ri V12,11ir V12,11ii V12,12rr V12,12ri V12,12ir V12,12ii V12,21rr V12,21ri V12,21ir " \
                 "V12,21ii V12,22rr V12,22ri V12,22ir V12,22ii (cont.)\n" \
                 "V21,11rr V21,11ri V21,11ir V21,11ii V21,12rr V21,12ri V21,12ir V21,12ii V21,21rr V21,21ri V21,21ir " \
                 "V21,21ii V21,22rr V21,22ri V21,22ir V21,22ii (cont.)\n" \
                 "V22,11rr V22,11ri V22,11ir V22,11ii V22,12rr V22,12ri V22,12ir V22,12ii V22,21rr V22,21ri V22,21ir " \
                 "V22,21ii V22,22rr V22,22ri V22,22ir V22,22ii (cont.)\n" \
                 "where 'r' stands for real and 'i' for imaginary parts\n" \
                 "One-port standards (Nports=1) have zeros in 12, 21 and 22 parameters\n"
        usnp = self.get_usnp()
        sparams = usnp.get_sparams()
        if usnp.get_ports() == 1:
            sparams = pad_to_2port(sparams)
        sparams = swap_s12_s21(sparams)
        with open(self.path, 'w') as f:
            f.write(title)
            f.write(filename_header)
            f.write(header)
            for ifreq in range(1, len(usnp.get_freqs())):
                f.write(str(usnp.get_freqs()[ifreq]/1e9) + " ")
                f.write(" ".join(" ".join([str(np.real(e)), str(np.imag(e))]) for e in sparams[ifreq]))
                f.write("\n")


