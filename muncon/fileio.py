from abc import ABCMeta, abstractmethod
from muncon.utils import Usnp


class AbstractFileHandler:
    __metaclass__ = ABCMeta

    def __init__(self, path):
        self.usnp = Usnp()
        self.path = path
        return

    @abstractmethod
    def read(self, **kwargs):
        return

    @abstractmethod
    def write(self, **kwargs):
        return


class MeasFile(AbstractFileHandler):
    """
    This class handles file operations on NIST MUF Measurement files.
    """
    def __init__(self, path):
        super(MeasFile, self).__init__(path)
        return

    def read(self, **kwargs):
        super(MeasFile, self).read()
        return

    def write(self, **kwargs):
        super(MeasFile, self).write()
        return


class SnpFile(AbstractFileHandler):
    """
    This class handles file operations on s-parameter files.
    """
    def __init__(self, path):
        super(SnpFile, self).__init__(path)
        self.column_header = ""
        self.format = ""
        self.f_multiplier = 1
        return

    def read(self, **kwargs):
        super(SnpFile, self).read()
        with open(self.path, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if line[0] == "!": #Comment
                    self.usnp.add_comment(line)
                if line[0] == "#": #Column header
                    self.column_header = line
                    self.parse_column_header()
        return

    def parse_column_header(self):
        lineparts = self.column_header.split()
        f_units = lineparts[2]
        f_unit_dict = {
            'Hz' : 1,
            'KHz' : 1e3,
            'MHz' : 1e6,
            'GHz' : 1e9
        }
        self.f_multiplier = f_unit_dict(f_units)
        self.format = lineparts[4]
        z0 = float(lineparts[6])
        self.usnp.set_z0(z0)

        return

    def write(self, **kwargs):
        super(SnpFile, self).write()
        return


class DsdFile(AbstractFileHandler):
    """
    This class handles file operations on Keysight DSD files.
    """
    def __init__(self, path):
        super(DsdFile, self).__init__(path)
        return

    def read(self, **kwargs):
        super(DsdFile, self).read()
        return

    def write(self, **kwargs):
        super(DsdFile, self).write()
        return
