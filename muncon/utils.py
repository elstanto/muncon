import numpy as np


class Usnp:
    """
    This class is the internal data structure used for storing s-parameter uncertainties.
    It contains the full covariance matrix.
    """

    def __init__(self):
        self.p = 0
        self.z0 = []
        self.f = []
        self.s = []
        self.v = []

    def set_ports(self, ports):
        self.p = ports
        return

    def get_ports(self):
        return self.p

    def set_z0(self, z0):
        if isinstance(z0, (list, tuple, np.ndarray)):
            if len(z0) == self.p:
                self.z0 = z0
            else:
                self.z0 = np.linspace(z0[0], z0[0], self.p)
        else:
            self.z0 = np.linspace(z0, z0, self.p)
        return
