from numpy import ndarray, linspace


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
        self.comments = []

    def set_ports(self, ports):
        self.p = ports

    def get_ports(self):
        return self.p

    def set_z0(self, z0):
        if isinstance(z0, (list, tuple, ndarray)):
            if len(z0) == self.p:
                self.z0 = z0
            else:
                self.z0 = linspace(z0[0], z0[0], self.p)
        else:
            self.z0 = linspace(z0, z0, self.p)

    def get_z0(self):
        return self.z0

    def set_freqs(self, freqs):
        self.f = freqs

    def get_freqs(self):
        return self.f

    def set_sparams(self, sparams):
        self.s = sparams

    def set_covariance(self, covariance):
        self.v = covariance

    def get_covariance(self):
        return self.v

    def set_comments(self, comments):
        self.comments = comments

    def get_comments(self):
        return self.comments

    def add_comment(self, comment):
        self.comments.append(comment)
