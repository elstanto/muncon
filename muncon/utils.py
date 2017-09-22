from numpy import ndarray, linspace, zeros, vectorize, sin, cos, radians


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

    def get_sparams(self):
        return self.s

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

    @staticmethod
    def pad_to_2port(sparams):
        sparams_2p = zeros((len(sparams), 8))
        sparams_2p[:, 0:2] = sparams[:, 0:2]
        return sparams_2p

    @staticmethod
    def swap_s12_s21(sparams):
        # s21 = sparams[:, 2:4].copy()
        # sparams[:, 2:4] = sparams[:, 4:6]  # S21 > S12
        # sparams[:, 4:6] = s21  # S12 > S21
        sparams[:, [2, 4]] = sparams[:, [4, 2]]
        sparams[:, [3, 5]] = sparams[:, [5, 3]]
        return sparams

    @staticmethod
    def swap_v12_v21(covariance):
        covariance[:, [2, 4]] = covariance[:, [4, 2]]
        covariance[:, [3, 5]] = covariance[:, [5, 3]]
        covariance[[2, 4], :] = covariance[[4, 2], :]
        covariance[[3, 5], :] = covariance[[5, 3], :]
        return covariance

    def format_sparam(sparam, old_format, new_format):
        #TODO
        if old_format == new_format:
            return sparam
        elif old_format == "RI":
            if new_format == "MA":
                return
            elif new_format == "DB":
                return
        elif old_format == "MA":
            if new_format == "RI":
                ri_sparam = [sparam[0]*cos(radians(sparam[1])), sparam[0]*sin(radians(sparam[1]))]
                return ri_sparam
            elif new_format == "DB":
                return
        elif old_format == "DB":
            if new_format == "RI":
                linear_mag = 10**(sparam[0]/20)
                ri_sparam = [linear_mag * cos(radians(sparam[1])), linear_mag * sin(radians(sparam[1]))]
                return ri_sparam
            elif new_format == "MA":
                return
        else:  # Handles no column header
            return sparam
