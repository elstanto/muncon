from muncon.fileio import MeasFile, DsdFile, SnpFile
import numpy as np

# s = SnpFile('C:\\Users\\helgr\\PycharmProjects\\muncon\\tests\\data\\out.s2p')
# s.read()
# d = DsdFile('C:\\Users\\helgr\\PycharmProjects\\muncon\\tests\\data\\out.dsd')
# d.set_usnp(s.get_usnp())
# d.write()
#
m = MeasFile('C:\\Users\\helgr\\PycharmProjects\\muncon\\tests\\data\\MUF\\interseries_1.meas')
m.read()
m.load_cv_samples()
m.build_cv_covariance()
m.mc_from_cv()
