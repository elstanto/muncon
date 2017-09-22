from muncon.fileio import MeasFile, DsdFile

m = MeasFile('C:\\Users\\helgr\\PycharmProjects\\muncon\\tests\\data\\MUF\\interseries_1.meas')
m.read()
m.load_mc_samples()
m.build_mc_covariance(True)
m.mc_from_cv()
