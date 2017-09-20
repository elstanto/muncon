from muncon.fileio import DsdFile, SnpFile

s = SnpFile('/home/lstant/Documents/Stant/code/python/muncon/tests/data/test.s2p')
d = SnpFile('/home/lstant/Documents/Stant/code/python/muncon/tests/out.s2p')
s.read()
d.set_usnp(s.get_usnp())
d.write()
