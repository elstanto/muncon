from muncon.fileio import DsdFile, SnpFile

s = SnpFile("C:\\Users\\helgr\\PycharmProjects\\muncon\\tests\\data\\test.s2p")
d = SnpFile("C:\\Users\\helgr\\PycharmProjects\\muncon\\tests\\data\\out.s2p")
s.read()
d.set_usnp(s.get_usnp())
d.write()
