from muncon.fileio import DsdFile, SnpFile

s = SnpFile("C:\\Users\\helgr\\PycharmProjects\\muncon\\tests\\data\\test.s2p")
d = DsdFile("C:\\Users\\helgr\\PycharmProjects\\muncon\\tests\\data\\out.dsd")
s.read()
d.set_usnp(s.get_usnp())
d.write()
