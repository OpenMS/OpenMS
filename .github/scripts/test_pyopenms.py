import pyopenms

seq = pyopenms.AASequence.fromString("DFPIANGER")
print(seq.toString())
seq.getMonoWeight()
