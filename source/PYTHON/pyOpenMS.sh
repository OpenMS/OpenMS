#!/share/opt/bin/python2.3 -i
import sys
sys.path.append("/home/kohlbach/Projects/OpenMS/open-ms/lib/Linux-1386-g++_3.2.2/:/share/opt/qt-x11-free-3.3.3/lib")
print sys.path
from OpenMS import *
print VersionInfo.getVersion()
f = DTAFile("test.dta")
spec = DDiscreteSpectrum1()
f >> spec
print "Number of Peaks:", spec.size()
