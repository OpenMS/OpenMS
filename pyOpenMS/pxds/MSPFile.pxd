from String cimport *

from MSSpectrum cimport *
from Peak1D cimport *
from FeatureMap cimport *
from Feature cimport *

cdef extern from "<OpenMS/FORMAT/MSPFile.h>" namespace "OpenMS":

    cdef cppclass MSPFile:

        MSPFile() nogil except +

        # TODO RichPeakMaps are not really supported ...
        # void store(String filename, RichPeakMap & exp)
        # void load(String filename, std::vector<PeptideIdentification> & ids, RichPeakMap & exp)
