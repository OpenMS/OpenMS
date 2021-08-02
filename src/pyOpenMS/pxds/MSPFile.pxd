from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from Peak1D cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/MSPFile.h>" namespace "OpenMS":

    cdef cppclass MSPFile:

        MSPFile() nogil except + # wrap-doc:File adapter for MSP files (NIST spectra library)
        MSPFile(MSPFile &) nogil except +

        void store(String filename, MSExperiment & exp) nogil except + # wrap-doc:Stores a map in a MSPFile file
        void load(String filename, libcpp_vector[PeptideIdentification] & ids, MSExperiment & exp) nogil except +
            # wrap-doc:
                #   Loads a map from a MSPFile file
                #   -----
                #   :param exp: PeakMap which contains the spectra after reading
                #   :param filename: The filename of the experiment
                #   :param ids: Output parameter which contains the peptide identifications from the spectra annotations
