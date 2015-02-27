from Types cimport *
from MzTab cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MzTabFile.h>" namespace "OpenMS":

    cdef cppclass MzTabFile:

        MzTabFile() nogil except +

        void store(String filename, MzTab & mz_tab) nogil except +
        void load(String filename, MzTab & mz_tab) nogil except +

        # Does not exist
        # void storeProteinReliabilityColumn(bool store) nogil except +
        # void storePeptideReliabilityColumn(bool store) nogil except +
        # void storePSMReliabilityColumn(bool store) nogil except +
        # void storeSmallMoleculeReliabilityColumn(bool store) nogil except +
        # void storeProteinUriColumn(bool store) nogil except +
        # void storePeptideUriColumn(bool store) nogil except +
        # void storePSMUriColumn(bool store) nogil except +
        # void storeSmallMoleculeUriColumn(bool store) nogil except +
        # void storeProteinGoTerms(bool store) nogil except +

