from Types cimport *
from MzTab cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MzTabFile.h>" namespace "OpenMS":

    cdef cppclass MzTabFile:

        MzTabFile() except + nogil 
        MzTabFile(MzTabFile &) except + nogil  # compiler

        void store(String filename, MzTab & mz_tab) except + nogil  # wrap-doc:Stores MzTab file
        void load(String filename, MzTab & mz_tab) except + nogil  # wrap-doc:Loads MzTab file

        # Does not exist
        # void storeProteinReliabilityColumn(bool store) except + nogil 
        # void storePeptideReliabilityColumn(bool store) except + nogil 
        # void storePSMReliabilityColumn(bool store) except + nogil 
        # void storeSmallMoleculeReliabilityColumn(bool store) except + nogil 
        # void storeProteinUriColumn(bool store) except + nogil 
        # void storePeptideUriColumn(bool store) except + nogil 
        # void storePSMUriColumn(bool store) except + nogil 
        # void storeSmallMoleculeUriColumn(bool store) except + nogil 
        # void storeProteinGoTerms(bool store) except + nogil 

