from libcpp.vector cimport vector as libcpp_vector
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

from libcpp cimport bool

cdef extern from "<OpenMS/FORMAT/ProtXMLFile.h>" namespace "OpenMS":

    cdef cppclass ProtXMLFile:

        ProtXMLFile() nogil except +
        # copy constructor of 'ProtXMLFile' is implicitly deleted because base class 'Internal::XMLHandler' has a deleted copy constructor protected Internal::XMLHandler
        ProtXMLFile(ProtXMLFile &) nogil except + # wrap-ignore

        void load(String filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids) nogil except +

        # Not implemented
        void store(String filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids, String document_id) nogil except +

