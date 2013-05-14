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
        ProtXMLFile(ProtXMLFile) nogil except + # wrap-ignore

        void load(String filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids)
        void store(String filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids, String document_id)

