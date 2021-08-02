from libcpp.vector cimport vector as libcpp_vector
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *
from ModificationDefinitionsSet cimport *

cdef extern from "<OpenMS/FORMAT/XTandemXMLFile.h>" namespace "OpenMS":

    cdef cppclass XTandemXMLFile:

        XTandemXMLFile() nogil except +
        # protected
        XTandemXMLFile(XTandemXMLFile &) nogil except + # wrap-ignore

        void load(String filename, ProteinIdentification & protein_identification,
                  libcpp_vector[PeptideIdentification] & id_data,
                  ModificationDefinitionsSet& mod_def_set) nogil except +
