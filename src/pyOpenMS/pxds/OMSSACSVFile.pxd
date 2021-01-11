from Types cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/FORMAT/OMSSACSVFile.h>" namespace "OpenMS":
    
    cdef cppclass OMSSACSVFile "OpenMS::OMSSACSVFile":
        OMSSACSVFile() nogil except +
        OMSSACSVFile(OMSSACSVFile) nogil except + #wrap-ignore

        void load(const String & filename,
                  ProteinIdentification & protein_identification,
                  libcpp_vector[ PeptideIdentification ] & id_data) nogil except +


