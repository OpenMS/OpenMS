from Types cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/FORMAT/OMSSACSVFile.h>" namespace "OpenMS":
    
    cdef cppclass OMSSACSVFile "OpenMS::OMSSACSVFile":
        # wrap-doc:
                #   File adapter for OMSSACSV files
                #   -----
                #   The files contain the results of the OMSSA algorithm in a comma separated manner. This file adapter is able to
                #   load the data from such a file into the structures of OpenMS

        OMSSACSVFile() nogil except +
        OMSSACSVFile(OMSSACSVFile &) nogil except + # compiler

        void load(const String & filename,
                  ProteinIdentification & protein_identification,
                  libcpp_vector[ PeptideIdentification ] & id_data) nogil except +
            # wrap-doc:
                #   Loads a OMSSA file
                #   -----
                #   :param filename: The name of the file to read from
                #   :param protein_identification: The protein ProteinIdentification data
                #   :param id_data: The peptide ids of the file
                #   -----
                #   The content of the file is stored in `features`
                #   -----
                #   :raises:
                #     Exception: FileNotFound is thrown if the file could not be opened
                #   :raises:
                #     Exception: ParseError is thrown if an error occurs during parsing

