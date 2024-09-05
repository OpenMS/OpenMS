from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from ProteinIdentification cimport *
from PeptideIdentification cimport *

from String cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FORMAT/MzIdentMLFile.h>" namespace "OpenMS":

    cdef cppclass MzIdentMLFile(ProgressLogger):
        # wrap-inherits:
        #  ProgressLogger

        MzIdentMLFile() except + nogil 
        MzIdentMLFile(MzIdentMLFile &) except + nogil 

        void load(String filename, libcpp_vector[ProteinIdentification] & poid, libcpp_vector[PeptideIdentification] & peid) except + nogil 
            # wrap-doc:
                #  Loads the identifications from a MzIdentML file
                #  
                #  
                #  :param filename: File name of the file to be checked
                #  :raises:
                #    Exception: FileNotFound is thrown if the file could not be opened
                #  :raises:
                #    Exception: ParseError is thrown if an error occurs during parsin

        void store(String filename, libcpp_vector[ProteinIdentification] & poid, libcpp_vector[PeptideIdentification] & peid) except + nogil 
            # wrap-doc:
                #  Stores the identifications in a MzIdentML file
                #  
                #  
                #  :raises:
                #    Exception: UnableToCreateFile is thrown if the file could not be created

        bool isSemanticallyValid(String filename, StringList errors, StringList warnings) except + nogil 
            # wrap-doc:
                #  Checks if a file is valid with respect to the mapping file and the controlled vocabulary
                #  
                #  
                #  :param filename: File name of the file to be checked
                #  :param errors: Errors during the validation are returned in this output parameter
                #  :param warnings: Warnings during the validation are returned in this output parameter
                #  :raises:
                #    Exception: FileNotFound is thrown if the file could not be opened

