from String cimport *
from StringList cimport *
from TargetedExperiment cimport *

cdef extern from "<OpenMS/FORMAT/TraMLFile.h>" namespace "OpenMS":

    cdef cppclass TraMLFile:

        TraMLFile() except + nogil 
        TraMLFile(TraMLFile &) except + nogil  # compiler

        void load(String filename,
                  TargetedExperiment & id) except + nogil  # wrap-doc:Loads a map from a TraML file

        void store(String filename,
                  TargetedExperiment & id) except + nogil  # wrap-doc:Stores a map in a TraML file

        bool isSemanticallyValid(String filename, StringList & errors,
                                 StringList & warnings) except + nogil 
        # wrap-doc:
                #  Checks if a file is valid with respect to the mapping file and the controlled vocabulary
                #  
                #  :param filename: File name of the file to be checked
                #  :param errors: Errors during the validation are returned in this output parameter
                #  :param warnings: Warnings during the validation are returned in this output parameter
