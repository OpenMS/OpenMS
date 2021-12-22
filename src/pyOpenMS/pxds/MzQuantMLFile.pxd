from MSQuantifications cimport *
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/MzQuantMLFile.h>" namespace "OpenMS":

    cdef cppclass MzQuantMLFile:
        MzQuantMLFile() nogil except +
        MzQuantMLFile(MzQuantMLFile &) nogil except +

        void load(String filename, MSQuantifications & msq) nogil except +
            # wrap-doc:
                #   Loads a map from a MzQuantML file
                #   -----
                #   :raises:
                #     Exception: FileNotFound is thrown if the file could not be opened
                #   :raises:
                #     Exception: ParseError is thrown if an error occurs during parsing

        void store(String filename, MSQuantifications & msq) nogil except +
            # wrap-doc:
                #   Stores a map in a MzQuantML file
                #   -----
                #   :raises:
                #     Exception: UnableToCreateFile is thrown if the file could not be created

        bool isSemanticallyValid(String filename,
                                 StringList & errors,
                                 StringList & warnings) nogil except +
            # wrap-doc:
                #   Checks if a file is valid with respect to the mapping file and the controlled vocabulary
                #   -----
                #   :param filename: File name of the file to be checked
                #   :param errors: Errors during the validation are returned in this output parameter
                #   :param warnings: Warnings during the validation are returned in this output parameter
                #   -----
                #   :raises:
                #     Exception: UnableToCreateFile is thrown if the file could not be created


