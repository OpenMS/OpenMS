from Param cimport *

cdef extern from "<OpenMS/FORMAT/ParamXMLFile.h>" namespace "OpenMS":

    cdef cppclass ParamXMLFile:

        # wrap-doc:
        #  The file pendant of the Param class used to load and store the param
        #  datastructure as paramXML
        ParamXMLFile() except + nogil 
        ParamXMLFile(ParamXMLFile &) except + nogil  # compiler

        void load(String filename, Param& param) except + nogil
            # wrap-doc:
                #  Read XML file
                #  
                #  
                #  :param filename: The file from where to read the Param object
                #  :param param: The param object where the read data should be stored
                #  :raises:
                #    Exception: FileNotFound is thrown if the file could not be found
                #  :raises:
                #    Exception: ParseError is thrown if an error occurs during parsing

        void store(String filename, Param& param) except + nogil
            # wrap-doc:
                #  Write XML file
                #  
                #  
                #  :param filename: The filename where the param data structure should be stored
                #  :param param: The Param class that should be stored in the file

