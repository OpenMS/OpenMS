from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/TextFile.h>" namespace "OpenMS":
    
    cdef cppclass TextFile "OpenMS::TextFile":
        TextFile() except + nogil  # wrap-doc:This class provides some basic file handling methods for text files
        TextFile(TextFile &) except + nogil  # compiler
        TextFile(const String &filename, bool trim_linesalse, Int first_n1) except + nogil 
        void load(const String &filename, bool trim_linesalse, Int first_n1) except + nogil  
        # wrap-doc:
                #  Loads data from a text file
                #  
                #  :param filename: The input file name
                #  :param trim_lines: Whether or not the lines are trimmed when reading them from file
                #  :param first_n: If set, only `first_n` lines the lines from the beginning of the file are read
                #  :param skip_empty_lines: Should empty lines be skipped? If used in conjunction with `trim_lines`, also lines with only whitespace will be skipped. Skipped lines do not count towards the total number of read lines

        void store(const String &filename) except + nogil  # wrap-doc:Writes the data to a file
        void addLine(const String line) except + nogil 
