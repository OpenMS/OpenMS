from Types cimport *
from libcpp cimport bool
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/CsvFile.h>" namespace "OpenMS":

  cdef cppclass CsvFile "OpenMS::CsvFile":
    CsvFile() except + nogil 
    CsvFile(CsvFile &) except + nogil  # compiler

    void load(const String& filename, char is_, bool ie_, int first_n) except + nogil  # wrap-doc:Loads data from a text file
    void store(const String& filename) except + nogil  # wrap-doc:Stores the buffer's content into a file
    void addRow(const StringList& list) except + nogil  # wrap-doc:Add a row to the buffer
    void clear() except + nogil  # wrap-doc:Clears the buffer
    bool getRow(int row, StringList& list) except + nogil  # wrap-doc:Writes all items from a row to list
