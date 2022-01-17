from Types cimport *
from libcpp cimport bool
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/CsvFile.h>" namespace "OpenMS":

  cdef cppclass CsvFile "OpenMS::CsvFile":
    CsvFile() nogil except +
    CsvFile(CsvFile &) nogil except + # compiler

    void load(const String& filename, char is_, bool ie_, int first_n) nogil except + # wrap-doc:Loads data from a text file
    void store(const String& filename) nogil except + # wrap-doc:Stores the buffer's content into a file
    void addRow(const StringList& list) nogil except + # wrap-doc:Add a row to the buffer
    void clear() nogil except + # wrap-doc:Clears the buffer
    bool getRow(int row, StringList& list) nogil except + # wrap-doc:Writes all items from a row to list
