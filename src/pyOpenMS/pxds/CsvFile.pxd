from Types cimport *
from libcpp cimport bool
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/CsvFile.h>" namespace "OpenMS":

  cdef cppclass CsvFile "OpenMS::CsvFile":
    CsvFile() nogil except +
    CsvFile(CsvFile) nogil except + #wrap-ignore

    void load(const String& filename, char is_, bool ie_, int first_n) nogil except +
    void store(const String& filename) nogil except +
    void addRow(const StringList& list) nogil except +
    void clear() nogil except +
    bool getRow(int row, StringList& list) nogil except +
