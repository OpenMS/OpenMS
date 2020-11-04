from Types cimport *
from libcpp cimport bool
from String cimport *
from StringList cimport *

cdef extern from "<OpenMS/FORMAT/CsvFile.h>" namespace "OpenMS":

  cdef cppclass CsvFile "OpenMS::CsvFile":
    # wrap-inherits:
    #   TextFile
      CsvFile() nogil except +
      #CsvFile(const String& filename, char is, bool ie, Int first_n) nogil except + #wrap-ignore
      CsvFile(CsvFile) nogil except +

      #void load(const String& filename, char is, bool ie, Int first_n) nogil except + #wrap-ignore
      void store(const String& filename) nogil except +
      void addRow(const StringList& list) nogil except +
      void clear() nogil except +
      bool getRow(int row, StringList& list) nogil except +
