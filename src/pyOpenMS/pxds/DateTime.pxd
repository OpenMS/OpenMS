from String cimport *
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DateTime.h>" namespace "OpenMS":

    cdef cppclass DateTime:

        DateTime()   nogil except +
        DateTime(DateTime) nogil except + # wrap-ignore

        void setDate(String date) nogil except +

        void setTime(String date) nogil except +

        # void setDate(UInt month, UInt day, UInt year) nogil except +

        # void setTime(UInt hour, UInt minute, UInt second) nogil except +

        # void set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second) nogil except +

        # void get(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second) nogil except +

        # void getDate(UInt month, UInt day, UInt year) nogil except +

        String getDate() nogil except +

        # void getTime(UInt hour, UInt minute, UInt second) nogil except +

        String getTime() nogil except +

        # Returns the current date and time
        DateTime now() nogil except +

        #Sets the undefined date: 00/00/0000 00:00:00
        void clear() nogil except +

        #   @brief Returns a string representation of the date and time
        #   The format of the string will be yyyy-MM-dd hh:mm:ss
        String get() nogil except +

        #     @brief Sets date and time
        #     The following formats are supported:
        #     - MM/dd/yyyy hh:mm:ss
        #     - dd.MM.yyyy hh:mm:ss
        #     - yyyy-MM-dd hh:mm:ss
        #     - yyyy-MM-ddThh:mm:ss (ISO 8601 format)
        #     - yyyy-MM-ddZ (ISO 8601 format)
        #     - yyyy-MM-dd+hh:mm (ISO 8601 format)
        void set(String date) nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/DateTime.h>" namespace "OpenMS::DateTime":

    DateTime now() # wrap-attach:DateTime




