from String cimport *
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DateTime.h>" namespace "OpenMS":

    cdef cppclass DateTime:

        DateTime()   except + nogil 
        DateTime(DateTime &) except + nogil 

        void setDate(String date) except + nogil 

        void setTime(String date) except + nogil 

        # void setDate(UInt month, UInt day, UInt year) except + nogil 

        # void setTime(UInt hour, UInt minute, UInt second) except + nogil 

        # void set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second) except + nogil 

        # void get(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second) except + nogil 

        # void getDate(UInt month, UInt day, UInt year) except + nogil 

        String getDate() except + nogil 

        # void getTime(UInt hour, UInt minute, UInt second) except + nogil 

        String getTime() except + nogil 

        # Returns the current date and time
        DateTime now() except + nogil 

        #Sets the undefined date: 00/00/0000 00:00:00
        void clear() except + nogil 

        #  @brief Returns a string representation of the date and time
        #  The format of the string will be yyyy-MM-dd hh:mm:ss
        String get() except + nogil 

        #    @brief Sets date and time
        #    The following formats are supported:
        #    - MM/dd/yyyy hh:mm:ss
        #    - dd.MM.yyyy hh:mm:ss
        #    - yyyy-MM-dd hh:mm:ss
        #    - yyyy-MM-ddThh:mm:ss (ISO 8601 format)
        #    - yyyy-MM-ddZ (ISO 8601 format)
        #    - yyyy-MM-dd+hh:mm (ISO 8601 format)
        void set(String date) except + nogil 

cdef extern from "<OpenMS/DATASTRUCTURES/DateTime.h>" namespace "OpenMS::DateTime":

    DateTime now() # wrap-attach:DateTime




