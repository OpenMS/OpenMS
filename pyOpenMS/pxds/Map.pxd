from libcpp cimport bool

cdef extern from "<OpenMS/DATASTRUCTURES/Map.h>" namespace "OpenMS":

    cdef cppclass Map[U, V]:
        # we have converters for this, do not wrap the class itself !
        # wrap-ignore

        U first
        V second

        Map()
        Map(Map[U,V])

        V& operator[](U&)  nogil

        cppclass iterator:
            bool operator!=(Map.iterator) nogil
            Map operator*() nogil
            Map operator++() nogil

        Map.iterator begin() nogil
        Map.iterator end() nogil
