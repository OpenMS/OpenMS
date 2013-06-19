from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/AAIndex.h>" namespace "OpenMS":
    
    cdef cppclass AAIndex "OpenMS::AAIndex":
        AAIndex(AAIndex) nogil except + #wrap-ignore
        DoubleReal aliphatic(char aa) nogil except +
        DoubleReal acidic(char aa) nogil except +
        DoubleReal basic(char aa) nogil except +
        DoubleReal polar(char aa) nogil except +
        DoubleReal getKHAG800101(char aa) nogil except +
        DoubleReal getVASM830103(char aa) nogil except +
        DoubleReal getNADH010106(char aa) nogil except +
        DoubleReal getNADH010107(char aa) nogil except +
        DoubleReal getWILM950102(char aa) nogil except +
        DoubleReal getROBB760107(char aa) nogil except +
        DoubleReal getOOBM850104(char aa) nogil except +
        DoubleReal getFAUJ880111(char aa) nogil except +
        DoubleReal getFINA770101(char aa) nogil except +
        DoubleReal getARGP820102(char aa) nogil except +
        DoubleReal calculateGB(AASequence &seq, DoubleReal T) nogil except +

