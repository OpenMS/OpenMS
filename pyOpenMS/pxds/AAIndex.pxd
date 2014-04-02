from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/AAIndex.h>" namespace "OpenMS":
    
    cdef cppclass AAIndex "OpenMS::AAIndex":
        AAIndex(AAIndex) nogil except + #wrap-ignore
        double aliphatic(char aa) nogil except +
        double acidic(char aa) nogil except +
        double basic(char aa) nogil except +
        double polar(char aa) nogil except +
        double getKHAG800101(char aa) nogil except +
        double getVASM830103(char aa) nogil except +
        double getNADH010106(char aa) nogil except +
        double getNADH010107(char aa) nogil except +
        double getWILM950102(char aa) nogil except +
        double getROBB760107(char aa) nogil except +
        double getOOBM850104(char aa) nogil except +
        double getFAUJ880111(char aa) nogil except +
        double getFINA770101(char aa) nogil except +
        double getARGP820102(char aa) nogil except +
        double calculateGB(AASequence &seq, double T) nogil except +

