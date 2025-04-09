from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/AAIndex.h>" namespace "OpenMS":
    
    cdef cppclass AAIndex "OpenMS::AAIndex":
        # private
        AAIndex() except + nogil  # wrap-ignore
        # private
        AAIndex(AAIndex &) except + nogil  # wrap-ignore
        double aliphatic(char aa) except + nogil 
        double acidic(char aa) except + nogil 
        double basic(char aa) except + nogil 
        double polar(char aa) except + nogil 
        double getKHAG800101(char aa) except + nogil 
        double getVASM830103(char aa) except + nogil 
        double getNADH010106(char aa) except + nogil 
        double getNADH010107(char aa) except + nogil 
        double getWILM950102(char aa) except + nogil 
        double getROBB760107(char aa) except + nogil 
        double getOOBM850104(char aa) except + nogil 
        double getFAUJ880111(char aa) except + nogil 
        double getFINA770101(char aa) except + nogil 
        double getARGP820102(char aa) except + nogil 
        double calculateGB(AASequence &seq, double T) except + nogil 

