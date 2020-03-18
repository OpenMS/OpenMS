from String cimport *



cdef extern from "<OpenMS/CHEMISTRY/DecoyGenerator.h>" namespace "OpenMS":

    cdef cppclass DecoyGenerator:


# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/DecoyGenerator.h>" namespace "OpenMS::DecoyGenerator":
        
    # static members
    AASequence reverseProtein(const AASequence& protein) nogil except +  # wrap-attach:DecoyGenerator
    
    AASequence reversePeptides(const AASequence& protein, const String& protease) nogil except +  # wrap-attach:DecoyGenerator

    AASequence shufflePeptide(const AASequence& aas, const String& protease, const int max_attempts, int seed) nogil except +  # wrap-attach:DecoyGenerator
