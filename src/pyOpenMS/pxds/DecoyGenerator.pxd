from String cimport *
from AASequence cimport *



cdef extern from "<OpenMS/CHEMISTRY/DecoyGenerator.h>" namespace "OpenMS":

    cdef cppclass DecoyGenerator:
        void setSeed(UInt64) nogil except +

        AASequence reverseProtein(const AASequence& protein) nogil except +
    
        AASequence reversePeptides(const AASequence& protein, const String& protease) nogil except +

        AASequence shufflePeptides(const AASequence& aas, const String& protease, const int max_attempts) nogil except +
