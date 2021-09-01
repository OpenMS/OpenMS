from String cimport *
from AASequence cimport *



cdef extern from "<OpenMS/CHEMISTRY/DecoyGenerator.h>" namespace "OpenMS":

    cdef cppclass DecoyGenerator:
        DecoyGenerator() nogil except +

        DecoyGenerator(DecoyGenerator &) nogil except +

        void setSeed(UInt64) nogil except +

        AASequence reverseProtein(const AASequence& protein) nogil except + # wrap-doc:Reverses the protein sequence
    
        AASequence reversePeptides(const AASequence& protein, const String& protease) nogil except + # wrap-doc:Reverses the protein's peptide sequences between enzymatic cutting positions 

        AASequence shufflePeptides(const AASequence& aas, const String& protease, const int max_attempts) nogil except + # wrap-doc:Shuffle the protein's peptide sequences between enzymatic cutting positions, each peptide is shuffled @param max_attempts times to minimize sequence identity
