from String cimport *
from AASequence cimport *



cdef extern from "<OpenMS/CHEMISTRY/DecoyGenerator.h>" namespace "OpenMS":

    cdef cppclass DecoyGenerator:
        DecoyGenerator() except + nogil 

        DecoyGenerator(DecoyGenerator &) except + nogil 

        void setSeed(UInt64) except + nogil 

        AASequence reverseProtein(const AASequence& protein) except + nogil  # wrap-doc:Reverses the protein sequence
    
        AASequence reversePeptides(const AASequence& protein, const String& protease) except + nogil  # wrap-doc:Reverses the protein's peptide sequences between enzymatic cutting positions 

        AASequence shufflePeptides(const AASequence& aas, const String& protease, const int max_attempts) except + nogil  # wrap-doc:Shuffle the protein's peptide sequences between enzymatic cutting positions, each peptide is shuffled @param max_attempts times to minimize sequence identity
