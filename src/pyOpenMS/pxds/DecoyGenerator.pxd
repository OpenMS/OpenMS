from String cimport *
from AASequence cimport *
from libcpp.vector cimport vector as libcpp_vector



cdef extern from "<OpenMS/CHEMISTRY/DecoyGenerator.h>" namespace "OpenMS":

    cdef cppclass DecoyGenerator:
        DecoyGenerator() except + nogil 

        DecoyGenerator(DecoyGenerator &) except + nogil 

        void setSeed(UInt64) except + nogil 

        AASequence reverseProtein(const AASequence& protein) except + nogil  # wrap-doc:Reverses the protein sequence
    
        AASequence reversePeptides(const AASequence& protein, const String& protease) except + nogil  # wrap-doc:Reverses the protein's peptide sequences between enzymatic cutting positions 

        AASequence shufflePeptides(const AASequence& aas, const String& protease, const int max_attempts) except + nogil  # wrap-doc:Shuffle the protein's peptide sequences between enzymatic cutting positions, each peptide is shuffled @param max_attempts times to minimize sequence identity

        libcpp_vector[AASequence] shuffle(const AASequence& protein, const String& protease, int decoy_factor) except + nogil  # wrap-doc:Generate decoy protein sequences using shuffle algorithm. Digests protein using specified protease and shuffles each peptide. For top-down proteomics use "no cleavage". decoy_factor is the number of complete decoy proteins to generate. Returns vector of AASequence
