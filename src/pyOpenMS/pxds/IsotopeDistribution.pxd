from libcpp cimport bool
from Types cimport *
from String cimport *
from Peak1D cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>" namespace "OpenMS":

    cdef cppclass IsotopeDistribution:        

        IsotopeDistribution() nogil except +
        IsotopeDistribution(IsotopeDistribution) nogil except + # wrap-ignore

        # overwrites the container which holds the distribution using @p distribution
        void set(libcpp_vector[ Peak1D ]& distribution) nogil except +

        void insert(double mass, float intensity) nogil except +

        libcpp_vector[Peak1D].iterator begin() nogil except +  # wrap-iter-begin:__iter__(Peak1D)
        libcpp_vector[Peak1D].iterator end()   nogil except +  # wrap-iter-end:__iter__(Peak1D)

        # returns the container which holds the distribution
        libcpp_vector[ Peak1D ]& getContainer() nogil except +

        # returns the maximal weight isotope which is stored in the distribution
        Size getMax() nogil except +

        # returns the minimal weight isotope which is stored in the distribution
        Size getMin() nogil except +

        # returns the most abundant isotope which is stored in the distribution
        Peak1D getMostAbundant() nogil except +

        # returns the size of the distribution which is the number of isotopes in the distribution
        Size size() nogil except +

        # clears the distribution and resets max isotope to 0
        void clear() nogil except +

        # renormalizes the sum of the probabilities of the isotopes to 1
        void renormalize() nogil except +

        # Trims the right side of the isotope distribution to isotopes with a significant contribution.
        void trimRight(double cutoff) nogil except +

        # Trims the left side of the isotope distribution to isotopes with a significant contribution.
        void trimLeft(double cutoff) nogil except +
        
        void merge(double, double) nogil except +

        void resize(UInt size) nogil except +
        void trimIntensities(double cutoff) nogil except +
        void sortByIntensity() nogil except +
        void sortByMass() nogil except +
        double averageMass() nogil except +

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>" namespace "OpenMS::IsotopeDistribution":
    
    cdef enum Sorted "OpenMS::IsotopeDistribution::Sorted":
        # wrap-attach:
        #    IsotopeDistribution
        INTENSITY
        MASS
        UNDEFINED

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>" namespace "OpenMS":

    cdef cppclass FineIsotopePatternGenerator:

        FineIsotopePatternGenerator() nogil except + 
        FineIsotopePatternGenerator(double threshold) nogil except +
        FineIsotopePatternGenerator(double threshold, bool use_total_prob) nogil except +
        FineIsotopePatternGenerator(double threshold, bool use_total_prob, bool absolute) nogil except +

        void setThreshold(double threshold) nogil except +
        double getThreshold() nogil except +

        void setAbsolute(bool absolute) nogil except +
        bool getAbsolute() nogil except +

        void setTotalProbability(bool total) nogil except +
        bool getTotalProbability() nogil except +

        IsotopeDistribution run(EmpiricalFormula) nogil except +

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>" namespace "OpenMS":

    cdef cppclass CoarseIsotopePatternGenerator:

        CoarseIsotopePatternGenerator() nogil except + 
        CoarseIsotopePatternGenerator(Size max_isotope) nogil except +
        CoarseIsotopePatternGenerator(Size max_isotope, bool round_masses) nogil except +

        IsotopeDistribution run(EmpiricalFormula) nogil except +

        # returns the current value of the flag to round masses to integer values (true) or return accurate masses (false)
        bool getRoundMasses() nogil except +

        # sets the round_masses_ flag to round masses to integer values (true) or return accurate masses (false)
        void setRoundMasses(bool round_masses_) nogil except +

        # returns the currently set maximum isotope
        Size getMaxIsotope() nogil except +
        
        #  sets the maximal isotope with @p max_isotope
        void setMaxIsotope(Size max_isotope) nogil except +

        # @brief Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported
        #   "Determination of Monoisotopic Masses and Ion Populations for Large Biomolecules from Resolved Isotopic Distributions"
        IsotopeDistribution estimateFromPeptideWeight(double average_weight) nogil except +

        # Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs
        IsotopeDistribution estimateFromPeptideWeightAndS(double average_weight, UInt S) nogil except +

        # Estimate Nucleotide Isotopedistribution from weight
        IsotopeDistribution estimateFromRNAWeight(double average_weight) nogil except +

        # Estimate Nucleotide Isotopedistribution from weight
        IsotopeDistribution estimateFromDNAWeight(double average_weight) nogil except +

        IsotopeDistribution estimateFromWeightAndComp(double average_weight,
                                                      double C, double H,
                                                      double N, double O,
                                                      double S, double P) nogil except +

        # Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition
        IsotopeDistribution estimateFromWeightAndCompAndS(double average_weight,
                                                          UInt S, double C,
                                                          double H, double N,
                                                          double O, double P) nogil except +

        # Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromPeptideWeight(double average_weight_precursor,
                                                                 double average_weight_fragment,
                                                                 libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        # number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
        # and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor,
                                                                     UInt S_precursor,
                                                                     double average_weight_fragment,
                                                                     UInt S_fragment,
                                                                     libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # Estimate RNA fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromRNAWeight(double average_weight_precursor,
                                                             double average_weight_fragment,
                                                             libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # Estimate DNA fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromDNAWeight(double average_weight_precursor,
                                                             double average_weight_fragment,
                                                             libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # Estimate fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, a set of isolated precursor isotopes, and average composition
        IsotopeDistribution estimateForFragmentFromWeightAndComp(double average_weight_precursor,
                                                                 double average_weight_fragment,
                                                                 libcpp_set[ unsigned int ]& precursor_isotopes,
                                                                 double C,
                                                                 double H,
                                                                 double N,
                                                                 double O,
                                                                 double S,
                                                                 double P) nogil except +

        # Calculate isotopic distribution for a fragment molecule
        IsotopeDistribution calcFragmentIsotopeDist(IsotopeDistribution& fragment_isotope_dist,
                                                    IsotopeDistribution& comp_fragment_isotope_dist,
                                                    libcpp_set[ unsigned int ]& precursor_isotopes,
                                                    double fragment_mono_mass) nogil except +

