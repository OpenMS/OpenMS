from libcpp cimport bool
from Types cimport *
from String cimport *
from Peak1D cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>" namespace "OpenMS":

    cdef cppclass IsotopeDistribution:        

        IsotopeDistribution() nogil except +
            # wrap-doc:
            #   Isotope distribution class
            #   -----
            #   A container that holds an isotope distribution. It consists of mass values
            #   and their correspondent probabilities (stored in the intensity slot)
            #   -----
            #   Isotope distributions can be calculated using either the
            #   CoarseIsotopePatternGenerator for quantized atomic masses which group
            #   isotopes with the same atomic number. Alternatively, the
            #   FineIsotopePatternGenerator can be used that calculates hyperfine isotopic
            #   distributions. 
            #    -----
            #   This class only describes the container that holds the isotopic
            #   distribution, calculations are done using classes derived from
            #   IsotopePatternGenerator.

        IsotopeDistribution(IsotopeDistribution) nogil except + # wrap-ignore

        # overwrites the container which holds the distribution using @p distribution
        void set(libcpp_vector[ Peak1D ]& distribution) nogil except + # wrap-doc:Overwrites the container which holds the distribution using 'distribution'

        void insert(double mass, float intensity) nogil except +

        libcpp_vector[Peak1D].iterator begin() nogil except +  # wrap-iter-begin:__iter__(Peak1D)
        libcpp_vector[Peak1D].iterator end()   nogil except +  # wrap-iter-end:__iter__(Peak1D)

        # returns the container which holds the distribution
        libcpp_vector[ Peak1D ]& getContainer() nogil except + # wrap-doc:Returns the container which holds the distribution

        # returns the maximal weight isotope which is stored in the distribution
        Size getMax() nogil except + # wrap-doc:Returns the maximal weight isotope which is stored in the distribution

        # returns the minimal weight isotope which is stored in the distribution
        Size getMin() nogil except + # wrap-doc:Returns the minimal weight isotope which is stored in the distribution

        # returns the most abundant isotope which is stored in the distribution
        Peak1D getMostAbundant() nogil except + # wrap-doc:Returns the most abundant isotope which is stored in the distribution

        # returns the size of the distribution which is the number of isotopes in the distribution
        Size size() nogil except + # wrap-doc:Returns the size of the distribution which is the number of isotopes in the distribution

        # clears the distribution and resets max isotope to 0
        void clear() nogil except + # wrap-doc:Clears the distribution and resets max isotope to 0

        # renormalizes the sum of the probabilities of the isotopes to 1
        void renormalize() nogil except + # wrap-doc:Renormalizes the sum of the probabilities of the isotopes to 1

        # Trims the right side of the isotope distribution to isotopes with a significant contribution.
        void trimRight(double cutoff) nogil except + # wrap-doc:Trims the right side of the isotope distribution to isotopes with a significant contribution

        # Trims the left side of the isotope distribution to isotopes with a significant contribution.
        void trimLeft(double cutoff) nogil except + # wrap-doc:Trims the left side of the isotope distribution to isotopes with a significant contribution
        
        void merge(double, double) nogil except + # wrap-doc:Merges distributions of arbitrary data points with constant defined resolution

        void resize(UInt size) nogil except + # wrap-doc:Resizes distribution container
        void trimIntensities(double cutoff) nogil except + # wrap-doc:Remove intensities below the cutoff
        void sortByIntensity() nogil except + # wrap-doc:Sort isotope distribution by intensity
        void sortByMass() nogil except + # wrap-doc:Sort isotope distribution by mass
        double averageMass() nogil except + # wrap-doc:Compute average mass of isotope distribution (weighted average of all isotopes)

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>" namespace "OpenMS::IsotopeDistribution":
    
    cdef enum Sorted "OpenMS::IsotopeDistribution::Sorted":
        # wrap-attach:
        #    IsotopeDistribution
        INTENSITY
        MASS
        UNDEFINED

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>" namespace "OpenMS":

    cdef cppclass FineIsotopePatternGenerator:
        #
        # wrap-doc:
        #   Isotope pattern generator for fine isotope distributions.
        #   Generates isotopes until a stop condition (threshold) is reached,
        #   the lower the threshold the more isotopes are generated. The
        #   parameter use_total_prob defines whether the stop condition is
        #   interpreted as the total probability that the distribution should
        #   cover (default) or as a threshold for individual peaks. Finally,
        #   the absolute parameter specifies for individual peak thresholding
        #   if the threshold is absolute or relative.

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
        bool getRoundMasses() nogil except + # wrap-doc:Returns the current value of the flag to round masses to integer values (true) or return accurate masses (false)

        # sets the round_masses_ flag to round masses to integer values (true) or return accurate masses (false)
        void setRoundMasses(bool round_masses_) nogil except + # wrap-doc:Sets the round_masses_ flag to round masses to integer values (true) or return accurate masses (false)

        # returns the currently set maximum isotope
        Size getMaxIsotope() nogil except + # wrap-doc:Returns the currently set maximum isotope
        
        #  sets the maximal isotope with @p max_isotope
        void setMaxIsotope(Size max_isotope) nogil except + # wrap-doc:Sets the maximal isotope with 'max_isotope'

        # @brief Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported
        #   "Determination of Monoisotopic Masses and Ion Populations for Large Biomolecules from Resolved Isotopic Distributions"
        IsotopeDistribution estimateFromPeptideWeight(double average_weight) nogil except + # wrap-doc:Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported

        # Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs
        IsotopeDistribution estimateFromPeptideWeightAndS(double average_weight, UInt S) nogil except + # wrap-doc:Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs

        # Estimate Nucleotide Isotopedistribution from weight
        IsotopeDistribution estimateFromRNAWeight(double average_weight) nogil except + # wrap-doc:Estimate Nucleotide Isotopedistribution from weight

        # Estimate Nucleotide Isotopedistribution from weight
        IsotopeDistribution estimateFromDNAWeight(double average_weight) nogil except + # wrap-doc:Estimate Nucleotide Isotopedistribution from weight

        IsotopeDistribution estimateFromWeightAndComp(double average_weight,
                                                      double C, double H,
                                                      double N, double O,
                                                      double S, double P) nogil except +

        # Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition
        IsotopeDistribution estimateFromWeightAndCompAndS(double average_weight,
                                                          UInt S, double C,
                                                          double H, double N,
                                                          double O, double P) nogil except + # wrap-doc:Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition

        # Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromPeptideWeight(double average_weight_precursor,
                                                                 double average_weight_fragment,
                                                                 libcpp_set[ unsigned int ]& precursor_isotopes) nogil except + # wrap-doc:Estimate peptide fragment IsotopeDistribution from the precursor's average weight, fragment's average weight, and a set of isolated precursor isotopes

        # Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        # number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
        # and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor,
                                                                     UInt S_precursor,
                                                                     double average_weight_fragment,
                                                                     UInt S_fragment,
                                                                     libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +
            # wrap-doc:
                #   Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
                #   number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
                #   and a set of isolated precursor isotopes.

        # Estimate RNA fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromRNAWeight(double average_weight_precursor,
                                                             double average_weight_fragment,
                                                             libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +
            # wrap-doc:
                #   Estimate RNA fragment IsotopeDistribution from the precursor's average weight,
                #   fragment's average weight, and a set of isolated precursor isotopes.

        # Estimate DNA fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        IsotopeDistribution estimateForFragmentFromDNAWeight(double average_weight_precursor,
                                                             double average_weight_fragment,
                                                             libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +
            # wrap-doc:
                #   Estimate DNA fragment IsotopeDistribution from the precursor's average weight,
                #   fragment's average weight, and a set of isolated precursor isotopes.

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
            # wrap-doc:
                #   Estimate fragment IsotopeDistribution from the precursor's average weight,
                #   fragment's average weight, a set of isolated precursor isotopes, and average composition

        # Calculate isotopic distribution for a fragment molecule
        IsotopeDistribution calcFragmentIsotopeDist(IsotopeDistribution& fragment_isotope_dist,
                                                    IsotopeDistribution& comp_fragment_isotope_dist,
                                                    libcpp_set[ unsigned int ]& precursor_isotopes,
                                                    double fragment_mono_mass) nogil except + # wrap-doc:Calculate isotopic distribution for a fragment molecule

