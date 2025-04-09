from libcpp cimport bool
from Types cimport *
from String cimport *
from Peak1D cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>" namespace "OpenMS":

    cdef cppclass IsotopeDistribution:
        # wrap-doc:
            #  Isotope distribution class
            #  
            #  A container that holds an isotope distribution. It consists of mass values
            #  and their correspondent probabilities (stored in the intensity slot)
            #  
            #  Isotope distributions can be calculated using either the
            #  CoarseIsotopePatternGenerator for quantized atomic masses which group
            #  isotopes with the same atomic number. Alternatively, the
            #  FineIsotopePatternGenerator can be used that calculates hyperfine isotopic
            #  distributions 
            #   
            #  This class only describes the container that holds the isotopic
            #  distribution, calculations are done using classes derived from
            #  IsotopePatternGenerator        

        IsotopeDistribution() except + nogil 
        IsotopeDistribution(IsotopeDistribution &) except + nogil 

        void set(libcpp_vector[ Peak1D ]& distribution) except + nogil  # wrap-doc:Overwrites the container which holds the distribution using 'distribution'

        void insert(double mass, float intensity) except + nogil 

        libcpp_vector[Peak1D].iterator begin() except + nogil   # wrap-iter-begin:__iter__(Peak1D)
        libcpp_vector[Peak1D].iterator end()   except + nogil   # wrap-iter-end:__iter__(Peak1D)

        libcpp_vector[ Peak1D ]& getContainer() except + nogil  # wrap-doc:Returns the container which holds the distribution

        double getMax() except + nogil  # wrap-doc:Returns the maximal weight isotope which is stored in the distribution

        double getMin() except + nogil  # wrap-doc:Returns the minimal weight isotope which is stored in the distribution

        Peak1D getMostAbundant() except + nogil  # wrap-doc:Returns the most abundant isotope which is stored in the distribution

        Size size() except + nogil  # wrap-doc:Returns the size of the distribution which is the number of isotopes in the distribution

        void clear() except + nogil  # wrap-doc:Clears the distribution and resets max isotope to 0

        void renormalize() except + nogil  # wrap-doc:Renormalizes the sum of the probabilities of the isotopes to 1

        void trimRight(double cutoff) except + nogil  # wrap-doc:Trims the right side of the isotope distribution to isotopes with a significant contribution

        void trimLeft(double cutoff) except + nogil  # wrap-doc:Trims the left side of the isotope distribution to isotopes with a significant contribution
        
        void merge(double, double) except + nogil  # wrap-doc:Merges distributions of arbitrary data points with constant defined resolution

        void resize(UInt size) except + nogil  # wrap-doc:Resizes distribution container
        void trimIntensities(double cutoff) except + nogil  # wrap-doc:Remove intensities below the cutoff
        void sortByIntensity() except + nogil  # wrap-doc:Sort isotope distribution by intensity
        void sortByMass() except + nogil  # wrap-doc:Sort isotope distribution by mass
        double averageMass() except + nogil  # wrap-doc:Compute average mass of isotope distribution (weighted average of all isotopes)

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>" namespace "OpenMS::IsotopeDistribution":
    
    cdef enum Sorted "OpenMS::IsotopeDistribution::Sorted":
        # wrap-attach:
        #   IsotopeDistribution
        INTENSITY
        MASS
        UNDEFINED

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>" namespace "OpenMS":

    cdef cppclass FineIsotopePatternGenerator:
        #
        # wrap-doc:
        #  Isotope pattern generator for fine isotope distributions.
        #  Generates isotopes until a stop condition (threshold) is reached,
        #  the lower the threshold the more isotopes are generated. The
        #  parameter use_total_prob defines whether the stop condition is
        #  interpreted as the total probability that the distribution should
        #  cover (default) or as a threshold for individual peaks. Finally,
        #  the absolute parameter specifies for individual peak thresholding
        #  if the threshold is absolute or relative.

        FineIsotopePatternGenerator() except + nogil  
        FineIsotopePatternGenerator(double threshold) except + nogil 
        FineIsotopePatternGenerator(double threshold, bool use_total_prob) except + nogil 
        FineIsotopePatternGenerator(double threshold, bool use_total_prob, bool absolute) except + nogil 

        void setThreshold(double threshold) except + nogil 
        double getThreshold() except + nogil 

        void setAbsolute(bool absolute) except + nogil 
        bool getAbsolute() except + nogil 

        void setTotalProbability(bool total) except + nogil 
        bool getTotalProbability() except + nogil 

        IsotopeDistribution run(EmpiricalFormula) except + nogil 

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>" namespace "OpenMS":

    cdef cppclass CoarseIsotopePatternGenerator:

        CoarseIsotopePatternGenerator() except + nogil  
        CoarseIsotopePatternGenerator(Size max_isotope) except + nogil 
        CoarseIsotopePatternGenerator(Size max_isotope, bool round_masses) except + nogil 

        IsotopeDistribution run(EmpiricalFormula) except + nogil 

        bool getRoundMasses() except + nogil  # wrap-doc:Returns the current value of the flag to round masses to integer values (true) or return accurate masses (false)

        void setRoundMasses(bool round_masses_) except + nogil  # wrap-doc:Sets the round_masses_ flag to round masses to integer values (true) or return accurate masses (false)

        Size getMaxIsotope() except + nogil  # wrap-doc:Returns the currently set maximum isotope
        
        void setMaxIsotope(Size max_isotope) except + nogil  # wrap-doc:Sets the maximal isotope with 'max_isotope'

        #  "Determination of Monoisotopic Masses and Ion Populations for Large Biomolecules from Resolved Isotopic Distributions"
        IsotopeDistribution estimateFromPeptideWeight(double average_weight) except + nogil  # wrap-doc:Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported

        IsotopeDistribution estimateFromPeptideWeightAndS(double average_weight, UInt S) except + nogil  # wrap-doc:Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs

        IsotopeDistribution estimateFromRNAWeight(double average_weight) except + nogil  # wrap-doc:Estimate Nucleotide Isotopedistribution from weight

        IsotopeDistribution estimateFromDNAWeight(double average_weight) except + nogil  # wrap-doc:Estimate Nucleotide Isotopedistribution from weight

        IsotopeDistribution estimateFromWeightAndComp(double average_weight,
                                                      double C, double H,
                                                      double N, double O,
                                                      double S, double P) except + nogil 

        # Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition
        IsotopeDistribution estimateFromWeightAndCompAndS(double average_weight,
                                                          UInt S, double C,
                                                          double H, double N,
                                                          double O, double P) except + nogil  # wrap-doc:Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition

        IsotopeDistribution estimateForFragmentFromPeptideWeight(double average_weight_precursor,
                                                                 double average_weight_fragment,
                                                                 libcpp_set[ unsigned int ]& precursor_isotopes) except + nogil  # wrap-doc:Estimate peptide fragment IsotopeDistribution from the precursor's average weight, fragment's average weight, and a set of isolated precursor isotopes

        IsotopeDistribution estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor,
                                                                     UInt S_precursor,
                                                                     double average_weight_fragment,
                                                                     UInt S_fragment,
                                                                     libcpp_set[ unsigned int ]& precursor_isotopes) except + nogil 
            # wrap-doc:
                #  Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
                #  number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
                #  and a set of isolated precursor isotopes.

        IsotopeDistribution approximateFromPeptideWeight(double mass, UInt num_peaks, UInt charge) except + nogil 

            # wrap-doc:
                #  Roughly approximate peptide IsotopeDistribution from monoisotopic weight using Poisson distribution.
                #  m/z values approximated by adding one neutron mass (divided by charge) for every peak, starting at 
		#  the given monoisotopic weight. Foundation from: Bellew et al, https://dx.doi.org/10.1093/bioinformatics/btl276
		#  This method is around 50 times faster than estimateFromPeptideWeight, but only an approximation.
		#  The following are the intensities of the first 6 peaks generated for a monoisotopic mass of 1000:
		#  estimateFromPeptideWeight:    0.571133000;0.306181000;0.095811100;0.022036900;0.004092170;0.000644568
		#  approximateFromPeptideWeight: 0.573753000;0.318752000;0.088542200;0.016396700;0.002277320;0.000253036
		#  KL divergences of the first 20 intensities of estimateFromPeptideWeight and this approximation range from 4.97E-5 for a
		#  monoisotopic mass of 20 to 0.0144 for a mass of 2500. For comparison, when comparing an observed pattern with a 
		#  theoretical ground truth, the observed pattern is said to be an isotopic pattern if the KL between the two is below 0.05
		#  for 2 peaks and below 0.6 for >=6 peaks by Guo Ci Teo et al.

        libcpp_vector[ double ] approximateIntensities(double mass, UInt num_peaks) except + nogil 

            # wrap-doc:
                #  Roughly approximate peptidic isotope pattern intensities from monoisotopic weight using Poisson distribution.
                #  Foundation from: Bellew et al, https://dx.doi.org/10.1093/bioinformatics/btl276
		#  This method is around 100 times faster than estimateFromPeptideWeight, but only an approximation, see approximateFromPeptideWeight.

        IsotopeDistribution estimateForFragmentFromRNAWeight(double average_weight_precursor,
                                                             double average_weight_fragment,
                                                             libcpp_set[ unsigned int ]& precursor_isotopes) except + nogil 
            # wrap-doc:
                #  Estimate RNA fragment IsotopeDistribution from the precursor's average weight,
                #  fragment's average weight, and a set of isolated precursor isotopes

        IsotopeDistribution estimateForFragmentFromDNAWeight(double average_weight_precursor,
                                                             double average_weight_fragment,
                                                             libcpp_set[ unsigned int ]& precursor_isotopes) except + nogil 
            # wrap-doc:
                #  Estimate DNA fragment IsotopeDistribution from the precursor's average weight,
                #  fragment's average weight, and a set of isolated precursor isotopes.

        IsotopeDistribution estimateForFragmentFromWeightAndComp(double average_weight_precursor,
                                                                 double average_weight_fragment,
                                                                 libcpp_set[ unsigned int ]& precursor_isotopes,
                                                                 double C,
                                                                 double H,
                                                                 double N,
                                                                 double O,
                                                                 double S,
                                                                 double P) except + nogil 
            # wrap-doc:
                #  Estimate fragment IsotopeDistribution from the precursor's average weight,
                #  fragment's average weight, a set of isolated precursor isotopes, and average composition

        IsotopeDistribution calcFragmentIsotopeDist(IsotopeDistribution& fragment_isotope_dist,
                                                    IsotopeDistribution& comp_fragment_isotope_dist,
                                                    libcpp_set[ unsigned int ]& precursor_isotopes,
                                                    double fragment_mono_mass) except + nogil  # wrap-doc:Calculate isotopic distribution for a fragment molecule

