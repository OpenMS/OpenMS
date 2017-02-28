from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/IsotopeDistribution.h>" namespace "OpenMS":

    cdef cppclass IsotopeDistribution:

        IsotopeDistribution() nogil except +
        IsotopeDistribution(IsotopeDistribution) nogil except + # wrap-ignore

        #  sets the maximal isotope with @p max_isotope
        void setMaxIsotope(Size max_isotope) nogil except +

        # returns the currently set maximum isotope
        Size getMaxIsotope() nogil except +

        # overwrites the container which holds the distribution using @p distribution
        void set(libcpp_vector[ libcpp_pair[ size_t, double] ] & distribution) nogil except +

        # returns the container which holds the distribution
        libcpp_vector[ libcpp_pair[ size_t, double] ] getContainer() nogil except +

        # returns the maximal weight isotope which is stored in the distribution
        Size getMax() nogil except +

        # returns the minimal weight isotope which is stored in the distribution
        Size getMin() nogil except +

        # returns the size of the distribtion which is the number of isotopes in the distribution
        Size size() nogil except +

        # clears the distribution and resets max isotope to 0
        void clear() nogil except +

        # Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported
        #   Implementation using the averagine model proposed by Senko et al. in
        #   "Determination of Monoisotopic Masses and Ion Populations for Large Biomolecules from Resolved Isotopic Distributions"
        void estimateFromPeptideWeight(double average_weight) nogil except +

        # Estimate Nucleotide Isotopedistribution from weight
        void estimateFromRNAWeight(double average_weight) nogil except +

        # Estimate Nucleotide Isotopedistribution from weight

        void estimateFromDNAWeight(double average_weight) nogil except +

        # Estimate Isotopedistribution from weight, average composition, and number of isotopes that should be reported

        void estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P) nogil except +

        # Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        void estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # Estimate RNA fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        void estimateForFragmentFromRNAWeight(double average_weight_precursor, double average_weight_fragment, libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # Estimate DNA fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, and a set of isolated precursor isotopes.
        void estimateForFragmentFromDNAWeight(double average_weight_precursor, double average_weight_fragment, libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # Estimate fragment IsotopeDistribution from the precursor's average weight,
        # fragment's average weight, a set of isolated precursor isotopes, and average composition
        void estimateForFragmentFromWeightAndComp(double average_weight_precursor, double average_weight_fragment, libcpp_set[ unsigned int ]& precursor_isotopes, double C, double H, double N, double O, double S, double P) nogil except +

        # Calculate isotopic distribution for a fragment molecule

        void calcFragmentIsotopeDist(IsotopeDistribution& fragment_isotope_dist, IsotopeDistribution& comp_fragment_isotope_dist, libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # renormalizes the sum of the probabilities of the isotopes to 1
        void renormalize() nogil except +

        # Trims the right side of the isotope distribution to isotopes with a significant contribution.
        void trimRight(double cutoff) nogil except +

        # Trims the left side of the isotope distribution to isotopes with a significant contribution.
        void trimLeft(double cutoff) nogil except +

