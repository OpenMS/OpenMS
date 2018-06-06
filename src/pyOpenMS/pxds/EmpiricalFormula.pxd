from libcpp cimport bool
from Types cimport *
from String cimport *
from Element cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/EmpiricalFormula.h>" namespace "OpenMS":

    cdef cppclass EmpiricalFormula:

        EmpiricalFormula() nogil except +
        EmpiricalFormula(EmpiricalFormula) nogil except + # wrap-ignore

        # constructor from string
        EmpiricalFormula(String) nogil except +

        # constructor with element pointer and number
        EmpiricalFormula(SignedSize number, Element * element, SignedSize charge) nogil except +

        # returns the mono isotopic weight of the formula (includes proton charges)
        double getMonoWeight() nogil except +

        # returns the average weight of the formula (includes proton charges)
        double getAverageWeight() nogil except +

        # Fills this EmpiricalFormula with an approximate elemental composition for a given average weight and approximate elemental stoichiometry
        bool estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P) nogil except +

        # Fills this EmpiricalFormula with an approximate elemental composition for a given average weight,
        # exact number of sulfurs, and approximate elemental stoichiometry
        bool estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P) nogil except +


        #  Computes the isotope distribution of an empirical formula using the CoarseIsotopePatternGenerator method
        IsotopeDistribution getIsotopeDistribution(CoarseIsotopePatternGenerator) nogil except +

        # @brief returns the fragment isotope distribution of this conditioned
        # on a precursor formula and a list of isolated precursor isotopes.
        # @param precursor: the empirical formula of the precursor
        # @param precursor_isotopes: the set of precursor isotopes that were isolated
        IsotopeDistribution getConditionalFragmentIsotopeDist(EmpiricalFormula& precursor, libcpp_set[ unsigned int ]& precursor_isotopes) nogil except +

        # returns the number of atoms
        Size getNumberOf(Element * element) nogil except +

        # returns the atoms total
        Size getNumberOfAtoms() nogil except +

        # returns the charge
        SignedSize getCharge() nogil except +

        # sets the charge
        void setCharge(SignedSize charge) nogil except +

        # returns the formula as a string (charges are not included)
        String toString() nogil except +

        # returns true if the formula does not contain a element
        bool isEmpty() nogil except +

        # returns true if charge != 0
        bool isCharged() nogil except +

        # returns true if the formula contains the element
        bool hasElement(Element * element) nogil except +

        # returns true if all elements from @p ef are LESS abundant (negative allowed) than the corresponding elements of this EmpiricalFormula
        bool contains(EmpiricalFormula ef) nogil except +

        # returns true if the formulas contain equal elements in equal quantities
        bool operator==(EmpiricalFormula) nogil except +

        # returns true if the formulas differ in elements composition
        bool operator!=(EmpiricalFormula) nogil except +

        EmpiricalFormula operator+(EmpiricalFormula) nogil except +
        # EmpiricalFormula operator-(EmpiricalFormula) nogil except +
        # EmpiricalFormula operator*(EmpiricalFormula) nogil except +

        EmpiricalFormula iadd(EmpiricalFormula)   nogil except + # wrap-as:operator+=
        # EmpiricalFormula iminus(EmpiricalFormula)   nogil except + # wrap-as:operator-=


