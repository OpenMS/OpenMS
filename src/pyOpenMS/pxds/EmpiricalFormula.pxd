from Types cimport *
from String cimport *
from Element cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/EmpiricalFormula.h>" namespace "OpenMS":

    cdef cppclass EmpiricalFormula:

        EmpiricalFormula() nogil except + # wrap-doc:Representation of an empirical formula
        EmpiricalFormula(EmpiricalFormula &) nogil except +

        EmpiricalFormula(String) nogil except + # wrap-doc:EmpiricalFormula Constructor from string

        EmpiricalFormula(SignedSize number, Element * element, SignedSize charge) nogil except + # wrap-doc:EmpiricalFormula Constructor with element pointer and number

        double getMonoWeight() nogil except + # wrap-doc:Returns the mono isotopic weight of the formula (includes proton charges)

        double getAverageWeight() nogil except + # wrap-doc:Returns the average weight of the formula (includes proton charges)

        bool estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P) nogil except + # wrap-doc:Fills this EmpiricalFormula with an approximate elemental composition for a given average weight and approximate elemental stoichiometry

        # Fills this EmpiricalFormula with an approximate elemental composition for a given average weight,
        # exact number of sulfurs, and approximate elemental stoichiometry
        bool estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P) nogil except + # wrap-doc:Fills this EmpiricalFormula with an approximate elemental composition for a given average weight, exact number of sulfurs, and approximate elemental stoichiometry

        IsotopeDistribution getIsotopeDistribution(CoarseIsotopePatternGenerator) nogil except + # wrap-doc:Computes the isotope distribution of an empirical formula using the CoarseIsotopePatternGenerator or the FineIsotopePatternGenerator method
        IsotopeDistribution getIsotopeDistribution(FineIsotopePatternGenerator) nogil except +

        # returns the fragment isotope distribution of this conditioned
        # on a precursor formula and a list of isolated precursor isotopes.
        # precursor: the empirical formula of the precursor
        # precursor_isotopes: the set of precursor isotopes that were isolated
        # method: the method that will be used for the calculation of the IsotopeDistribution
        IsotopeDistribution getConditionalFragmentIsotopeDist(EmpiricalFormula& precursor, libcpp_set[ unsigned int ]& precursor_isotopes, CoarseIsotopePatternGenerator method) nogil except + # TODO

        # returns the number of atoms
        # doesn't work!
        ## Size getNumberOf(Element * element) nogil except +

        Size getNumberOfAtoms() nogil except + # wrap-doc:Returns the total number of atoms

        SignedSize getCharge() nogil except + # wrap-doc:Returns the total charge

        void setCharge(SignedSize charge) nogil except + # wrap-doc:Sets the charge

        String toString() nogil except + # wrap-doc:Returns the formula as a string (charges are not included)

        libcpp_map[libcpp_string, int] toMap() nogil except + #wrap-as:getElementalComposition wrap-doc:Get elemental composition as a hash {'Symbol' -> NrAtoms}


        bool isEmpty() nogil except + # wrap-doc:Returns true if the formula does not contain a element

        bool isCharged() nogil except + # wrap-doc:Returns true if charge is not equal to zero

        bool hasElement(Element * element) nogil except + # wrap-doc:Returns true if the formula contains the element

        bool contains(EmpiricalFormula ef) nogil except + # wrap-doc:Returns true if all elements from `ef` ( empirical formula ) are LESS abundant (negative allowed) than the corresponding elements of this EmpiricalFormula

        bool operator==(EmpiricalFormula) nogil except + # wrap-doc:Returns true if the formulas contain equal elements in equal quantities

        bool operator!=(EmpiricalFormula) nogil except + # wrap-doc:Returns true if the formulas differ in elements composition

        EmpiricalFormula operator+(EmpiricalFormula) nogil except +
        # EmpiricalFormula operator-(EmpiricalFormula) nogil except +
        # EmpiricalFormula operator*(EmpiricalFormula) nogil except +

        EmpiricalFormula iadd(EmpiricalFormula)   nogil except + # wrap-as:operator+=
        # EmpiricalFormula iminus(EmpiricalFormula)   nogil except + # wrap-as:operator-=

        double calculateTheoreticalIsotopesNumber() nogil except +
        
