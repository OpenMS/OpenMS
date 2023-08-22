from Types cimport *
from String cimport *
from Element cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/EmpiricalFormula.h>" namespace "OpenMS":

    cdef cppclass EmpiricalFormula:

        EmpiricalFormula() except + nogil  # wrap-doc:Representation of an empirical formula
        EmpiricalFormula(EmpiricalFormula &) except + nogil 

        EmpiricalFormula(String) except + nogil  # wrap-doc:EmpiricalFormula Constructor from string

        EmpiricalFormula(SignedSize number, Element * element, SignedSize charge) except + nogil  # wrap-doc:EmpiricalFormula Constructor with element pointer and number

        double getMonoWeight() except + nogil  # wrap-doc:Returns the mono isotopic weight of the formula (includes proton charges)

        double getAverageWeight() except + nogil  # wrap-doc:Returns the average weight of the formula (includes proton charges)

        bool estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P) except + nogil  # wrap-doc:Fills this EmpiricalFormula with an approximate elemental composition for a given average weight and approximate elemental stoichiometry

        # Fills this EmpiricalFormula with an approximate elemental composition for a given average weight,
        # exact number of sulfurs, and approximate elemental stoichiometry
        bool estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P) except + nogil  # wrap-doc:Fills this EmpiricalFormula with an approximate elemental composition for a given average weight, exact number of sulfurs, and approximate elemental stoichiometry

        IsotopeDistribution getIsotopeDistribution(CoarseIsotopePatternGenerator) except + nogil  # wrap-doc:Computes the isotope distribution of an empirical formula using the CoarseIsotopePatternGenerator or the FineIsotopePatternGenerator method
        IsotopeDistribution getIsotopeDistribution(FineIsotopePatternGenerator) except + nogil 

        # returns the fragment isotope distribution of this conditioned
        # on a precursor formula and a list of isolated precursor isotopes.
        # precursor: the empirical formula of the precursor
        # precursor_isotopes: the set of precursor isotopes that were isolated
        # method: the method that will be used for the calculation of the IsotopeDistribution
        IsotopeDistribution getConditionalFragmentIsotopeDist(EmpiricalFormula& precursor, libcpp_set[ unsigned int ]& precursor_isotopes, CoarseIsotopePatternGenerator method) except + nogil  # TODO

        # returns the number of atoms
        # doesn't work!
        ## Size getNumberOf(Element * element) except + nogil 

        Size getNumberOfAtoms() except + nogil  # wrap-doc:Returns the total number of atoms

        SignedSize getCharge() except + nogil  # wrap-doc:Returns the total charge

        void setCharge(SignedSize charge) except + nogil  # wrap-doc:Sets the charge

        String toString() except + nogil  # wrap-doc:Returns the formula as a string (charges are not included)

        libcpp_map[libcpp_string, int] toMap() except + nogil  #wrap-as:getElementalComposition wrap-doc:Get elemental composition as a hash {'Symbol' -> NrAtoms}


        bool isEmpty() except + nogil  # wrap-doc:Returns true if the formula does not contain a element

        bool isCharged() except + nogil  # wrap-doc:Returns true if charge is not equal to zero

        bool hasElement(Element * element) except + nogil  # wrap-doc:Returns true if the formula contains the element

        bool contains(EmpiricalFormula ef) except + nogil  # wrap-doc:Returns true if all elements from `ef` ( empirical formula ) are LESS abundant (negative allowed) than the corresponding elements of this EmpiricalFormula

        bool operator==(EmpiricalFormula) except + nogil  # wrap-doc:Returns true if the formulas contain equal elements in equal quantities

        bool operator!=(EmpiricalFormula) except + nogil  # wrap-doc:Returns true if the formulas differ in elements composition

        EmpiricalFormula operator+(EmpiricalFormula) except + nogil 
        EmpiricalFormula operator-(EmpiricalFormula) except + nogil 
        #EmpiricalFormula operator*(EmpiricalFormula) except + nogil 

        EmpiricalFormula iadd(EmpiricalFormula)   except + nogil  # wrap-as:operator+=
        EmpiricalFormula isub(EmpiricalFormula)   except + nogil  # wrap-as:operator-=

        # autowrap does not support overloaded operators. Only same type
        #EmpiricalFormula imul(unsigned int)   except + nogil  # wrap-as:operator*=

        double calculateTheoreticalIsotopesNumber() except + nogil 
        
