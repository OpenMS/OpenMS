from libcpp cimport bool
from Types cimport *
from String cimport *
from Element cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/EmpiricalFormula.h>" namespace "OpenMS":

    cdef cppclass EmpiricalFormula:

        EmpiricalFormula() nogil except +
        EmpiricalFormula(EmpiricalFormula) nogil except + # wrap-ignore

        # constructor with element pointer and number
        ## EmpiricalFormula(SignedSize number, Element * element, SignedSize charge = 0) nogil except +

        # returns the mono isotopic weight of the formula (includes proton charges)
        DoubleReal getMonoWeight() nogil except +

        # returns the average weight of the formula (includes proton charges)
        DoubleReal getAverageWeight() nogil except +

        # @brief returns the isotope distribution of the formula
        #   *	The details of the calculation of the isotope distribution
        #   * are described in the doc to the IsotopeDistribution class.
        #   *	@param max_depth: this parameter gives the max isotope which is considered, if 0 all are reported
        IsotopeDistribution getIsotopeDistribution(UInt max_depth) nogil except +

        ## # returns a pointer to the element with name or symbol or 0 if no such element is fount
        ## Element * getElement(String & name) nogil except +

        ## # returns a pointer to the element with given atomic number or 0 if none if found
        ## Element * getElement(UInt atomic_number) nogil except +

        ## # returns a pointer to the element db which is used with this class
        ## ElementDB * getElementDB() nogil except +

        # returns the number of atoms with the given atomic_number
        Size getNumberOf(UInt atomic_number) nogil except +

        # returns the number of atoms with the given name
        Size getNumberOf(String & name) nogil except +

        # returns the number of atoms
        Size getNumberOf(Element * element) nogil except +

        # returns the atoms total
        Size getNumberOfAtoms() nogil except +

        # returns the charge
        SignedSize getCharge() nogil except +

        # sets the charge
        void setCharge(SignedSize charge) nogil except +

        # returns the formula as a string (charges are not included)
        String getString() nogil except +

        # TargetedExperiment iadd(TargetedExperiment)   nogil except + 
        # adds the elements of the given formula
        EmpiricalFormula iadd(EmpiricalFormula) nogil except + # wrap-as:operator+=

"""
        /** adds the elements from the given formula, which is given as a OpenMS String

                @throw throws ParseError if the formula cannot be parsed
        */
        EmpiricalFormula & operator+=(String & rhs) nogil except +

        # multiplies the elements and charge with a factor
        EmpiricalFormula operator*(SignedSize & times) nogil except +

        # adds the elements of the given formula and returns a new formula
        EmpiricalFormula operator+(EmpiricalFormula & rhs) nogil except +

        /** adds the elements of the given formula (given as a String) and returns a new formula

                @throw throws ParseError if the formula cannot be parsed
        */
        EmpiricalFormula operator+(String & rhs) nogil except +

        # subtracts the elements of a formula
        EmpiricalFormula & operator-=(EmpiricalFormula & rhs) nogil except +

        /** subtracts the elements of a formula given as string

                @throw throws ParseError if the formula cannot be parsed
        */
        EmpiricalFormula & operator-=(String & rhs) nogil except +

        # subtracts the elements of a formula an returns a new formula
        EmpiricalFormula operator-(EmpiricalFormula & rhs) nogil except +

        /** subtracts the elements of a formula given as a String and returns a new formula

                @throw throws ParseError if the formula cannot be parsed
        */
        EmpiricalFormula operator-(String & rhs) nogil except +
        #@}

        /**@name Predicates
        */
        #@{
        # returns true if the formula does not contain a element
        bool isEmpty() nogil except +

        # returns true if charge != 0
        bool isCharged() nogil except +

        # returns true if the formula contains the element
        bool hasElement(Element * element) nogil except +

        # returns true if the formula contains the element, given with its name or symbol
        bool hasElement(String & name) nogil except +

        # returns true if the formula contains the element with the given atomic number
        bool hasElement(UInt atomic_number) nogil except +

        # returns true if the formulas contain equal elements in equal quantities
        bool operator==(EmpiricalFormula & rhs) nogil except +

        /** returns true if the formulas contain equal elements in equal quantities

                @throw throws ParseError if the formula cannot be parsed
        */
        bool operator==(String & rhs) nogil except +

        # returns true if the formulas differ in elements composition
        bool operator!=(EmpiricalFormula & rhs) nogil except +

        /** returns true if the formulas differ in elements composition

                @throw throws ParseError if the formula cannot be parsed
        */
        bool operator!=(String & rhs) nogil except +
        #@}

        # writes the formula to a stream
        friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, EmpiricalFormula & formula) nogil except +

        /** @name Iterators
        */
        #@{
        inline ConstIterator begin() { return formula_.begin() nogil except + }

        inline ConstIterator end() { return formula_.end() nogil except + }
        #@}

    protected:

        # remove elements with count 0
        void removeZeroedElements_() nogil except +

        Map<Element *, SignedSize] formula_ nogil except +

        SignedSize charge_ nogil except +

        void readElementsFromFile_(String & file_name) nogil except +

        SignedSize parseFormula_(Map<Element *, SignedSize] & ef, String & formula) nogil except +

        ElementDB * element_db_ nogil except +
      } nogil except +

"""
