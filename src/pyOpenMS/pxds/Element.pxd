from libcpp cimport bool
from Types cimport *
from String cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/Element.h>" namespace "OpenMS":

    cdef cppclass Element:

        Element() except + nogil 
        Element(Element &) except + nogil 

        # detailed constructor
        Element(String name,
                String symbol,
                UInt atomic_number,
                double average_weight,
                double mono_weight,
                IsotopeDistribution isotopes) except + nogil 

        # sets unique atomic number
        void setAtomicNumber(UInt atomic_number) except + nogil  # wrap-doc:Sets unique atomic number

        # returns the unique atomic number
        UInt getAtomicNumber() except + nogil  # wrap-doc:Returns the unique atomic number

        # sets the average weight of the element
        void setAverageWeight(double weight) except + nogil  # wrap-doc:Sets the average weight of the element

        # returns the average weight of the element
        double getAverageWeight() except + nogil  # wrap-doc:Returns the average weight of the element

        # sets the mono isotopic weight of the element
        void setMonoWeight(double weight) except + nogil  # wrap-doc:Sets the mono isotopic weight of the element

        # returns the mono isotopic weight of the element
        double getMonoWeight() except + nogil  # wrap-doc:Returns the mono isotopic weight of the element

        # sets the isotope distribution of the element
        void setIsotopeDistribution(IsotopeDistribution isotopes) except + nogil  # wrap-doc:Sets the isotope distribution of the element

        # returns the isotope distribution of the element
        IsotopeDistribution getIsotopeDistribution() except + nogil  # wrap-doc:Returns the isotope distribution of the element

        # set the name of the element
        void setName(String name) except + nogil  # wrap-doc:Sets the name of the element

        # returns the name of the element
        String getName() except + nogil  # wrap-doc:Returns the name of the element

        # sets symbol of the element
        void setSymbol(String symbol) except + nogil  # wrap-doc:Sets symbol of the element

        # returns symbol of the element
        String getSymbol() except + nogil  # wrap-doc:Returns symbol of the element

