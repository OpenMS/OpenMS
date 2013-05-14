from libcpp cimport bool
from Types cimport *
from String cimport *
from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/Element.h>" namespace "OpenMS":

    cdef cppclass Element:

        Element() nogil except +
        Element(Element) nogil except + # wrap-ignore

        # detailed constructor
        Element(String name,
                String symbol,
                UInt atomic_number,
                DoubleReal average_weight,
                DoubleReal mono_weight,
                IsotopeDistribution isotopes) nogil except +

        # sets unique atomic number
        void setAtomicNumber(UInt atomic_number) nogil except +

        # returns the unique atomic number
        UInt getAtomicNumber() nogil except +

        # sets the average weight of the element
        void setAverageWeight(DoubleReal weight) nogil except +

        # returns the average weight of the element
        DoubleReal getAverageWeight() nogil except +

        # sets the mono isotopic weight of the element
        void setMonoWeight(DoubleReal weight) nogil except +

        # returns the mono isotopic weight of the element
        DoubleReal getMonoWeight() nogil except +

        # sets the isotope distribution of the element
        void setIsotopeDistribution(IsotopeDistribution isotopes) nogil except +

        # returns the isotope distribution of the element
        IsotopeDistribution getIsotopeDistribution() nogil except +

        # set the name of the element
        void setName(String name) nogil except +

        # returns the name of the element
        String getName() nogil except +

        # sets symbol of the element
        void setSymbol(String symbol) nogil except +

        # returns symbol of the element
        String getSymbol() nogil except +

