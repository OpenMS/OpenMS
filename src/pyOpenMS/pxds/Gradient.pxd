from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/Gradient.h>" namespace "OpenMS":

    cdef cppclass Gradient:

        Gradient() nogil except +
        Gradient(Gradient) nogil except + # wrap-ignore

        #   @brief Adds an eluent at the end of the eluent array
        #
        #   @exception Exception::InvalidValue is thrown if the same eluent name is used twice.
        void addEluent(String eluent) nogil except +
        # removes all eluents
        void clearEluents() nogil except +
        # returns a reference to the list of eluents
        libcpp_vector[String] getEluents() nogil except +

        #   @brief Adds a timepoint at the end of the timepoint array
        #
        #   @exception Exception::OutOfRange is thrown if the new timpoint is before the last timepoint.
        void addTimepoint(Int timepoint) nogil except +
        # removes all timepoints
        void clearTimepoints() nogil except +
        # returns a reference to the list of timepoints
        libcpp_vector[Int] getTimepoints() nogil except +

        #   @brief sets the percentage of eluent @p eluent at timepoint @p timepoint
        #
        #   @exception Exception::InvalidValue is thrown if the eluent, timepoint or percentage is invalid.
        void setPercentage(String eluent, Int timepoint, UInt percentage) nogil except +
        UInt getPercentage(String eluent, Int timepoint) nogil except +

        # TODO 
        # libcpp_vector[ libcpp_vector[unsigned int] ] getPercentages() nogil except +

        # sets all percentage values to 0
        void clearPercentages() nogil except +
        # checks if the percentages of all timepoints add up to 100%
        bool isValid() nogil except +

