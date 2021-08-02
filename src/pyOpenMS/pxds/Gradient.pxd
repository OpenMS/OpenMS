from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/Gradient.h>" namespace "OpenMS":

    cdef cppclass Gradient:

        Gradient() nogil except + # wrap-doc:Representation of a HPLC gradient
        Gradient(Gradient &) nogil except +

        #   @brief Adds an eluent at the end of the eluent array
        #
        #   @exception Exception::InvalidValue is thrown if the same eluent name is used twice.
        void addEluent(String eluent) nogil except + # wrap-doc:Adds an eluent at the end of the eluent array
        # removes all eluents
        void clearEluents() nogil except + # wrap-doc:Removes all eluents
        # returns a reference to the list of eluents
        libcpp_vector[String] getEluents() nogil except + # wrap-doc:Returns a reference to the list of eluents

        #   @brief Adds a timepoint at the end of the timepoint array
        #   @exception Exception::OutOfRange is thrown if the new timpoint is before the last timepoint.
        void addTimepoint(Int timepoint) nogil except + # wrap-doc:Adds a timepoint at the end of the timepoint array
        # removes all timepoints
        void clearTimepoints() nogil except + # wrap-doc:Removes all timepoints
        # returns a reference to the list of timepoints
        libcpp_vector[Int] getTimepoints() nogil except + # wrap-doc:Returns a reference to the list of timepoints

        #   @brief sets the percentage of eluent @p eluent at timepoint @p timepoint
        #   @exception Exception::InvalidValue is thrown if the eluent, timepoint or percentage is invalid.
        void setPercentage(String eluent, Int timepoint, UInt percentage) nogil except + # wrap-doc:Sets the percentage of 'eluent' at 'timepoint'
        UInt getPercentage(String eluent, Int timepoint) nogil except + # wrap-doc:Returns a const reference to the percentages

        # TODO 
        # libcpp_vector[ libcpp_vector[unsigned int] ] getPercentages() nogil except +

        # sets all percentage values to 0
        void clearPercentages() nogil except + # wrap-doc:Sets all percentage values to 0
        # checks if the percentages of all timepoints add up to 100%
        bool isValid() nogil except + # wrap-doc:Checks if the percentages of all timepoints add up to 100%
