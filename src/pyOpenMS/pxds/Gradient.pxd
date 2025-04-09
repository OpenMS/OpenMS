from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/Gradient.h>" namespace "OpenMS":

    cdef cppclass Gradient:

        Gradient() except + nogil  # wrap-doc:Representation of a HPLC gradient
        Gradient(Gradient &) except + nogil 

        #  @brief Adds an eluent at the end of the eluent array
        #
        #  @exception Exception::InvalidValue is thrown if the same eluent name is used twice.
        void addEluent(String eluent) except + nogil  # wrap-doc:Adds an eluent at the end of the eluent array
        
        void clearEluents() except + nogil  # wrap-doc:Removes all eluents
        
        libcpp_vector[String] getEluents() except + nogil  # wrap-doc:Returns a reference to the list of eluents

        #  @brief Adds a timepoint at the end of the timepoint array
        #  @exception Exception::OutOfRange is thrown if the new timpoint is before the last timepoint.
        void addTimepoint(Int timepoint) except + nogil  # wrap-doc:Adds a timepoint at the end of the timepoint array
        
        void clearTimepoints() except + nogil  # wrap-doc:Removes all timepoints
        
        libcpp_vector[Int] getTimepoints() except + nogil  # wrap-doc:Returns a reference to the list of timepoints

        #  @brief sets the percentage of eluent @p eluent at timepoint @p timepoint
        #  @exception Exception::InvalidValue is thrown if the eluent, timepoint or percentage is invalid.
        void setPercentage(String eluent, Int timepoint, UInt percentage) except + nogil  # wrap-doc:Sets the percentage of 'eluent' at 'timepoint'
        UInt getPercentage(String eluent, Int timepoint) except + nogil  # wrap-doc:Returns a const reference to the percentages

        # TODO 
        # libcpp_vector[ libcpp_vector[unsigned int] ] getPercentages() except + nogil 

    
        void clearPercentages() except + nogil  # wrap-doc:Sets all percentage values to 0
    
        bool isValid() except + nogil  # wrap-doc:Checks if the percentages of all timepoints add up to 100%
