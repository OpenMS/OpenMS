from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from MobilityPeak1D cimport *
from String cimport *
from RangeManager cimport *
from DataArrays cimport *
from IMTypes cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/Mobilogram.h>" namespace "OpenMS":

    cdef cppclass Mobilogram (RangeManagerMobInt):
        # wrap-inherits:
        #  RangeManagerMobInt
        #
        # wrap-doc:
        #  The representation of a 1D ion mobilogram.
        #  Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
        #  Iterations yields access to underlying peak objects but is slower
        #  Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
        #  
        #  Usage:
        #
        #  .. code-block:: python
        #  
        #    rt = mobilogram.getRT()
        #    drift_time_unit = mobilogram.getDriftTimeUnit()
        #    mobility, intensities = mobilogram.get_peaks()
        #  

        Mobilogram() except + nogil 
        Mobilogram(Mobilogram &) except + nogil 

        double getRT() except + nogil  # wrap-doc:Returns the retention time (in seconds)
        void setRT(double) except + nogil  # wrap-doc:Sets the retention time (in seconds)

        DriftTimeUnit getDriftTimeUnit() except + nogil  # wrap-doc:Returns the ion mobility drift time unit
        String getDriftTimeUnitAsString() except + nogil  # wrap-doc:Returns the ion mobility drift time unit as string
        void setDriftTimeUnit(DriftTimeUnit dt) except + nogil  # wrap-doc:Sets the ion mobility drift time unit

        Size size() except + nogil  # wrap-doc:Returns the number of peaks in the mobilogram
        void reserve(size_t n) except + nogil  
        void resize(size_t n) except + nogil  # wrap-doc:Resize the peak array 

        MobilityPeak1D & operator[](size_t) except + nogil 

        void updateRanges() except + nogil 
        void clear() except + nogil 
            # wrap-doc:
                #  Clears all data and ranges
                #  
                #  Will delete (clear) all peaks contained in the mobilogram

        void push_back(MobilityPeak1D)  except + nogil  # wrap-doc:Append a peak

        bool isSorted() except + nogil  # wrap-doc:Checks if all peaks are sorted with respect to ascending mobility

        void sortByIntensity(bool reverse) except + nogil 
            # wrap-doc:
                #  Lexicographically sorts the peaks by their intensity
                #  
                #  
                #  Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly

        void sortByPosition() except + nogil 
            # wrap-doc:
                #  Lexicographically sorts the peaks by their position (mobility)
                #  
                #  
                #  The mobilogram is sorted with respect to position (mobility). Meta data arrays will be sorted accordingly

        int findNearest(double) except + nogil 
            # wrap-doc:
                #  Binary search for the peak nearest to a specific mobility
                #  :note: Make sure the mobilogram is sorted with respect to mobility! Otherwise the result is undefined
                #  
                #  
                #  :param mb: The searched for mobility value
                #  :return: Returns the index of the peak.
                #  :raises:
                #    Exception: Precondition is thrown if the mobilogram is empty (not only in debug mode)

        void assign(libcpp_vector[MobilityPeak1D].iterator, libcpp_vector[MobilityPeak1D].iterator) except + nogil  # wrap-ignore
        libcpp_vector[MobilityPeak1D].iterator begin() except + nogil   # wrap-iter-begin:__iter__(MobilityPeak1D)
        libcpp_vector[MobilityPeak1D].iterator end()   except + nogil   # wrap-iter-end:__iter__(MobilityPeak1D)

        libcpp_vector[FloatDataArray] getFloatDataArrays() except + nogil  # wrap-doc:Returns a reference to the float meta data arrays
        libcpp_vector[IntegerDataArray] getIntegerDataArrays() except + nogil  # wrap-doc:Returns a reference to the integer meta data arrays
        libcpp_vector[StringDataArray] getStringDataArrays() except + nogil  # wrap-doc:Returns a reference to the string meta data arrays

        void setFloatDataArrays(libcpp_vector[FloatDataArray] fda) except + nogil  # wrap-doc:Sets the float meta data arrays
        void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida) except + nogil  # wrap-doc:Sets the integer meta data arrays
        void setStringDataArrays(libcpp_vector[StringDataArray] sda) except + nogil  # wrap-doc:Sets the string meta data arrays
