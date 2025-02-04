from Types cimport *
from ChromatogramSettings cimport *
from ChromatogramPeak cimport *
from String cimport *
from RangeManager cimport *
from MetaInfoDescription cimport *
from DataArrays cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/MSChromatogram.h>" namespace "OpenMS":

    cdef cppclass MSChromatogram (ChromatogramSettings, RangeManagerRtInt):
        # wrap-inherits:
        #  ChromatogramSettings
        #  RangeManagerRtInt
        #
        # wrap-doc:
        #  The representation of a chromatogram.
        #  Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
        #  Iterations yields access to underlying peak objects but is slower
        #  Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
        #  See help(ChromatogramSettings) for information about meta-information
        #  
        #  Usage:
        #
        #  .. code-block:: python
        #  
        #    precursor = chromatogram.getPrecursor()
        #    product = chromatogram.getProduct()
        #    rt, intensities = chromatogram.get_peaks()
        #  

        MSChromatogram() except + nogil 
        MSChromatogram(MSChromatogram &) except + nogil 
        double getMZ() except + nogil  # wrap-doc:Returns the mz of the product entry, makes sense especially for MRM scans
        # void   setMZ(double) except + nogil 

        String getName() except + nogil  # wrap-doc:Returns the name
        void setName(String) except + nogil  # wrap-doc:Sets the name

        Size size() except + nogil 
        void reserve(size_t n) except + nogil  
        void resize(size_t n) except + nogil  # wrap-doc:Resize the peak array 

        ChromatogramPeak & operator[](size_t) except + nogil 

        void updateRanges() except + nogil 
        void clear(int) except + nogil 
            # wrap-doc:
                #  Clears all data and meta data
                #  
                #  
                #  :param clear_meta_data: If true, all meta data is cleared in addition to the data

        void push_back(ChromatogramPeak)  except + nogil  # wrap-doc:Append a peak

        bool isSorted() except + nogil  # wrap-doc:Checks if all peaks are sorted with respect to ascending RT

        void sortByIntensity(bool reverse) except + nogil 
            # wrap-doc:
                #  Lexicographically sorts the peaks by their intensity
                #  
                #  
                #  Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly

        void sortByPosition() except + nogil 
            # wrap-doc:
                #  Lexicographically sorts the peaks by their position
                #  
                #  
                #  The chromatogram is sorted with respect to position. Meta data arrays will be sorted accordingly

        int findNearest(double) except + nogil 
            # wrap-doc:
                #  Binary search for the peak nearest to a specific RT
                #  :note: Make sure the chromatogram is sorted with respect to RT! Otherwise the result is undefined
                #  
                #  
                #  :param rt: The searched for mass-to-charge ratio searched
                #  :return: Returns the index of the peak.
                #  :raises:
                #    Exception: Precondition is thrown if the chromatogram is empty (not only in debug mode)

        void assign(libcpp_vector[ChromatogramPeak].iterator, libcpp_vector[ChromatogramPeak].iterator) except + nogil  # wrap-ignore
        libcpp_vector[ChromatogramPeak].iterator begin() except + nogil   # wrap-iter-begin:__iter__(ChromatogramPeak)
        libcpp_vector[ChromatogramPeak].iterator end()   except + nogil   # wrap-iter-end:__iter__(ChromatogramPeak)

        libcpp_vector[FloatDataArray] getFloatDataArrays() except + nogil  # wrap-doc:Returns a reference to the float meta data arrays
        libcpp_vector[IntegerDataArray] getIntegerDataArrays() except + nogil  # wrap-doc:Returns a reference to the integer meta data arrays
        libcpp_vector[StringDataArray] getStringDataArrays() except + nogil  # wrap-doc:Returns a reference to the string meta data arrays

        void setFloatDataArrays(libcpp_vector[FloatDataArray] fda) except + nogil  # wrap-doc:Sets the float meta data arrays
        void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida) except + nogil  # wrap-doc:Sets the integer meta data arrays
        void setStringDataArrays(libcpp_vector[StringDataArray] sda) except + nogil  # wrap-doc:Sets the string meta data arrays
