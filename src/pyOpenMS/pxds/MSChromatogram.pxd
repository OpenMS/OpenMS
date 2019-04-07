from Types cimport *
from ChromatogramSettings cimport *
from ChromatogramPeak cimport *
from String cimport *
from RangeManager cimport *
from MetaInfoDescription cimport *
from DataArrays cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/MSChromatogram.h>" namespace "OpenMS":

    cdef cppclass MSChromatogram (ChromatogramSettings, RangeManager1):
        # wrap-inherits:
        #  ChromatogramSettings
        #  RangeManager1
        #
        # wrap-doc:
        #   The representation of a chromatogram.
        #   Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
        #   Iterations yields access to underlying peak objects but is slower
        #   Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
        #   See help(ChromatogramSettings) for information about meta-information

        MSChromatogram() nogil except +
        MSChromatogram(MSChromatogram &) nogil except +
        double getMZ() nogil except + #wrap-doc:returns the mz of the product entry, makes sense especially for MRM scans
        # void   setMZ(double) nogil except +

        String getName() nogil except +
        void setName(String) nogil except +

        Size size() nogil except +
        void reserve(size_t n) nogil except + 

        ChromatogramPeak operator[](int) nogil except +

        void updateRanges() nogil except +
        void clear(int) nogil except +
        void push_back(ChromatogramPeak)  nogil except + #wrap-doc:Append a peak

        bool isSorted() nogil except +

        void sortByIntensity(bool reverse) nogil except +
        void sortByPosition() nogil except +

        int findNearest(double) nogil except+

        void assign(libcpp_vector[ChromatogramPeak].iterator, libcpp_vector[ChromatogramPeak].iterator) nogil except + # wrap-ignore
        libcpp_vector[ChromatogramPeak].iterator begin() nogil except +  # wrap-iter-begin:__iter__(ChromatogramPeak)
        libcpp_vector[ChromatogramPeak].iterator end()   nogil except +  # wrap-iter-end:__iter__(ChromatogramPeak)

        libcpp_vector[FloatDataArray] getFloatDataArrays() nogil except +
        libcpp_vector[IntegerDataArray] getIntegerDataArrays() nogil except +
        libcpp_vector[StringDataArray] getStringDataArrays() nogil except +

        void setFloatDataArrays(libcpp_vector[FloatDataArray] fda) nogil except +
        void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida) nogil except +
        void setStringDataArrays(libcpp_vector[StringDataArray] sda) nogil except +

