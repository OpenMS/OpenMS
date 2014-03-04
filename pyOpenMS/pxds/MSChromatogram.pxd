from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from ChromatogramSettings cimport *
from MetaInfoInterface cimport *
from ChromatogramPeak cimport *
from String cimport *
from RangeManager cimport *

cdef extern from "<OpenMS/KERNEL/MSChromatogram.h>" namespace "OpenMS":

    cdef cppclass MSChromatogram[ChromatogramPeakT] (ChromatogramSettings, MetaInfoInterface, RangeManager1):
        # wrap-inherits:
        #  ChromatogramSettings
        #  MetaInfoInterface
        #  RangeManager1

        # wrap-instances:
        #   MSChromatogram := MSChromatogram[ChromatogramPeak]

        MSChromatogram() nogil except +
        MSChromatogram(MSChromatogram) nogil except +
        double getMZ() nogil except +
        # void   setMZ(double) nogil except +

        libcpp_string getName() nogil except +
        void setName(libcpp_string) nogil except +

        Size size() nogil except +
        ChromatogramPeakT operator[](int) nogil except +

        void updateRanges() nogil except +
        void clear(int) nogil except +
        void push_back(ChromatogramPeakT)  nogil except +

        bool isSorted() nogil except +

        void sortByIntensity(bool reverse) nogil except +
        void sortByPosition() nogil except +

        int findNearest(double) nogil except+

        void assign(libcpp_vector[ChromatogramPeak].iterator, libcpp_vector[ChromatogramPeak].iterator) nogil except + # wrap-ignore
        libcpp_vector[ChromatogramPeakT].iterator begin() nogil except +  # wrap-iter-begin:__iter__(ChromatogramPeakT)
        libcpp_vector[ChromatogramPeakT].iterator end()   nogil except +  # wrap-iter-end:__iter__(ChromatogramPeakT)

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

