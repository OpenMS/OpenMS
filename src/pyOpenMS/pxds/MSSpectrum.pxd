from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from SpectrumSettings cimport *
from MetaInfoInterface cimport *
from Peak1D cimport *
from String cimport *
from RangeManager cimport *
from DataArrays cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/MSSpectrum.h>" namespace "OpenMS":

    cdef cppclass MSSpectrum[PeakT](SpectrumSettings, MetaInfoInterface, RangeManager1):
        # wrap-inherits:
        #  SpectrumSettings
        #  MetaInfoInterface
        #  RangeManager1

        # wrap-instances:
        #   MSSpectrum := MSSpectrum[Peak1D]
        #   ChromatogramSpectrum := MSSpectrum[ChromatogramPeak]

        # COMMENT: get raw data through get_peaks or by iterating through peaks
        # COMMENT: set raw data through set_peaks

        MSSpectrum() nogil except +
        MSSpectrum(MSSpectrum[PeakT] &) nogil except + 
        double getRT() nogil except +
        void setRT(double) nogil except +
        double getDriftTime() nogil except +
        void setDriftTime(double) nogil except +
        unsigned int getMSLevel() nogil except +
        void setMSLevel(unsigned int) nogil except +

        libcpp_string getName() nogil except +
        void setName(libcpp_string) nogil except +

        Size size() nogil except +

        PeakT operator[](int) nogil except + # wrap-upper-limit:size()

        void updateRanges() nogil except +
        void clear(int) nogil except +
        void push_back(PeakT)  nogil except +

        bool isSorted() nogil except +

        int findNearest(double) nogil except+
        int findNearest(double, double) nogil except+
        int findNearest(double, double, double) nogil except+

        MSSpectrum[PeakT] select(libcpp_vector[ size_t ] & indices) nogil except +

        void assign(libcpp_vector[Peak1D].iterator, libcpp_vector[Peak1D].iterator) nogil except + # wrap-ignore
        libcpp_vector[PeakT].iterator begin() nogil except +  # wrap-iter-begin:__iter__(PeakT)
        libcpp_vector[PeakT].iterator end()   nogil except +  # wrap-iter-end:__iter__(PeakT)

        bool operator==(MSSpectrum[PeakT]) nogil except +
        bool operator!=(MSSpectrum[PeakT]) nogil except +

        void sortByIntensity(bool reverse) nogil except +
        void sortByPosition() nogil except +

        libcpp_vector[FloatDataArray] getFloatDataArrays() nogil except +
        libcpp_vector[IntegerDataArray] getIntegerDataArrays() nogil except +
        libcpp_vector[StringDataArray] getStringDataArrays() nogil except +

        void setFloatDataArrays(libcpp_vector[FloatDataArray] fda) nogil except +
        void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida) nogil except +
        void setStringDataArrays(libcpp_vector[StringDataArray] sda) nogil except +

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

