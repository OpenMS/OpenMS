from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CalibrationData.h>" namespace "OpenMS":

    cdef cppclass CalibrationData:

        CalibrationData()  nogil except +
        CalibrationData(CalibrationData &) nogil except +
        double getMZ(Size) nogil except +
        double getRT(Size) nogil except +
        double getIntensity(Size) nogil except +
        Size size() nogil except +
        bool empty() nogil except +
        void clear() nogil except +
        void setUsePPM(bool) nogil except +
        bool usePPM() nogil except +
        void insertCalibrationPoint(double rt, double mz_obs, float intensity, double mz_ref, double weight, int group) nogil except +
        Size getNrOfGroups() nogil except +
        double getError(Size) nogil except +
        double getRefMZ(Size) nogil except +
        double getWeight(Size) nogil except +
        int getGroup(Size i) nogil except +
        CalibrationData median(double, double) nogil except +
        void sortByRT() nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/CalibrationData.h>" namespace "OpenMS::CalibrationData":

        # static members
        StringList getMetaValues() nogil except + # wrap-attach:CalibrationData

