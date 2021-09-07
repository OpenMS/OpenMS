from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CalibrationData.h>" namespace "OpenMS":

    cdef cppclass CalibrationData:

        CalibrationData()  nogil except +

        CalibrationData(CalibrationData &) nogil except + # compiler
        double getMZ(Size) nogil except + # wrap-doc:Retrieve the observed m/z of the i'th calibration point
        double getRT(Size) nogil except + # wrap-doc:Retrieve the observed RT of the i'th calibration point
        double getIntensity(Size) nogil except + # wrap-doc:Retrieve the intensity of the i'th calibration point
        Size size() nogil except + # wrap-doc:Number of calibration points
        bool empty() nogil except + # wrap-doc:Returns `True` if there are no peaks
        void clear() nogil except + # wrap-doc:Remove all calibration points

        void setUsePPM(bool) nogil except +
        bool usePPM() nogil except + # wrap-doc:Current error unit (ppm or Th)
        void insertCalibrationPoint(double rt, double mz_obs, float intensity, double mz_ref, double weight, int group) nogil except +
        Size getNrOfGroups() nogil except + # wrap-doc:Number of peak groups (can be 0)
        double getError(Size) nogil except + # wrap-doc:Retrieve the error for i'th calibrant in either ppm or Th (depending on usePPM())
        double getRefMZ(Size) nogil except + # wrap-doc:Retrieve the theoretical m/z of the i'th calibration point
        double getWeight(Size) nogil except + # wrap-doc:Retrieve the weight of the i'th calibration point
        int getGroup(Size i) nogil except + # wrap-doc:Retrieve the group of the i'th calibration point
        CalibrationData median(double, double) nogil except + # wrap-doc:Compute the median in the given RT range for every peak group
        void sortByRT() nogil except + # wrap-doc:Sort calibration points by RT, to allow for valid RT chunking

cdef extern from "<OpenMS/DATASTRUCTURES/CalibrationData.h>" namespace "OpenMS::CalibrationData":

        # static members
        StringList getMetaValues() nogil except + # wrap-attach:CalibrationData
