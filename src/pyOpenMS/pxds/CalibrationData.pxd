from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CalibrationData.h>" namespace "OpenMS":

    cdef cppclass CalibrationData:

        CalibrationData()  except + nogil 

        CalibrationData(CalibrationData &) except + nogil  # compiler
        double getMZ(Size) except + nogil  # wrap-doc:Retrieve the observed m/z of the i'th calibration point
        double getRT(Size) except + nogil  # wrap-doc:Retrieve the observed RT of the i'th calibration point
        double getIntensity(Size) except + nogil  # wrap-doc:Retrieve the intensity of the i'th calibration point
        Size size() except + nogil  # wrap-doc:Number of calibration points
        bool empty() except + nogil  # wrap-doc:Returns `True` if there are no peaks
        void clear() except + nogil  # wrap-doc:Remove all calibration points

        void setUsePPM(bool) except + nogil 
        bool usePPM() except + nogil  # wrap-doc:Current error unit (ppm or Th)
        void insertCalibrationPoint(double rt, double mz_obs, float intensity, double mz_ref, double weight, int group) except + nogil 
        Size getNrOfGroups() except + nogil  # wrap-doc:Number of peak groups (can be 0)
        double getError(Size) except + nogil  # wrap-doc:Retrieve the error for i'th calibrant in either ppm or Th (depending on usePPM())
        double getRefMZ(Size) except + nogil  # wrap-doc:Retrieve the theoretical m/z of the i'th calibration point
        double getWeight(Size) except + nogil  # wrap-doc:Retrieve the weight of the i'th calibration point
        int getGroup(Size i) except + nogil  # wrap-doc:Retrieve the group of the i'th calibration point
        CalibrationData median(double, double) except + nogil  # wrap-doc:Compute the median in the given RT range for every peak group
        void sortByRT() except + nogil  # wrap-doc:Sort calibration points by RT, to allow for valid RT chunking

cdef extern from "<OpenMS/DATASTRUCTURES/CalibrationData.h>" namespace "OpenMS::CalibrationData":

        # static members
        StringList getMetaValues() except + nogil  # wrap-attach:CalibrationData
