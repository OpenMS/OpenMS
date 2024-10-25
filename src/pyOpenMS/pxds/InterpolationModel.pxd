# from BaseModel cimport *
from LinearInterpolation cimport *
from DPosition cimport *

ctypedef double IntensityType
ctypedef DPosition1 PositionType
ctypedef double CoordinateType
# ctypedef Math::LinearInterpolation<double> LinearInterpolation

cdef extern from "<OpenMS/FEATUREFINDER/InterpolationModel.h>" namespace "OpenMS":
    
    cdef cppclass InterpolationModel "OpenMS::InterpolationModel":
      
        InterpolationModel() except + nogil  # wrap-doc:Abstract class for 1D-models that are approximated using linear interpolation
        InterpolationModel(InterpolationModel &) except + nogil 
        
        # double getIntensity(DPosition1 &pos) except + nogil 
        double getIntensity(double coord) except + nogil  # wrap-doc:Access model predicted intensity at position 'pos'
        double getScalingFactor() except + nogil  # wrap-doc:Returns the interpolation class
        void setOffset(double offset) except + nogil  # wrap-doc:Sets the offset of the model
        double getCenter() except + nogil  # wrap-doc:Returns the "center" of the model, particular definition (depends on the derived model)
        void setSamples() except + nogil  # wrap-doc:Sets sample/supporting points of interpolation wrt params
        void setInterpolationStep(double interpolation_step) except + nogil  # wrap-doc:Sets the interpolation step for the linear interpolation of the model
        void setScalingFactor(double scaling) except + nogil  # wrap-doc:Sets the scaling factor of the model

        LinearInterpolation[double,double] getInterpolation() except + nogil  # wrap-doc:Returns the interpolation class
        # void getSamples(SamplesType &cont) except + nogil 
        # typedef typename DPeak<D>::Type PeakType;
        # typedef std::vector<PeakType> SamplesType;
