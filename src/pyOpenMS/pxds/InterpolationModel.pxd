# from BaseModel cimport *
from LinearInterpolation cimport *
from DPosition cimport *

ctypedef double IntensityType
ctypedef DPosition1 PositionType
ctypedef double CoordinateType
# ctypedef Math::LinearInterpolation<double> LinearInterpolation

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>" namespace "OpenMS":
    
    cdef cppclass InterpolationModel "OpenMS::InterpolationModel":
      
        InterpolationModel() nogil except + # wrap-doc:Abstract class for 1D-models that are approximated using linear interpolation
        InterpolationModel(InterpolationModel &) nogil except +
        
        # double getIntensity(DPosition1 &pos) nogil except +
        double getIntensity(double coord) nogil except + # wrap-doc:Access model predicted intensity at position 'pos'
        double getScalingFactor() nogil except + # wrap-doc:Returns the interpolation class
        void setOffset(double offset) nogil except + # wrap-doc:Sets the offset of the model
        double getCenter() nogil except + # wrap-doc:Returns the "center" of the model, particular definition (depends on the derived model)
        void setSamples() nogil except + # wrap-doc:Sets sample/supporting points of interpolation wrt params
        void setInterpolationStep(double interpolation_step) nogil except + # wrap-doc:Sets the interpolation step for the linear interpolation of the model
        void setScalingFactor(double scaling) nogil except + # wrap-doc:Sets the scaling factor of the model

        LinearInterpolation[double,double] getInterpolation() nogil except + # wrap-doc:Returns the interpolation class
        # void getSamples(SamplesType &cont) nogil except +
        # typedef typename DPeak<D>::Type PeakType;
        # typedef std::vector<PeakType> SamplesType;
