from Types cimport *
from Param cimport *

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>" namespace "OpenMS":

    cdef cppclass TransformationModelInterpolated(TransformationModel):
        # wrap-inherits:
        #   TransformationModel

        # copy constructor of 'TransformationModelInterpolated' is implicitly deleted because base class 'OpenMS::TransformationModel' has an inaccessible copy constructor public TransformationModel
        TransformationModelInterpolated(TransformationModelInterpolated &) nogil except + #wrap-ignore
        TransformationModelInterpolated(libcpp_vector[TM_DataPoint]& data, Param& params) nogil except +

        void getDefaultParameters(Param &) nogil except + # wrap-doc:Gets the default parameters
        double evaluate(double value) nogil except +
        # wrap-doc:
                #   Evaluate the interpolation model at the given value
                #   -----
                #   :param value: The position where the interpolation should be evaluated
                #   :returns: The interpolated value

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>" namespace "OpenMS::TransformationModelInterpolated":

    cdef cppclass Interpolator:
        # wrap-ignore
        # ABSTRACT
        # no-pxd-import

      void init(libcpp_vector[double] x, libcpp_vector[double] y)

      double eval(double x)
