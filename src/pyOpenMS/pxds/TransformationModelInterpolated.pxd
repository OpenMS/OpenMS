from Types cimport *
from Param cimport *

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>" namespace "OpenMS":

    cdef cppclass TransformationModelInterpolated(TransformationModel):
        # wrap-inherits:
        #   TransformationModel

        TransformationModelInterpolated(TransformationModelInterpolated) nogil except + #wrap-ignore
        TransformationModelInterpolated(libcpp_vector[TM_DataPoint]& data, Param& params) nogil except +

        void getDefaultParameters(Param &)
        double evaluate(double value) nogil except +

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>" namespace "OpenMS:TransformationModelInterpolated":

    cdef cppclass Interpolator:
        # wrap-ignore
        # ABSTRACT
        # no-pxd-import

      void init(libcpp_vector[double] x, libcpp_vector[double] y)

      double eval(double x)
