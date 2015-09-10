from Types cimport *
from Param cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>" namespace "OpenMS":

    cdef cppclass TransformationModelInterpolated(TransformationModel):
        # wrap-inherits:
        #   TransformationModel

        void getDefaultParameters(Param &)
        double evaluate(double value) nogil except +

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>" namespace "OpenMS:TransformationModelInterpolated":

    cdef cppclass Interpolator:
        # wrap-ignore
        # ABSTRACT

      void init(libcpp_vector[double] x, libcpp_vector[double] y)

      double eval(double x) 

