from Types cimport *
from TransformationModel cimport *
from TransformationModelInterpolated cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>" namespace "OpenMS":

    cdef cppclass TransformationModelLowess(TransformationModel) :
        # wrap-inherits:
        #  TransformationModel

        # copy constructor of 'TransformationModelLowess' is implicitly deleted because base class 'OpenMS::TransformationModel' has an inaccessible copy constructor public TransformationModel
        TransformationModelLowess(TransformationModelLowess &) except + nogil  # wrap-ignore
        TransformationModelLowess(libcpp_vector[TM_DataPoint]& data, Param& params) except + nogil
        double evaluate(double value) except + nogil

        @staticmethod
        void getDefaultParameters(Param& params) except + nogil
