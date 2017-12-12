from Types cimport *
from TransformationModel cimport *
from TransformationModelInterpolated cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>" namespace "OpenMS":

    cdef cppclass TransformationModelLowess(TransformationModel) :
        # wrap-inherits:
        #  TransformationModel

        TransformationModelLowess() nogil except + #wrap-ignore
        TransformationModelLowess(TransformationModelLowess) nogil except + #wrap-ignore
        TransformationModelLowess(libcpp_vector[TM_DataPoint]& data, Param& params) nogil except +

        double evaluate(double value) nogil except +
        void getDefaultParameters(Param& params) nogil except +
