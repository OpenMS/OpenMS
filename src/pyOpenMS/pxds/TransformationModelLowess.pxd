from Types cimport *
from TransformationModel cimport *
from TransformationModelInterpolated cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>" namespace "OpenMS":

    cdef cppclass TransformationModelLowess(TransformationModel) :
        # wrap-inherits:
        #  TransformationModel

        # copy constructor of 'TransformationModelLowess' is implicitly deleted because base class 'OpenMS::TransformationModel' has an inaccessible copy constructor public TransformationModel
        TransformationModelLowess(TransformationModelLowess &) nogil except + # wrap-ignore
        TransformationModelLowess(libcpp_vector[TM_DataPoint]& data, Param& params) nogil except +
        void getDefaultParameters(Param &) nogil except +
        double evaluate(double value) nogil except +

# COMMENT: wrap static methods
#cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>" namespace "OpenMS::TransformationModelLowess":     
#    void getDefaultParameters(Param& params) nogil except + #wrap-attach:TransformationModelLowess
