from Types cimport *
from Param cimport *
from String cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as String

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>" namespace "OpenMS":


    cdef cppclass TransformationModelLinear(TransformationModel):
        # wrap-inherits:
        #  TransformationModel

        # copy constructor of 'TransformationModelLinear' is implicitly deleted because base class 'OpenMS::TransformationModel' has an inaccessible copy constructor public TransformationModel
        TransformationModelLinear(TransformationModelLinear &) except + nogil  # wrap-ignore
        TransformationModelLinear(libcpp_vector[TM_DataPoint]& data, Param& params) except + nogil 

        double evaluate(double value) except + nogil 
        # void getParameters(double & slope, double & intercept, String& x_weight, String& y_weight, double & x_datum_min, double & x_datum_max, double & y_datum_min, double & y_datum_max) except + nogil 
        void invert() except + nogil 

# COMMENT: wrap static methods
cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>" namespace "OpenMS::TransformationModelLinear":
        
        # static members
        void getDefaultParameters(Param &) except + nogil  # wrap-attach:TransformationModelLinear

