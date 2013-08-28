from Types cimport *
from MassTrace cimport *
from Feature cimport *
from FeatureMap cimport *
from MzTab cimport *
from DefaultParamHandler cimport *
from String cimport *
from ProgressLogger cimport *
from AccurateMassSearchResult cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>" namespace "OpenMS":
    
    cdef cppclass AccurateMassSearchEngine(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        AccurateMassSearchEngine() nogil except +
        AccurateMassSearchEngine(AccurateMassSearchEngine) nogil except + #wrap-ignore
        # TODO : not implemented
        # AccurateMassSearchEngine(String & map_fname) nogil except +
        # void queryByMass(DoubleReal & , DoubleReal & , libcpp_vector[ AccurateMassSearchResult ] & ) nogil except +
        # void queryByFeature(Feature & , libcpp_vector[ AccurateMassSearchResult ] & ) nogil except +
        void run(FeatureMap[Feature] & , MzTab & ) nogil except +

