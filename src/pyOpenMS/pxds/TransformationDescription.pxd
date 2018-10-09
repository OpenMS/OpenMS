from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from Param cimport *
from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS::TransformationDescription":
    
    cdef cppclass TransformationStatistics "OpenMS::TransformationDescription::TransformationStatistics":
        TransformationStatistics(TransformationStatistics) nogil except + #wrap-ignore
        libcpp_vector[ size_t ] percents
        double xmin
        double xmax
        double ymin
        double ymax
        libcpp_map[ Size, double ] percentiles_before
        libcpp_map[ Size, double ] percentiles_after

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS":
    
    cdef cppclass TransformationDescription "OpenMS::TransformationDescription":
        TransformationDescription() nogil except +
        TransformationDescription(TransformationDescription) nogil except +
        TransformationDescription(const DataPoints & data) nogil except +
        void fitModel(const String & model_type, const Param & params) nogil except +
        double apply(double value) nogil except +
        String  getModelType() nogil except +
        void setDataPoints(const DataPoints & data) nogil except +
        void setDataPoints(const libcpp_vector[ libcpp_pair[ double, double ] ] & data) nogil except +
        DataPoints  getDataPoints() nogil except +
        Param  getModelParameters() nogil except +
        void invert() nogil except +
        void getDeviations(libcpp_vector[ double ] & diffs, bool do_apply, bool do_sort) nogil except +
        TransformationStatistics getStatistics() nogil except +
        # NAMESPACE # void printSummary(std::ostream & os) nogil except +
        void getModelTypes(StringList & result) nogil except +














from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair   cimport pair as libcpp_pair
from Param cimport *
from String cimport *
from StringList cimport *
from TransformationModel cimport TM_DataPoint

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS":

    cdef cppclass TransformationDescription:
        TransformationDescription() nogil except +
        TransformationDescription(TransformationDescription) nogil except + # wrap-ignore
        libcpp_vector[TM_DataPoint] getDataPoints() nogil except +
        void setDataPoints(libcpp_vector[TM_DataPoint]& data) nogil except +
        # void setDataPoints(libcpp_vector[libcpp_pair[double,double]]& data) nogil except +
        double apply(double) nogil except +

        void fitModel(String model_type, Param params)  nogil except +
        void fitModel(String model_type)  nogil except +

        String getModelType()  nogil except +
        Param getModelParameters() nogil except +

        void invert() nogil except +

        void getDeviations(libcpp_vector[double]& diffs, bool do_apply, bool do_sort) nogil except +
        TransformationStatistics getStatistics() nogil except +

        # NAMESPACE # void printSummary(std::ostream & os) nogil except +

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS::TransformationDescription":

    void getModelTypes(StringList result) nogil except + # wrap-attach:TransformationDescription

