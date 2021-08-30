from Types cimport *
from Param cimport *
from TransformationModel cimport *
from String cimport *
from StringList cimport *
from TransformationModel cimport TM_DataPoint

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS::TransformationDescription":
    
    cdef cppclass TransformationStatistics "OpenMS::TransformationDescription::TransformationStatistics":
        TransformationStatistics() nogil except +
        TransformationStatistics(TransformationStatistics &) nogil except +
        # libcpp_vector[ size_t ] percents # const
        double xmin
        double xmax
        double ymin
        double ymax
        libcpp_map[size_t, double ] percentiles_before
        libcpp_map[size_t, double ] percentiles_after

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS":

    cdef cppclass TransformationDescription:
        TransformationDescription() nogil except +
        TransformationDescription(TransformationDescription &) nogil except +

        libcpp_vector[TM_DataPoint] getDataPoints() nogil except + # wrap-doc:Returns the data points
        void setDataPoints(libcpp_vector[TM_DataPoint]& data) nogil except + # wrap-doc:Sets the data points. Removes the model that was previously fitted to the data (if any)
        void setDataPoints(libcpp_vector[libcpp_pair[double,double]]& data) nogil except + # wrap-doc:Sets the data points (backwards-compatible overload). Removes the model that was previously fitted to the data (if any)
        double apply(double) nogil except + # wrap-doc:Applies the transformation to `value`

        void fitModel(String model_type, Param params)  nogil except + # wrap-doc:Fits a model to the data
        void fitModel(String model_type)  nogil except + # wrap-doc:Fits a model to the data

        String getModelType()  nogil except + # wrap-doc:Gets the type of the fitted model
        Param getModelParameters() nogil except + # wrap-doc:Returns the model parameters

        void invert() nogil except + # wrap-doc:Computes an (approximate) inverse of the transformation

        void getDeviations(libcpp_vector[double]& diffs, bool do_apply, bool do_sort) nogil except +
        # wrap-doc:
                #   Get the deviations between the data pairs
                #   -----
                #   :param diffs: Output
                #   :param do_apply: Get deviations after applying the model?
                #   :param do_sort: Sort `diffs` before returning?

        TransformationStatistics getStatistics() nogil except +

        # NAMESPACE # void printSummary(std::ostream & os) nogil except +

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS::TransformationDescription":

    void getModelTypes(StringList result) nogil except + # wrap-attach:TransformationDescription
