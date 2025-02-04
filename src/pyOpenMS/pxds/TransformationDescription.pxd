from Types cimport *
from Param cimport *
from TransformationModel cimport *
from String cimport *
from StringList cimport *
from TransformationModel cimport TM_DataPoint

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS::TransformationDescription":
    
    cdef cppclass TransformationStatistics "OpenMS::TransformationDescription::TransformationStatistics":
        TransformationStatistics() except + nogil 
        TransformationStatistics(TransformationStatistics &) except + nogil 
        # libcpp_vector[ size_t ] percents # const
        double xmin
        double xmax
        double ymin
        double ymax
        libcpp_map[size_t, double ] percentiles_before
        libcpp_map[size_t, double ] percentiles_after

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS":

    cdef cppclass TransformationDescription:
        TransformationDescription() except + nogil 
        TransformationDescription(TransformationDescription &) except + nogil 

        libcpp_vector[TM_DataPoint] getDataPoints() except + nogil  # wrap-doc:Returns the data points
        void setDataPoints(libcpp_vector[TM_DataPoint]& data) except + nogil  # wrap-doc:Sets the data points. Removes the model that was previously fitted to the data (if any)
        void setDataPoints(libcpp_vector[libcpp_pair[double,double]]& data) except + nogil  # wrap-doc:Sets the data points (backwards-compatible overload). Removes the model that was previously fitted to the data (if any)
        double apply(double) except + nogil  # wrap-doc:Applies the transformation to `value`

        void fitModel(String model_type, Param params)  except + nogil  # wrap-doc:Fits a model to the data
        void fitModel(String model_type)  except + nogil  # wrap-doc:Fits a model to the data

        String getModelType()  except + nogil  # wrap-doc:Gets the type of the fitted model
        Param getModelParameters() except + nogil  # wrap-doc:Returns the model parameters

        void invert() except + nogil  # wrap-doc:Computes an (approximate) inverse of the transformation

        void getDeviations(libcpp_vector[double]& diffs, bool do_apply, bool do_sort) except + nogil 
        # wrap-doc:
                #  Get the deviations between the data pairs
                #  
                #  :param diffs: Output
                #  :param do_apply: Get deviations after applying the model?
                #  :param do_sort: Sort `diffs` before returning?

        TransformationStatistics getStatistics() except + nogil 

        # NAMESPACE # void printSummary(std::ostream & os) except + nogil 

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>" namespace "OpenMS::TransformationDescription":

    void getModelTypes(StringList result) except + nogil  # wrap-attach:TransformationDescription
