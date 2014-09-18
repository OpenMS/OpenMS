from Types cimport *
from ProgressLogger cimport *
from ISpectrumAccess cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>" namespace "OpenMS":

    cdef cppclass ChromatogramExtractorAlgorithm(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        ChromatogramExtractorAlgorithm() nogil except +
        ChromatogramExtractorAlgorithm(ChromatogramExtractorAlgorithm) nogil except + 

        ## wont work because of vector<shared_ptr<...> >
        # void extractChromatograms(SpectrumAccessPtr input_, 
        #     libcpp_vector[OSChromatogramPtr]& output, 
        #     libcpp_vector[ExtractionCoordinates] extraction_coordinates,
        #     double mz_extraction_window, bool ppm, String filter_) nogil except +

# no support for nested classes yet in Cython
cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>" namespace "OpenMS::ChromatogramExtractorAlgorithm":

    cdef cppclass ExtractionCoordinates:

        ExtractionCoordinates() nogil except +
        ExtractionCoordinates(ExtractionCoordinates) nogil except +

        # members
        double mz
        double rt_start
        double rt_end
        #libcpp_string id_ # TODO how to name the field

