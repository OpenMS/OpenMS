from smart_ptr cimport shared_ptr
from ProgressLogger cimport *
from OpenSwathDataStructures cimport *
from SpectrumAccessOpenMS cimport *
from ISpectrumAccess cimport *
from libcpp cimport bool

# typedef boost::shared_ptr<ISpectrumAccess> SpectrumAccessPtr;

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>" namespace "OpenMS":

    cdef cppclass ChromatogramExtractorAlgorithm(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        ChromatogramExtractorAlgorithm() nogil except +
        ChromatogramExtractorAlgorithm(ChromatogramExtractorAlgorithm) nogil except + 

        void extractChromatograms(
            shared_ptr[ ISpectrumAccess ] input,
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output, 
            libcpp_vector[ ExtractionCoordinates ] extraction_coordinates, 
            double mz_extraction_window,
            bool ppm, String filter) nogil except + # wrap-ignore

        # void extract_value_tophat # -> uses iterators

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>" namespace "OpenMS::ChromatogramExtractorAlgorithm":

    cdef cppclass ExtractionCoordinates:

        ExtractionCoordinates() nogil except +
        ExtractionCoordinates(ExtractionCoordinates) nogil except + 

        double mz # mz around which should be extracted
        double rt_start # rt start of extraction (in seconds)
        double rt_end # rt end of extraction (in seconds)
        libcpp_string id # identifier

