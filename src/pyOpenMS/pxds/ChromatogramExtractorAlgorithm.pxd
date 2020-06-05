from Types cimport *
from String cimport *
from OpenSwathDataStructures cimport *
from SpectrumAccessOpenMS cimport *
from SpectrumAccessOpenMSCached cimport *
from SpectrumAccessOpenMSInMemory cimport *
from SpectrumAccessQuadMZTransforming cimport *
from ISpectrumAccess cimport *
from ProgressLogger cimport *

# typedef boost::shared_ptr<ISpectrumAccess> SpectrumAccessPtr;

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>" namespace "OpenMS":

    cdef cppclass ChromatogramExtractorAlgorithm(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        ChromatogramExtractorAlgorithm() nogil except +
        ChromatogramExtractorAlgorithm(ChromatogramExtractorAlgorithm) nogil except +

        # abstract base class ISpectrumAccess given as first input arg
        void extractChromatograms(
            shared_ptr[ SpectrumAccessOpenMS ] input,
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output,
            libcpp_vector[ ExtractionCoordinates ] extraction_coordinates,
            double mz_extraction_window,
            bool ppm,
            double im_extraction_window,
            String filter) nogil except +

        void extractChromatograms(
            shared_ptr[ SpectrumAccessOpenMSCached ] input,
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output,
            libcpp_vector[ ExtractionCoordinates ] extraction_coordinates,
            double mz_extraction_window,
            bool ppm,
            double im_extraction_window,
            String filter) nogil except +

        void extractChromatograms(
            shared_ptr[ SpectrumAccessOpenMSInMemory ] input,
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output,
            libcpp_vector[ ExtractionCoordinates ] extraction_coordinates,
            double mz_extraction_window,
            bool ppm,
            double im_extraction_window,
            String filter) nogil except +

        void extractChromatograms(
            shared_ptr[ SpectrumAccessQuadMZTransforming ] input,
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output,
            libcpp_vector[ ExtractionCoordinates ] extraction_coordinates,
            double mz_extraction_window,
            bool ppm,
            double im_extraction_window,
            String filter) nogil except +


        # void extract_value_tophat # -> uses iterators

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>" namespace "OpenMS::ChromatogramExtractorAlgorithm":

    cdef cppclass ExtractionCoordinates:

        ExtractionCoordinates() nogil except +
        ExtractionCoordinates(ExtractionCoordinates) nogil except +

        double mz # mz around which should be extracted
        double mz_precursor # precursor m/z value (is currently ignored by the algorithm)
        double rt_start # rt start of extraction (in seconds)
        double rt_end # rt end of extraction (in seconds)
        double ion_mobility # ion mobility value around which should be extracted
        libcpp_string id # identifier

