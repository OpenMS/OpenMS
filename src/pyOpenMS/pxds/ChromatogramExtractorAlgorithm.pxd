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
        ChromatogramExtractorAlgorithm(ChromatogramExtractorAlgorithm &) nogil except + # compiler

        # abstract base class ISpectrumAccess given as first input arg
        void extractChromatograms(
            shared_ptr[ SpectrumAccessOpenMS ] input,
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output,
            libcpp_vector[ ExtractionCoordinates ] extraction_coordinates,
            double mz_extraction_window,
            bool ppm,
            double im_extraction_window,
            String filter) nogil except +
            # wrap-doc:
            #     Extract chromatograms at the m/z and RT defined by the ExtractionCoordinates
            #     -----
            #     :param: input Input spectral map
            #     :param output: Output chromatograms (XICs)
            #     :param extraction_coordinates: Extracts around these coordinates (from
            #      rt_start to rt_end in seconds - extracts the whole chromatogram if
            #      rt_end - rt_start < 0).
            #     :param mz_extraction_window: Extracts a window of this size in m/z
            #     dimension in Th or ppm (e.g. a window of 50 ppm means an extraction of
            #     25 ppm on either side)
            #     :param ppm: Whether mz_extraction_window is in ppm or in Th
            #     :param filter: Which function to apply in m/z space (currently "tophat" only)

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
