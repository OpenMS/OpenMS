from Types cimport *
from String cimport *
from OpenSwathDataStructures cimport *
from SpectrumAccessOpenMS cimport *
from SpectrumAccessOpenMSCached cimport *
from SpectrumAccessOpenMSInMemory cimport *
from SpectrumAccessQuadMZTransforming cimport *
from ISpectrumAccess cimport *
from ProgressLogger cimport *

# typedef std::shared_ptr<ISpectrumAccess> SpectrumAccessPtr;

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>" namespace "OpenMS":

    cdef cppclass ChromatogramExtractorAlgorithm(ProgressLogger):
        # wrap-doc:
        #  The ChromatogramExtractorAlgorithm extracts chromatograms from a MS
        #  data. * * It will take as input a set of transitions coordinates and
        #  will extract * the signal of the provided map at the product ion m/z
        #  and retention time * (rt) values specified by the extraction
        #  coordinates. This interface only * expects a set of coordinates which
        #  are up to the user to fill but a * convenient prepare_coordinates
        #  function is provided (in the * ChromatogramExtractor class) to create
        #  the coordinates for the most common * case of an MS2 and MS1
        #  extraction. * * In the case of MS2 extraction, the map is assumed to
        #  originate from a SWATH * (data-independent acquisition or DIA)
        #  experiment. *
        # wrap-inherits:
        #   ProgressLogger

        ChromatogramExtractorAlgorithm() except + nogil 
        ChromatogramExtractorAlgorithm(ChromatogramExtractorAlgorithm &) except + nogil  # compiler

        # abstract base class ISpectrumAccess given as first input arg
        void extractChromatograms(
            shared_ptr[ SpectrumAccessOpenMS ] input,
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output,
            libcpp_vector[ ExtractionCoordinates ] extraction_coordinates,
            double mz_extraction_window,
            bool ppm,
            double im_extraction_window,
            String filter) except + nogil 
            # wrap-doc:
            #    Extract chromatograms at the m/z and RT defined by the ExtractionCoordinates
            #      
            #  
            #  :param input: Input spectral map
            #  :param output: Output chromatograms (XICs)
            #  :param extraction_coordinates: Extracts around these coordinates (from
            #   rt_start to rt_end in seconds - extracts the whole chromatogram if
            #   rt_end - rt_start < 0).
            #  :param mz_extraction_window: Extracts a window of this size in m/z
            #    dimension in Th or ppm (e.g. a window of 50 ppm means an extraction of
            #    25 ppm on either side)
            #  :param ppm: Whether mz_extraction_window is in ppm or in Th
            #  :param filter: Which function to apply in m/z space (currently "tophat" only)

            #  void extract_value_tophat # -> uses iterators

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>" namespace "OpenMS::ChromatogramExtractorAlgorithm":

    cdef cppclass ExtractionCoordinates:

        ExtractionCoordinates() except + nogil 
        ExtractionCoordinates(ExtractionCoordinates) except + nogil 

        double mz # mz around which should be extracted
        double mz_precursor # precursor m/z value (is currently ignored by the algorithm)
        double rt_start # rt start of extraction (in seconds)
        double rt_end # rt end of extraction (in seconds)
        double ion_mobility # ion mobility value around which should be extracted
        libcpp_string id # identifier
