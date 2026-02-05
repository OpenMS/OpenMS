from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ExperimentalSettings cimport *
from SwathMap cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>" namespace "OpenMS":

    cdef cppclass FullSwathFileConsumer:
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import

        FullSwathFileConsumer() except + nogil  #wrap-ignore
        FullSwathFileConsumer(FullSwathFileConsumer &) except + nogil  # compiler
        FullSwathFileConsumer(libcpp_vector[SwathMap] swath_boundaries) except + nogil 

        void setExpectedSize(Size s, Size c) except + nogil 
        void setExperimentalSettings(ExperimentalSettings exp) except + nogil 

        void retrieveSwathMaps(libcpp_vector[SwathMap]& maps) except + nogil 

        void consumeSpectrum(MSSpectrum & s) except + nogil 
        void consumeChromatogram(MSChromatogram & c) except + nogil 

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>" namespace "OpenMS":

    cdef cppclass RegularSwathFileConsumer(FullSwathFileConsumer):
        # wrap-doc:
        #  Abstract base class which can consume spectra coming from SWATH
        #  experiment stored in a single file. * * The class consumes spectra
        #  which are coming from a complete SWATH * experiment. It will group MS2
        #  spectra by their precursor m/z, assuming * that they correspond to the
        #  same SWATH window. For example, the spectra * could be arranged in the
        #  following fashion: * * - MS1 Spectrum (no precursor) * - MS2 Spectrum
        #  (precursor = [400,425]) * - MS2 Spectrum (precursor = [425,450]) * -
        #  [...] * - MS2 Spectrum (precursor = [1175,1200]) * - MS1 Spectrum (no
        #  precursor) * - MS2 Spectrum (precursor = [400,425]) * - MS2 Spectrum
        #  (precursor = [425,450]) * - [...] * * Base classes are expected to
        #  implement functions consuming a spectrum coming * from a specific
        #  SWATH or an MS1 spectrum and a final function * ensureMapsAreFilled_
        #  after which the swath_maps_ vector needs to contain * valid pointers
        #  to MSExperiment. * * In addition it is possible to provide the swath
        #  boundaries and the read in * spectra will be matched by their
        #  precursor m/z to the "center" attribute * of the provided Swath maps.
        #  * * Usage: * * @code * FullSwathFileConsumer * dataConsumer; * //
        #  assign dataConsumer to an implementation of FullSwathFileConsumer *
        #  MzMLFile().transform(file, dataConsumer); *
        #  dataConsumer->retrieveSwathMaps(maps); * @endcode * */ class
        #  OPENMS_DLLAPI FullSwathFileConsumer : public
        #  Interfaces::IMSDataConsumer {
        # wrap-inherits:
        #   FullSwathFileConsumer

        RegularSwathFileConsumer() except + nogil 
        RegularSwathFileConsumer(RegularSwathFileConsumer &) except + nogil  # compiler


cdef extern from "<OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>" namespace "OpenMS":

    cdef cppclass CachedSwathFileConsumer(FullSwathFileConsumer):
        # wrap-inherits:
        #   FullSwathFileConsumer

        CachedSwathFileConsumer() except + nogil  #wrap-ignore
        CachedSwathFileConsumer(CachedSwathFileConsumer &) except + nogil  # compiler
        CachedSwathFileConsumer(String cachedir, String basename, 
                                Size nr_ms1_spectra, libcpp_vector[int] nr_ms2_spectra)

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MzMLSwathFileConsumer(FullSwathFileConsumer) :
        # wrap-inherits:
        #  FullSwathFileConsumer

        MzMLSwathFileConsumer() except + nogil  # wrap-ignore
        MzMLSwathFileConsumer(MzMLSwathFileConsumer) except + nogil  # compiler
        MzMLSwathFileConsumer(String cachedir,
                              String basename, 
                              Size nr_ms1_spectra, 
                              libcpp_vector[ int ] nr_ms2_spectra) except + nogil 
        MzMLSwathFileConsumer(libcpp_vector[ SwathMap ] known_window_boundaries, 
                              String cachedir,
                              String basename,
                              Size nr_ms1_spectra, 
                              libcpp_vector[ int ] nr_ms2_spectra) except + nogil 
