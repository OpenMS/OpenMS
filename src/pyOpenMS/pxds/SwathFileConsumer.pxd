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
        #wrap-ignore
        #no-pxd-import

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
