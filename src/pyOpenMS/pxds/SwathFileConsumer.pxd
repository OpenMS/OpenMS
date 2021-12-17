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

        FullSwathFileConsumer() nogil except + #wrap-ignore
        FullSwathFileConsumer(FullSwathFileConsumer &) nogil except + # compiler
        FullSwathFileConsumer(libcpp_vector[SwathMap] swath_boundaries) nogil except +

        void setExpectedSize(Size s, Size c) nogil except +
        void setExperimentalSettings(ExperimentalSettings exp) nogil except +

        void retrieveSwathMaps(libcpp_vector[SwathMap]& maps) nogil except +

        void consumeSpectrum(MSSpectrum & s) nogil except +
        void consumeChromatogram(MSChromatogram & c) nogil except +

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>" namespace "OpenMS":

    cdef cppclass RegularSwathFileConsumer(FullSwathFileConsumer):
        # wrap-inherits:
        #    FullSwathFileConsumer

        RegularSwathFileConsumer() nogil except +
        RegularSwathFileConsumer(RegularSwathFileConsumer &) nogil except + # compiler


cdef extern from "<OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>" namespace "OpenMS":

    cdef cppclass CachedSwathFileConsumer(FullSwathFileConsumer):
        # wrap-inherits:
        #    FullSwathFileConsumer

        CachedSwathFileConsumer() nogil except + #wrap-ignore
        CachedSwathFileConsumer(CachedSwathFileConsumer &) nogil except + # compiler
        CachedSwathFileConsumer(String cachedir, String basename, 
                                Size nr_ms1_spectra, libcpp_vector[int] nr_ms2_spectra)

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MzMLSwathFileConsumer(FullSwathFileConsumer) :
        # wrap-inherits:
        #  FullSwathFileConsumer

        MzMLSwathFileConsumer() nogil except + # wrap-ignore
        MzMLSwathFileConsumer(MzMLSwathFileConsumer) nogil except + # compiler
        MzMLSwathFileConsumer(String cachedir,
                              String basename, 
                              Size nr_ms1_spectra, 
                              libcpp_vector[ int ] nr_ms2_spectra) nogil except +
        MzMLSwathFileConsumer(libcpp_vector[ SwathMap ] known_window_boundaries, 
                              String cachedir,
                              String basename,
                              Size nr_ms1_spectra, 
                              libcpp_vector[ int ] nr_ms2_spectra) nogil except +
