from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from SavitzkyGolayFilter cimport *
from PeakPickerHiRes cimport *
from MzTab cimport *
from MzMLFile cimport *
from FeatureMapping cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>" namespace "OpenMS":

    cdef cppclass FIAMSDataProcessor(DefaultParamHandler) :
        #
        # wrap-doc:
        #     ADD PYTHON DOCUMENTATION HERE
        #
        # wrap-inherits:
        #  DefaultParamHandler
        FIAMSDataProcessor() nogil except +
        FIAMSDataProcessor(FIAMSDataProcessor) nogil except +
        bool run(MSExperiment & experiment, float & n_seconds, MzTab & output, bool load_cached_spectrum) nogil except +
        # void cutForTime(MSExperiment & experiment, float & n_seconds, libcpp_vector[ MSSpectrum ] & output) nogil except +
        # NAMESPACE # MSSpectrum mergeAlongTime(libcpp_vector[ OpenMS::MSSpectrum ] & input_) nogil except +
        MSSpectrum extractPeaks(MSSpectrum & input_) nogil except +
        FeatureMap convertToFeatureMap(MSSpectrum & input_) nogil except +
        MSSpectrum trackNoise(MSSpectrum & input_) nogil except +
        # NAMESPACE # void runAccurateMassSearch(FeatureMap & input_, OpenMS::MzTab & output) nogil except +
        # libcpp_vector[ float ] getMZs() nogil except +
        # libcpp_vector[ float ] getBinSizes() nogil except +

