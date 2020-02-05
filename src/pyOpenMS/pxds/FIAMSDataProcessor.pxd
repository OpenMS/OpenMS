from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from SavitzkyGolayFilter cimport *
from PeakPickerHiRes cimport *
from MzTab cimport *
from MzMLFile cimport *

cdef extern from "</home/svegal/OpenMS/src/openms/include/OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>" namespace "OpenMS":

    cdef cppclass FIAMSDataProcessor(DefaultParamHandler) :
        #
        # wrap-doc:
        #     ADD PYTHON DOCUMENTATION HERE
        #
        # wrap-inherits:
        #  DefaultParamHandler
        FIAMSDataProcessor() nogil except +
        FIAMSDataProcessor(FIAMSDataProcessor) nogil except +
        # NAMESPACE # void run(MSExperiment & experiment, float & n_seconds, OpenMS::MzTab & output) nogil except +
        void cutForTime(MSExperiment & experiment, float & n_seconds, libcpp_vector[ MSSpectrum ] & output) nogil except +
        # NAMESPACE # MSSpectrum mergeAlongTime(libcpp_vector[ OpenMS::MSSpectrum ] & input_) nogil except +
        MSSpectrum extractPeaks(MSSpectrum & input_) nogil except +
        FeatureMap convertToFeatureMap(MSSpectrum & input_) nogil except +
        MSSpectrum trackNoise(MSSpectrum & input_) nogil except +
        # NAMESPACE # void runAccurateMassSearch(FeatureMap & input_, OpenMS::MzTab & output) nogil except +
        libcpp_vector[ float ] getMZs() nogil except +
        libcpp_vector[ float ] getBinSizes() nogil except +

