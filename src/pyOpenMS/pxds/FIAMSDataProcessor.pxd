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

        FIAMSDataProcessor() nogil except + # wrap-doc:Data processing for FIA-MS data
        FIAMSDataProcessor(FIAMSDataProcessor &) nogil except +

        bool run(MSExperiment & experiment, float & n_seconds, MzTab & output, bool load_cached_spectrum) nogil except +
            # wrap-doc:
                #   Run the full analysis for the experiment for the given time interval
                #   -----
                #   The workflow steps are:
                #   - the time axis of the experiment is cut to the interval from 0 to n_seconds
                #   - the spectra are summed into one along the time axis with the bin size determined by mz and instrument resolution
                #   - data is smoothed by applying the Savitzky-Golay filter
                #   - peaks are picked
                #   - the accurate mass search for all the picked peaks is performed
                #
                #   The intermediate summed spectra and picked peaks can be saved to the filesystem. 
                #   Also, the results of the accurate mass search and the signal-to-noise information 
                #   of the resulting spectrum is saved.
                #   -----
                #   :param experiment: Input MSExperiment
                #   :param n_seconds: Input number of seconds
                #   :param load_cached_spectrum: Load the cached picked spectrum if exists
                #   :param output: Output of the accurate mass search results
                #   :returns: A boolean indicating if the picked spectrum was loaded from the cached file
        # void cutForTime(MSExperiment & experiment, float & n_seconds, libcpp_vector[ MSSpectrum ] & output) nogil except +
        # NAMESPACE # MSSpectrum mergeAlongTime(libcpp_vector[ OpenMS::MSSpectrum ] & input_) nogil except +
        MSSpectrum extractPeaks(MSSpectrum & input_) nogil except +
            # wrap-doc:
                #   Pick peaks from the summed spectrum
                #   -----
                #   :param input: Input vector of spectra
                #   :returns: A spectrum with picked peaks
        FeatureMap convertToFeatureMap(MSSpectrum & input_) nogil except +
            # wrap-doc:
                #   Convert a spectrum to a feature map with the corresponding polarity
                #   -----
                #   Applies `SavitzkyGolayFilter` and `PeakPickerHiRes`
                #   -----
                #   :param input: Input a picked spectrum
                #   :returns: A feature map with the peaks converted to features and polarity from the parameters
        MSSpectrum trackNoise(MSSpectrum & input_) nogil except +
            # wrap-doc:
                #   Estimate noise for each peak
                #   -----
                #   Uses `SignalToNoiseEstimatorMedianRapid`
                #   -----
                #   :param input: Input a picked spectrum
                #   :returns: A spectrum object storing logSN information
        # NAMESPACE # void runAccurateMassSearch(FeatureMap & input_, OpenMS::MzTab & output) nogil except +
        # libcpp_vector[ float ] getMZs() nogil except +
        # libcpp_vector[ float ] getBinSizes() nogil except +
