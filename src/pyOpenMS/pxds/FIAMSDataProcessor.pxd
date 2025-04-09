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
        #    ADD PYTHON DOCUMENTATION HERE
        #
        # wrap-inherits:
        #  DefaultParamHandler

        FIAMSDataProcessor() except + nogil  # wrap-doc:Data processing for FIA-MS data
        FIAMSDataProcessor(FIAMSDataProcessor &) except + nogil 

        bool run(MSExperiment & experiment, float & n_seconds, MzTab & output, bool load_cached_spectrum) except + nogil 
            # wrap-doc:
                #  Run the full analysis for the experiment for the given time interval\n
                #  
                #  The workflow steps are:
                #  - the time axis of the experiment is cut to the interval from 0 to n_seconds
                #  - the spectra are summed into one along the time axis with the bin size determined by mz and instrument resolution
                #  - data is smoothed by applying the Savitzky-Golay filter
                #  - peaks are picked
                #  - the accurate mass search for all the picked peaks is performed
                #  
                #  The intermediate summed spectra and picked peaks can be saved to the filesystem. 
                #  Also, the results of the accurate mass search and the signal-to-noise information 
                #  of the resulting spectrum is saved.
                #  
                #  
                #  :param experiment: Input MSExperiment
                #  :param n_seconds: Input number of seconds
                #  :param load_cached_spectrum: Load the cached picked spectrum if exists
                #  :param output: Output of the accurate mass search results
                #  :return: A boolean indicating if the picked spectrum was loaded from the cached file

        # void cutForTime(MSExperiment & experiment, float & n_seconds, libcpp_vector[ MSSpectrum ] & output) except + nogil 
        # NAMESPACE # MSSpectrum mergeAlongTime(libcpp_vector[ OpenMS::MSSpectrum ] & input_) except + nogil 
        MSSpectrum extractPeaks(MSSpectrum & input_) except + nogil 
            # wrap-doc:
                #  Pick peaks from the summed spectrum
                #  
                #  
                #  :param input: Input vector of spectra
                #  :return: A spectrum with picked peaks

        FeatureMap convertToFeatureMap(MSSpectrum & input_) except + nogil 
            # wrap-doc:
                #  Convert a spectrum to a feature map with the corresponding polarity\n
                #  
                #  Applies `SavitzkyGolayFilter` and `PeakPickerHiRes`
                #  
                #  
                #  :param input: Input a picked spectrum
                #  :return: A feature map with the peaks converted to features and polarity from the parameters

        MSSpectrum trackNoise(MSSpectrum & input_) except + nogil 
            # wrap-doc:
                #  Estimate noise for each peak\n
                #  
                #  Uses `SignalToNoiseEstimatorMedianRapid`
                #  
                #  
                #  :param input: Input a picked spectrum
                #  :return: A spectrum object storing logSN information
        # NAMESPACE # void runAccurateMassSearch(FeatureMap & input_, OpenMS::MzTab & output) except + nogil 
        # libcpp_vector[ float ] getMZs() except + nogil 
        # libcpp_vector[ float ] getBinSizes() except + nogil 
