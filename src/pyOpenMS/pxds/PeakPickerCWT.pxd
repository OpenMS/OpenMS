from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>" namespace "OpenMS":

    cdef cppclass PeakPickerCWT(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger
        # wrap-doc:
        #   This class implements a peak picking algorithm using wavelet techniques
        #   -----
        #   The algorithm is described in detail in Lange et al. (2006) Proc. PSB-06
        #   -----
        #   This peak picking algorithm uses the continuous wavelet transform of a profile data signal to detect mass peaks.
        #   Afterwards a given asymmetric peak function is fitted to the profile data and important peak parameters (e.g. fwhm)
        #   are extracted
        #   In an optional step these parameters can be optimized using a non-linear optimization method
        #   -----
        #   The peak parameters are stored in the meta data arrays of the spectra (see MSSpectrum) in this order:
        #   - rValue
        #   - area
        #   - fwhm
        #   - leftWidth
        #   - rightWidth
        #   - peakShape
        #   - SignalToNoise

        PeakPickerCWT() nogil except +
        PeakPickerCWT(PeakPickerCWT &) nogil except + # compiler

        void pick(MSSpectrum & input, MSSpectrum & output) nogil except +
            # wrap-doc:
                #   Applies the peak picking algorithm to a single spectrum
                #   -----
                #   Picks the peaks in the input spectrum and writes the resulting peaks to the output container

        void pickExperiment(MSExperiment & input, MSExperiment & output) nogil except +
             # wrap-doc:
                #   Picks the peaks in an MSExperiment
                #   -----
                #   Picks the peaks successive in every scan in the spectrum range. The detected peaks are stored in the output MSExperiment

        double estimatePeakWidth(MSExperiment & input) nogil except +
            # wrap-doc:
                #   Estimates average peak width that can then be used for peak picking
                #   -----
                #   The spectra with the highest TICs are used to estimate an average peak width that
                #   can be used as the peak_width parameter for picking the complete data set.
                #   Typically, the number of peaks increases with decreasing peak width until a plateau
                #   is reached. The beginning of this plateau is our estimate for the peak width.
                #   This estimate is averaged over several spectra
