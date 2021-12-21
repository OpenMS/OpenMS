from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>" namespace "OpenMS":

    cdef cppclass PeakPickerHiRes(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        PeakPickerHiRes() nogil except +
        PeakPickerHiRes(PeakPickerHiRes &) nogil except + # compiler

        void pick(MSSpectrum & input,
                  MSSpectrum & output
                 ) nogil except +
        void pick(MSChromatogram & input,
                  MSChromatogram & output
                 ) nogil except +

        void pickExperiment(MSExperiment & input,
                            MSExperiment & output,
                            bool check_spectrum_type
                           ) nogil except +
            # wrap-doc:
                #   Applies the peak-picking algorithm to a map (MSExperiment). This
                #   method picks peaks for each scan in the map consecutively. The resulting
                #   picked peaks are written to the output map
                #   -----
                #   :param input: Input map in profile mode
                #   :param output: Output map with picked peaks
                #   :param check_spectrum_type: If set, checks spectrum type and throws an exception if a centroided spectrum is passed 

        void pickExperiment(MSExperiment & input,
                            MSExperiment & output,
                            libcpp_vector[libcpp_vector[PeakBoundary] ]& boundaries_spec,
                            libcpp_vector[libcpp_vector[PeakBoundary] ]& boundaries_chrom,
                            bool check_spectrum_type
                           ) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>" namespace "OpenMS::PeakPickerHiRes":
    
    cdef cppclass PeakBoundary "OpenMS::PeakPickerHiRes::PeakBoundary":
        PeakBoundary() nogil except + # compiler
        PeakBoundary(PeakBoundary &) nogil except + # compiler
        double mz_min
        double mz_max
