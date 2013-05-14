from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Residue cimport *
from DefaultParamHandler cimport *
from AASequence cimport *
from MSSpectrum cimport *
from RichPeak1D cimport *

cdef extern from "<OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>" namespace "OpenMS":
    
    cdef cppclass TheoreticalSpectrumGenerator "OpenMS::TheoreticalSpectrumGenerator":
        TheoreticalSpectrumGenerator() nogil except +
        TheoreticalSpectrumGenerator(TheoreticalSpectrumGenerator) nogil except +
        void getSpectrum(MSSpectrum[RichPeak1D] &spec, AASequence &peptide, Int charge) nogil except +
        void addPeaks(MSSpectrum[RichPeak1D] &spectrum, AASequence &peptide, ResidueType res_type, Int charge) nogil except +
        void addPrecursorPeaks(MSSpectrum[RichPeak1D] &spec, AASequence &peptide, Int charge) nogil except +
        void addAbundantImmoniumIons(MSSpectrum[RichPeak1D] &spec) nogil except +

