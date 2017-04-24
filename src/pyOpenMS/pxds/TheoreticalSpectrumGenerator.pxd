from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Residue cimport *
from DefaultParamHandler cimport *
from AASequence cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>" namespace "OpenMS":
    
    cdef cppclass TheoreticalSpectrumGenerator(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        TheoreticalSpectrumGenerator() nogil except +
        TheoreticalSpectrumGenerator(TheoreticalSpectrumGenerator) nogil except +
        void getSpectrum(MSSpectrum[Peak1D] &spec, AASequence &peptide, Int min_charge, Int max_charge) nogil except +

        # these functions are now protected, equal functionality is preserved by using getSpectrum and adapting Params
        # void addPeaks(MSSpectrum[Peak1D] &spectrum, AASequence &peptide, ResidueType res_type, Int charge) nogil except +
        # void addPrecursorPeaks(MSSpectrum[Peak1D] &spec, AASequence &peptide, Int charge) nogil except +
        # void addAbundantImmoniumIons(MSSpectrum[RPeak1D] &spec, AASequence &peptide) nogil except +

