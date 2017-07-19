from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Residue cimport *
from DefaultParamHandler cimport *
from AASequence cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>" namespace "OpenMS":
    
    cdef cppclass TheoreticalSpectrumGenerator(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        TheoreticalSpectrumGenerator() nogil except +
        TheoreticalSpectrumGenerator(TheoreticalSpectrumGenerator) nogil except +
        void getSpectrum(MSSpectrum &spec, AASequence &peptide, Int min_charge, Int max_charge) nogil except +
