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
        #   DefaultParamHandler

        TheoreticalSpectrumGenerator() except + nogil 
        TheoreticalSpectrumGenerator(TheoreticalSpectrumGenerator &) except + nogil 
        void getSpectrum(MSSpectrum &spec, AASequence &peptide, Int min_charge, Int max_charge) except + nogil  # wrap-doc:Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters. If precursor_charge is set to 0 max_charge + 1 will be used
