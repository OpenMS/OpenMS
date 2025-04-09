from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from NASequence cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>" namespace "OpenMS":
    
    cdef cppclass NucleicAcidSpectrumGenerator(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        NucleicAcidSpectrumGenerator() except + nogil 
        NucleicAcidSpectrumGenerator(NucleicAcidSpectrumGenerator &) except + nogil 
        void getSpectrum(MSSpectrum &spec, NASequence &oligo, Int min_charge, Int max_charge) except + nogil  # wrap-doc:Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters. If precursor_charge is set to 0 max_charge + 1 will be used
