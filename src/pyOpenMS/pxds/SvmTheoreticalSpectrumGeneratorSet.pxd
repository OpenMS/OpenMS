from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from SimTypes cimport *
from SvmTheoreticalSpectrumGenerator cimport *

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>" namespace "OpenMS":
    
    cdef cppclass SvmTheoreticalSpectrumGeneratorSet "OpenMS::SvmTheoreticalSpectrumGeneratorSet":
        SvmTheoreticalSpectrumGeneratorSet() except + nogil 
        SvmTheoreticalSpectrumGeneratorSet(SvmTheoreticalSpectrumGeneratorSet &) except + nogil 

        # void simulate(MSSpectrum &spectrum, AASequence &peptide, boost::random::mt19937_64& rng, Size precursor_charge) except + nogil 
        void load(String) except + nogil 
        void getSupportedCharges(libcpp_set[ size_t ] &charges) except + nogil 
        SvmTheoreticalSpectrumGenerator  getSvmModel(Size) except + nogil 
