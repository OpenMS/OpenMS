from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from SimTypes cimport *
from SvmTheoreticalSpectrumGenerator cimport *

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>" namespace "OpenMS":
    
    cdef cppclass SvmTheoreticalSpectrumGeneratorSet "OpenMS::SvmTheoreticalSpectrumGeneratorSet":
        SvmTheoreticalSpectrumGeneratorSet() nogil except +
        SvmTheoreticalSpectrumGeneratorSet(SvmTheoreticalSpectrumGeneratorSet &) nogil except +

        # void simulate(MSSpectrum &spectrum, AASequence &peptide, boost::random::mt19937_64& rng, Size precursor_charge) nogil except +
        void load(String) nogil except +
        void getSupportedCharges(libcpp_set[ size_t ] &charges) nogil except +
        SvmTheoreticalSpectrumGenerator  getSvmModel(Size) nogil except +
