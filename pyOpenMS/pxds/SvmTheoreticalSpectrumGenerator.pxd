from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from SimTypes cimport *
from SVMWrapper cimport *
from AASequence cimport *
from IonType cimport *

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>" namespace "OpenMS":
    
    cdef cppclass SvmTheoreticalSpectrumGenerator "OpenMS::SvmTheoreticalSpectrumGenerator":
        SvmTheoreticalSpectrumGenerator() nogil except +
        SvmTheoreticalSpectrumGenerator(SvmTheoreticalSpectrumGenerator) nogil except +
        # TODO pointer to gsl_rng
        # void simulate(MSSpectrum[RichPeak1D] &spectrum, AASequence &peptide, gsl_rng *rng, Size precursor_charge) nogil except +
        void load() nogil except +
        libcpp_vector[ IonType ]  getIonTypes() nogil except +

