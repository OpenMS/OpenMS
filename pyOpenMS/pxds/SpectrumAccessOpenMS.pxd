from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>" namespace "OpenMS":

  cdef cppclass SpectrumAccessOpenMS:
        SpectrumAccessOpenMS(SpectrumAccessOpenMS) # wrap-ignore
        SpectrumAccessOpenMS(shared_ptr[ MSExperiment[Peak1D, ChromatogramPeak] ] & ms_experiment)

        shared_ptr[Spectrum] getSpectrumById(int id)  #wrap-ignore
        libcpp_vector[size_t] getSpectraByRT(double RT, double deltaRT)
        size_t getNrSpectra()

        shared_ptr[Chromatogram] getChromatogramById(int id)  #wrap-ignore
        size_t getNrChromatograms()
        libcpp_string getChromatogramNativeID(int id_)

