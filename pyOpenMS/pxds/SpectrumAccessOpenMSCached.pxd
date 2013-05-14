from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h>" namespace "OpenMS":

  cdef cppclass SpectrumAccessOpenMSCached:
        SpectrumAccessOpenMSCached(String filename)
        # SpectrumAccessOpenMSCached(SpectrumAccessOpenMSCached)

        shared_ptr[Spectrum] getSpectrumById(int id) nogil except + # wrap-ignore
        libcpp_vector[size_t] getSpectraByRT(double RT, double deltaRT) nogil except +
        size_t getNrSpectra() nogil except +

        shared_ptr[Chromatogram] getChromatogramById(int id) nogil except + # wrap-ignore
        size_t getNrChromatograms() nogil except +

        libcpp_string getChromatogramNativeID(int id_) nogil except +

#  -- TODO missing function in PXD:  # NAMESPACE # OpenSwath::SpectrumMeta getSpectrumMetaById(int id_)
#  -- TODO missing function in PXD:  SpectrumSettings getSpectraMetaInfo(int id_)
#  -- TODO missing function in PXD:  ChromatogramSettings getChromatogramMetaInfo(int id_)
#  -- TODO missing function in PXD:  libcpp_string getChromatogramNativeID(int id_)
