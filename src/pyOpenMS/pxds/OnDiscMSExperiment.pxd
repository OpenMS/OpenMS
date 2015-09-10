from Types cimport *
from MSExperiment cimport *
from ExperimentalSettings cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from InterfaceDataStructures cimport *

cdef extern from "<OpenMS/KERNEL/OnDiscMSExperiment.h>" namespace "OpenMS":

    cdef cppclass OnDiscMSExperiment[PeakT, ChromoPeakT](ExperimentalSettings):
        # wrap-instances:
        #   OnDiscMSExperiment := OnDiscMSExperiment[Peak1D, ChromatogramPeak]

        OnDiscMSExperiment() nogil except +
        OnDiscMSExperiment(OnDiscMSExperiment &) nogil except +

        bool openFile(String filename) nogil except +
        Size getNrSpectra() nogil except +
        Size getNrChromatograms() nogil except +

        # TODO const shared ptr
        # shared_ptr[const ExperimentalSettings] getExperimentalSettings() nogil except + # wrap-ignore

        MSSpectrum[PeakT] getSpectrum(Size id) nogil except +
        MSChromatogram[ChromoPeakT] getChromatogram(Size id) nogil except +

        # TODO decide for 1.12 whether to include those ... 
        shared_ptr[Spectrum] getSpectrumById(int id_) nogil except +
        shared_ptr[Chromatogram] getChromatogramById(int id_) nogil except +


