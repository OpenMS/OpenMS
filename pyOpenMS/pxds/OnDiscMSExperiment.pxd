from Types cimport *
from MSExperiment cimport *
from ExperimentalSettings cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/KERNEL/OnDiscMSExperiment.h>" namespace "OpenMS":

    cdef cppclass OnDiscMSExperiment[PeakT, ChromoPeakT](ExperimentalSettings):
        # wrap-instances:
        #   OnDiscMSExperiment := OnDiscMSExperiment[Peak1D, ChromatogramPeak]

        OnDiscMSExperiment() nogil except +
        bool openFile(String filename) nogil except +
        Size getNrSpectra() nogil except +
        Size getNrChromatograms() nogil except +
        # shared_ptr[ExperimentalSettings] getExperimentalSettings() nogil except +
        MSSpectrum[PeakT] getSpectrum(Size id) nogil except +
        MSChromatogram[ChromoPeakT] getChromatogram(Size id) nogil except +
