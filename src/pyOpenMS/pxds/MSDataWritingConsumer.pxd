from MSSpectrum cimport *
from MSChromatogram cimport *
from ExperimentalSettings cimport *
from DataProcessing cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>" namespace "OpenMS":

    cdef cppclass PlainMSDataWritingConsumer:

        PlainMSDataWritingConsumer(String filename) nogil except +
        PlainMSDataWritingConsumer(PlainMSDataWritingConsumer) nogil except + #wrap-ignore

        void consumeSpectrum(MSSpectrum[Peak1D] & s) nogil except + 
        void consumeChromatogram(MSChromatogram[ChromatogramPeak] & c) nogil except + 

        void setExperimentalSettings(ExperimentalSettings& exp) nogil except +
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +

        void addDataProcessing(DataProcessing d) nogil except +
        Size getNrSpectraWritten()  nogil except +
        Size getNrChromatogramsWritten() nogil except +


    cdef cppclass NoopMSDataWritingConsumer:

        NoopMSDataWritingConsumer(String filename) nogil except +
        NoopMSDataWritingConsumer(NoopMSDataWritingConsumer) nogil except + #wrap-ignore

        void consumeSpectrum(MSSpectrum[Peak1D] & s) nogil except + 
        void consumeChromatogram(MSChromatogram[ChromatogramPeak] & c) nogil except + 

        void setExperimentalSettings(ExperimentalSettings& exp) nogil except +
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +

        void addDataProcessing(DataProcessing d) nogil except +
        Size getNrSpectraWritten()  nogil except +
        Size getNrChromatogramsWritten() nogil except +

