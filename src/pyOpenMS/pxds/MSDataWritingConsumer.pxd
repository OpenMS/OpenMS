from MSSpectrum cimport *
from MSChromatogram cimport *
from ExperimentalSettings cimport *
from DataProcessing cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>" namespace "OpenMS":

    cdef cppclass PlainMSDataWritingConsumer:

        PlainMSDataWritingConsumer(String filename) nogil except +
        # copy constructor of 'PlainMSDataWritingConsumer' is implicitly deleted because base class 'OpenMS::MSDataWritingConsumer' has a deleted copy constructor (see XMLHandler)
        PlainMSDataWritingConsumer(PlainMSDataWritingConsumer &) nogil except + # wrap-ignore

        void consumeSpectrum(MSSpectrum & s) nogil except +
        void consumeChromatogram(MSChromatogram & c) nogil except +

        void setExperimentalSettings(ExperimentalSettings& exp) nogil except +
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +

        void addDataProcessing(DataProcessing d) nogil except +
        Size getNrSpectraWritten()  nogil except +
        Size getNrChromatogramsWritten() nogil except +

        void setOptions(PeakFileOptions opt) nogil except +
        PeakFileOptions getOptions() nogil except +

    cdef cppclass NoopMSDataWritingConsumer:

        NoopMSDataWritingConsumer(String filename) nogil except +
        # copy constructor of 'NoopMSDataWritingConsumer' is implicitly deleted because base class 'OpenMS::MSDataWritingConsumer' has a deleted copy constructor (see XMLHandler)
        NoopMSDataWritingConsumer(NoopMSDataWritingConsumer &) nogil except + # wrap-ignore

        void consumeSpectrum(MSSpectrum & s) nogil except +
        void consumeChromatogram(MSChromatogram & c) nogil except +

        void setExperimentalSettings(ExperimentalSettings& exp) nogil except +
        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +

        void addDataProcessing(DataProcessing d) nogil except +
        Size getNrSpectraWritten()  nogil except +
        Size getNrChromatogramsWritten() nogil except +

