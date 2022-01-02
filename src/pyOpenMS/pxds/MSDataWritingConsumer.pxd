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
            # wrap-doc:
                #   Set experimental settings for the whole file
                #   -----
                #   :param exp: Experimental settings to be used for this file (from this
                #   and the first spectrum/chromatogram, the class will deduce most of
                #   the header of the mzML file)

        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +
            # wrap-doc:
                #   Set expected size of spectra and chromatograms to be written
                #   -----
                #   These numbers will be written in the spectrumList and chromatogramList
                #   tag in the mzML file. Therefore, these will contain wrong numbers if
                #   the expected size is not set correctly
                #   -----
                #   :param expectedSpectra: Number of spectra expected
                #   :param expectedChromatograms: Number of chromatograms expected

        void addDataProcessing(DataProcessing d) nogil except +
            # wrap-doc:
                #   Optionally add a data processing method to each chromatogram and spectrum
                #   -----
                #   The provided DataProcessing object will be added to each chromatogram
                #   and spectrum written to to the mzML file
                #   -----
                #   :param d: The DataProcessing object to be added

        Size getNrSpectraWritten()  nogil except + # wrap-doc:Returns the number of spectra written
        Size getNrChromatogramsWritten() nogil except + # wrap-doc:Returns the number of chromatograms written

        void setOptions(PeakFileOptions opt) nogil except +
        PeakFileOptions getOptions() nogil except +

    cdef cppclass NoopMSDataWritingConsumer:
        # wrap-doc:
                #   Consumer class that perform no operation
                #   -----
                #   This is sometimes necessary to fulfill the requirement of passing an
                #   valid MSDataWritingConsumer object or pointer but no operation is
                #   required

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
