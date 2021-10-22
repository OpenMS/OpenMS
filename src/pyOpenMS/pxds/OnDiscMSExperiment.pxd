from Types cimport *
from MSExperiment cimport *
from ExperimentalSettings cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from InterfaceDataStructures cimport *

cdef extern from "<OpenMS/KERNEL/OnDiscMSExperiment.h>" namespace "OpenMS":

    cdef cppclass OnDiscMSExperiment(ExperimentalSettings):
        # wrap-doc:
        #   Representation of a mass spectrometry experiment on disk.
        OnDiscMSExperiment() nogil except +
        OnDiscMSExperiment(OnDiscMSExperiment &) nogil except +

        bool openFile(String filename) nogil except +
        bool openFile(String filename, bool skipLoadingMetaData) nogil except +
            # wrap-doc:
                #   Open a specific file on disk
                #   -----
                #   This tries to read the indexed mzML by parsing the index and then reading the meta information into memory
                #   -----
                #   returns: Whether the parsing of the file was successful (if false, the file most likely was not an indexed mzML file)

        Size getNrSpectra() nogil except + # wrap-doc:Returns the total number of spectra available
        Size getNrChromatograms() nogil except + # wrap-doc:Returns the total number of chromatograms available

        # COMMENT: only retrieves experiment meta data (no actual data in spectra/chromatograms)
        # COMMENT: useful for filtering by attributes to then retrieve data
        shared_ptr[const ExperimentalSettings] getExperimentalSettings() nogil except + # wrap-doc:Returns the meta information of this experiment (const access)
        shared_ptr[MSExperiment] getMetaData() nogil except + # wrap-doc:Returns the meta information of this experiment

        MSSpectrum getSpectrum(Size id) nogil except +
            # wrap-doc:
                #   Returns a single spectrum
                #   -----
                #   :param id: The index of the spectrum

        MSSpectrum getSpectrumByNativeId(String id) nogil except +
            # wrap-doc:
                #   Returns a single spectrum
                #   -----
                #   :param id: The native identifier of the spectrum

        MSChromatogram getChromatogram(Size id) nogil except +
            # wrap-doc:
                #   Returns a single chromatogram
                #   -----
                #   :param id: The index of the chromatogram
                
        MSChromatogram getChromatogramByNativeId(String id) nogil except +
            # wrap-doc:
                #   Returns a single chromatogram
                #   -----
                #   :param id: The native identifier of the chromatogram

        # TODO decide for 1.12 whether to include those ... 
        shared_ptr[Spectrum] getSpectrumById(int id_) nogil except + # wrap-doc:Returns a single spectrum
        shared_ptr[Chromatogram] getChromatogramById(int id_) nogil except + # wrap-doc:Returns a single chromatogram

        void setSkipXMLChecks(bool skip) nogil except + # wrap-doc:Sets whether to skip some XML checks and be fast instead

