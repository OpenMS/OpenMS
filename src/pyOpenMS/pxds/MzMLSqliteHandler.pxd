from MSExperiment  cimport *
from MSSpectrum  cimport *
from Peak1D  cimport *
from ChromatogramPeak  cimport *
from MSChromatogram  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from PeakFileOptions cimport *
from IMSDataConsumer cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>" namespace "OpenMS::Internal":

    cdef cppclass MzMLSqliteHandler:

        MzMLSqliteHandler(String filename, UInt64 run_id) nogil except +
        MzMLSqliteHandler(MzMLSqliteHandler &) nogil except + # compiler

        void readExperiment(MSExperiment & exp, bool meta_only )  nogil except +
            # wrap-doc:
                #   Read an experiment into an MSExperiment structure
                #   -----
                #   :param exp: The result data structure
                #   :param meta_only: Only read the meta data
  
        void readSpectra(libcpp_vector[MSSpectrum] & exp, libcpp_vector[int] indices, bool meta_only ) nogil except +
            # wrap-doc:
                #   Read a set of spectra (potentially restricted to a subset)
                #   -----
                #   :param exp: The result data structure
                #   :param indices: A list of indices restricting the resulting spectra only to those specified here
                #   :param meta_only: Only read the meta data

        void readChromatograms(libcpp_vector[MSChromatogram] & exp, libcpp_vector[int] indices, bool meta_only ) nogil except +
            # wrap-doc:
                #   Read a set of chromatograms (potentially restricted to a subset)
                #   -----
                #   :param exp: The result data structure
                #   :param indices: A list of indices restricting the resulting spectra only to those specified here
                #   :param meta_only: Only read the meta data
  
        Size getNrSpectra() nogil except + # wrap-doc:Returns number of spectra in the file, reutrns the number of spectra
  
        Size getNrChromatograms() nogil except + # wrap-doc:Returns the number of chromatograms in the file
  
        void setConfig(bool write_full_meta, bool use_lossy_compression, double linear_abs_mass_acc)  nogil except +
            # wrap-doc:
                #   Sets file configuration
                #   -----
                #   :param write_full_meta: Whether to write a complete mzML meta data structure into the RUN_EXTRA field (allows complete recovery of the input file)
                #   :param use_lossy_compression: Whether to use lossy compression (ms numpress)
                #   :param linear_abs_mass_acc: Accepted loss in mass accuracy (absolute m/z, in Th)
  
        libcpp_vector[size_t] getSpectraIndicesbyRT(double RT, double deltaRT, libcpp_vector[int] indices) nogil except +
            # wrap-doc:
                #   Returns spectral indices around a specific retention time
                #   -----
                #   :param RT: The retention time
                #   :param deltaRT: Tolerance window around RT (if less or equal than zero, only the first spectrum *after* RT is returned)
                #   :param indices: Spectra to consider (if empty, all spectra are considered)
                #   :returns: The indices of the spectra within RT +/- deltaRT
  
        void writeExperiment(MSExperiment exp) nogil except + # wrap-doc:Write an MSExperiment to disk
  
        void createTables() nogil except + # wrap-doc:Create data tables for a new file
  
        void writeSpectra(libcpp_vector[MSSpectrum] spectra) nogil except + # wrap-doc:Writes a set of spectra to disk
  
        void writeChromatograms(libcpp_vector[MSChromatogram] chroms) nogil except + # wrap-doc:Writes a set of chromatograms to disk
  
        void writeRunLevelInformation(MSExperiment exp, bool write_full_meta) nogil except +
            # wrap-doc:
                #   Write the run-level information for an experiment into tables
                #   -----
                #   This is a low level function, do not call this function unless you know what you are doing
                #   -----
                #   :param exp: The result data structure
                #   :param meta_only: Only read the meta data
        
        UInt64 getRunID() nogil except + # wrap-doc:Extract the `RUN` ID from the sqMass file
        
