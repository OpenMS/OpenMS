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

cdef extern from "<OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>" namespace "OpenMS::Internal":

    cdef cppclass MzMLSqliteHandler:

        # MzMLSqliteHandler() nogil except +

        MzMLSqliteHandler(String filename) nogil except +
        MzMLSqliteHandler(MzMLSqliteHandler h) nogil except +

        void readExperiment(MSExperiment & exp, bool meta_only )  nogil except +
  
        void readSpectra(libcpp_vector[MSSpectrum] & exp, libcpp_vector[int] indices, bool meta_only ) nogil except +

        void readChromatograms(libcpp_vector[MSChromatogram] & exp, libcpp_vector[int] indices, bool meta_only ) nogil except +
  
        Size getNrSpectra() nogil except +
  
        Size getNrChromatograms() nogil except +
  
        void setConfig(bool write_full_meta, bool use_lossy_compression, double linear_abs_mass_acc)  nogil except +
  
        libcpp_vector[size_t] getSpectraIndicesbyRT(double RT, double deltaRT, libcpp_vector[int] indices) nogil except +
  
        void writeExperiment(MSExperiment exp) nogil except +
  
        void createTables() nogil except +
  
        void writeSpectra(libcpp_vector[MSSpectrum] spectra) nogil except +
  
        void writeChromatograms(libcpp_vector[MSChromatogram] chroms) nogil except +
  
        void writeRunLevelInformation(MSExperiment exp, bool write_full_meta, int run_id) nogil except +

