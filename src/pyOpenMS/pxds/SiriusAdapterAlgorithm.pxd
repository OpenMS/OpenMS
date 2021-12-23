from Types cimport *
from String cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from FeatureMapping cimport * 
from KDTreeFeatureMaps cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair

from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>" namespace "OpenMS":

    cdef cppclass SiriusAdapterAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        SiriusAdapterAlgorithm() nogil except +
        SiriusAdapterAlgorithm(SiriusAdapterAlgorithm &) nogil except + # compiler

        bool isFeatureOnly() nogil except +
        UInt getFilterByNumMassTraces() nogil except +
        double getPrecursorMzTolerance() nogil except +
        double getPrecursorRtTolerance() nogil except +
        bool precursorMzToleranceUnitIsPPM() nogil except +
        bool isNoMasstraceInfoIsotopePattern() nogil except +
        int getIsotopePatternIterations() nogil except +
        int getNumberOfSiriusCandidates() nogil except +
        int getNumberOfCSIFingerIDCandidates() nogil except +

        String determineSiriusExecutable(String & executable) nogil except +
        # wrap-doc:
                #   Checks if the provided String points to a valid SIRIUS executable, otherwise tries
                #   to select the executable from the environment
                #   -----
                #   :param executable: Path to the potential executable
                #   :returns: Path to SIRIUS executable

        void preprocessingSirius(const String& featureinfo,
                                 MSExperiment& spectra,
                                 FeatureMapping_FeatureMappingInfo& fm_info,
                                 FeatureMapping_FeatureToMs2Indices& feature_mapping) nogil except +
        # wrap-doc:
                #   Preprocessing needed for SIRIUS
                #   -----
                #   Filter number of masstraces and perform feature mapping
                #   -----
                #   :param featureinfo: Path to featureXML
                #   :param spectra: Input of MSExperiment with spectra information
                #   :param fm_info: Emtpy - stores FeatureMaps and KDTreeMaps internally 
                #   :param feature_mapping: Empty FeatureToMs2Indices

        void logFeatureSpectraNumber(const String& featureinfo,
                                     FeatureMapping_FeatureToMs2Indices& feature_mapping,
                                     MSExperiment& spectra) nogil except +
        # wrap-doc:
                #   Logs number of features and spectra used
                #   -----
                #   Prints the number of features and spectra used (OPENMS_LOG_INFO)
                #   -----
                #   :param featureinfo: Path to featureXML
                #   :param feature_mapping: FeatureToMs2Indices with feature mapping
                #   :param spectra: Input of MSExperiment with spectra information

        libcpp_vector[String] callSiriusQProcess(const String& tmp_ms_file,
                                                 const String& tmp_out_dir,
                                                 String& executable,
                                                 const String& out_csifingerid,
                                                 bool decoy_generation) nogil except +
        # wrap-doc:
                #   Call SIRIUS with QProcess
                #   -----
                #   :param tmp_ms_file: Path to temporary .ms file
                #   :param tmp_out_dir: Path to temporary output folder
                #   :param executable: Path to executable
                #   :param out_csifingerid: Path to CSI:FingerID output (can be empty)

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>" namespace "OpenMS::SiriusAdapterAlgorithm":

    cdef cppclass SiriusTemporaryFileSystemObjects "OpenMS::SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects":
        SiriusTemporaryFileSystemObjects(int debug_level) nogil except +
        SiriusTemporaryFileSystemObjects(SiriusTemporaryFileSystemObjects &) nogil except + # compiler
        
        String getTmpDir() nogil except +
        String getTmpOutDir() nogil except +
        String getTmpMsFile() nogil except + 

# wrap static method:
cdef extern from "<OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>" namespace "OpenMS::SiriusAdapterAlgorithm":

        void  sortSiriusWorkspacePathsByScanIndex(libcpp_vector[ String ]& subdirs) nogil except + # wrap-attach:SiriusAdapterAlgorithm
