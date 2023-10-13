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
        #   DefaultParamHandler

        SiriusAdapterAlgorithm() except + nogil 
        SiriusAdapterAlgorithm(SiriusAdapterAlgorithm &) except + nogil  # compiler

        bool isFeatureOnly() except + nogil 
        UInt getFilterByNumMassTraces() except + nogil 
        double getPrecursorMzTolerance() except + nogil 
        double getPrecursorRtTolerance() except + nogil 
        bool precursorMzToleranceUnitIsPPM() except + nogil 
        bool isNoMasstraceInfoIsotopePattern() except + nogil 
        int getIsotopePatternIterations() except + nogil 
        int getNumberOfSiriusCandidates() except + nogil 

        String determineSiriusExecutable(String & executable) except + nogil 
        # wrap-doc:
                #  Checks if the provided String points to a valid SIRIUS executable, otherwise tries
                #  to select the executable from the environment
                #  
                #  :param executable: Path to the potential executable
                #  :returns: Path to SIRIUS executable

        void preprocessingSirius(const String& featureinfo,
                                 MSExperiment& spectra,
                                 FeatureMapping_FeatureMappingInfo& fm_info,
                                 FeatureMapping_FeatureToMs2Indices& feature_mapping) except + nogil 
        # wrap-doc:
                #  Preprocessing needed for SIRIUS
                #  
                #  Filter number of masstraces and perform feature mapping
                #  
                #  :param featureinfo: Path to featureXML
                #  :param spectra: Input of MSExperiment with spectra information
                #  :param fm_info: Emtpy - stores FeatureMaps and KDTreeMaps internally 
                #  :param feature_mapping: Empty FeatureToMs2Indices

        void logFeatureSpectraNumber(const String& featureinfo,
                                     FeatureMapping_FeatureToMs2Indices& feature_mapping,
                                     MSExperiment& spectra) except + nogil 
        # wrap-doc:
                #  Logs number of features and spectra used
                #  
                #  Prints the number of features and spectra used (OPENMS_LOG_INFO)
                #  
                #  :param featureinfo: Path to featureXML
                #  :param feature_mapping: FeatureToMs2Indices with feature mapping
                #  :param spectra: Input of MSExperiment with spectra information

        void logInSiriusAccount(String& executable,
                                const String& email,
                                const String& password) except + nogil 
        # wrap-doc:
                #  Log in to SIRIUS using your personal account
                #  
                #  :param executable: Path to executable.
                #  :param email: User account E-Mail.
                #  :param password: User account password.
     
        libcpp_vector[String] callSiriusQProcess(const String& tmp_ms_file,
                                                 const String& tmp_out_dir,
                                                 String& executable,
                                                 const String& out_csifingerid,
                                                 bool decoy_generation) except + nogil 
        # wrap-doc:
                #  Call SIRIUS with QProcess
                #  
                #  :param tmp_ms_file: Path to temporary .ms file
                #  :param tmp_out_dir: Path to temporary output folder
                #  :param executable: Path to executable
                #  :param out_csifingerid: Path to CSI:FingerID output (can be empty)

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>" namespace "OpenMS::SiriusAdapterAlgorithm":

    cdef cppclass SiriusTemporaryFileSystemObjects "OpenMS::SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects":
        SiriusTemporaryFileSystemObjects(int debug_level) except + nogil 
        SiriusTemporaryFileSystemObjects(SiriusTemporaryFileSystemObjects &) except + nogil  # compiler
        
        String getTmpDir() except + nogil 
        String getTmpOutDir() except + nogil 
        String getTmpMsFile() except + nogil  

# wrap static method:
cdef extern from "<OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>" namespace "OpenMS::SiriusAdapterAlgorithm":

        void  sortSiriusWorkspacePathsByScanIndex(libcpp_vector[ String ]& subdirs) except + nogil  # wrap-attach:SiriusAdapterAlgorithm
