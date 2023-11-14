from Types cimport *
from String cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from FeatureMapping cimport * 
from KDTreeFeatureMaps cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair

from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h>" namespace "OpenMS":

    cdef cppclass SiriusExportAlgorithm(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        SiriusExportAlgorithm() except + nogil 
        SiriusExportAlgorithm(SiriusExportAlgorithm &) except + nogil  # compiler

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


# wrap static method:
cdef extern from "<OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h>" namespace "OpenMS::SiriusExportAlgorithm":

        void  sortSiriusWorkspacePathsByScanIndex(libcpp_vector[ String ]& subdirs) except + nogil  # wrap-attach:SiriusExportAlgorithm
