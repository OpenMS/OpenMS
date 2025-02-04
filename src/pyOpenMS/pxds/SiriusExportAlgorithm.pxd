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

        void preprocessing(const String& featureXML_path,
                                 MSExperiment& spectra,
                                 FeatureMapping_FeatureMappingInfo& feature_mapping_info,
                                 FeatureMapping_FeatureToMs2Indices& feature_ms2_indices) except + nogil 
        # wrap-doc:
                #  Preprocessing needed for SIRIUS
                #  
                #  Filter number of masstraces and perform feature mapping
                #  
                #  :param featureXML_path: Path to featureXML
                #  :param spectra: Input of MSExperiment with spectra information
                #  :param feature_mapping_info: Emtpy - stores FeatureMaps and KDTreeMaps internally 
                #  :param feature_ms2_indices: Empty FeatureToMs2Indices

        void logFeatureSpectraNumber(const String& featureXML_path,
                                     FeatureMapping_FeatureToMs2Indices& feature_ms2_indices,
                                     MSExperiment& spectra) except + nogil 
        # wrap-doc:
                #  Logs number of features and spectra used
                #  
                #  Prints the number of features and spectra used (OPENMS_LOG_INFO)
                #  
                #  :param featureXML_path: Path to featureXML
                #  :param feature_ms2_indices: FeatureToMs2Indices with feature mapping
                #  :param spectra: Input of MSExperiment with spectra information

        void run(const StringList& mzML_files,
                 const StringList& featureXML_files,
                 const String& out_ms,
                 const String& out_compoundinfo) except + nogil

        # wrap-doc:
                #  Runs SiriusExport with mzML and featureXML (optional) files as input.
                #  
                #  Generates a SIRIUS .ms file and compound info table (optional).
                #  
                #  :param mzML_files: List with paths to mzML files
                #  :param featureXML_files: List with paths to featureXML files
                #  :param out_ms: Output file name for SIRIUS .ms file
                #  :param out_compoundinfo: Output file name for tsv file with compound info