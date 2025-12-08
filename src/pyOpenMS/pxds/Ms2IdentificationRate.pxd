from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/QC/Ms2IdentificationRate.h>" namespace "OpenMS":
    
    cdef cppclass Ms2IdentificationRate(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  This class is a metric for the QualityControl-ToppTool.
        #  
        #  This class computes the MS2 Identification Rate (as #identified PSMs 
        #  divided by total number of MS2 scans) given a FeatureMap and an MSExperiment.

        Ms2IdentificationRate() except + nogil 
        Ms2IdentificationRate(Ms2IdentificationRate &) except + nogil 

        void compute(FeatureMap & feature_map, MSExperiment & exp, bool assume_all_target) except + nogil 
            # wrap-doc:
            #  Computes MS2 Identification Rate with FeatureMap
            #  
            #  Stores results as a struct in a vector. Only pep-ids with target/decoy 
            #  annotation as 'target' are counted, unless assume_all_target flag is set.
            #  
            #  :param feature_map: Input FeatureMap with target/decoy annotation
            #  :param exp: MSExperiment for counting number of MS2 spectra
            #  :param assume_all_target: Count all PepIDs towards number of identified MS2 spectra

        void compute(PeptideIdentificationList & pep_ids, MSExperiment & exp, bool assume_all_target) except + nogil 
            # wrap-doc:
            #  Computes MS2 Identification Rate with PeptideIdentifications
            #  
            #  Stores results as a struct in a vector. Only pep-ids with target/decoy 
            #  annotation as 'target' are counted, unless assume_all_target flag is set.
            #  
            #  :param pep_ids: Input PeptideIdentifications with target/decoy annotation
            #  :param exp: MSExperiment for counting number of MS2 spectra
            #  :param assume_all_target: Count all PepIDs towards number of identified MS2 spectra

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

        libcpp_vector[IdentificationRateData] getResults() except + nogil  # wrap-doc:Returns results

cdef extern from "<OpenMS/QC/Ms2IdentificationRate.h>" namespace "OpenMS::Ms2IdentificationRate":
    
    cdef cppclass IdentificationRateData "OpenMS::Ms2IdentificationRate::IdentificationRateData":
        # wrap-doc:
        #  Structure for storing identification rate results
        
        IdentificationRateData() except + nogil 
        IdentificationRateData(IdentificationRateData &) except + nogil 
        
        Size num_peptide_identification  # wrap-doc:Number of peptide identifications
        Size num_ms2_spectra  # wrap-doc:Number of MS2 spectra
        double identification_rate  # wrap-doc:Identification rate
