from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from ProteinIdentification cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/QC/PSMExplainedIonCurrent.h>" namespace "OpenMS":
    
    cdef cppclass PSMExplainedIonCurrent(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  QC metric for PSM explained ion current
        #  
        #  Computes PSM explained ion current (only of the first PeptideHit of each PepID).
        #  To calculate PSMExplainedIonCurrent the theoretical spectrum is generated and 
        #  matched with the original one.

        PSMExplainedIonCurrent() except + nogil 
        PSMExplainedIonCurrent(PSMExplainedIonCurrent &) except + nogil 

        void compute(FeatureMap & fmap, 
                     MSExperiment & exp, 
                     SpectraMap & map_to_spectrum,
                     ToleranceUnit tolerance_unit,
                     double tolerance) except + nogil 
            # wrap-doc:
            #  Computes PSMExplainedIonCurrent (only of the first PeptideHit of each PepID)
            #  
            #  Stores average and variance of PSMExplainedIonCurrent as a struct and 
            #  stores it in the results vector. Each PSMExplainedIonCurrent is also 
            #  stored in the first PeptideHit of the corresponding PeptideIdentification 
            #  as metavalue "PSM_correctness".

        void compute(PeptideIdentificationList & pep_ids, 
                     SearchParameters & search_params, 
                     MSExperiment & exp, 
                     SpectraMap & map_to_spectrum,
                     ToleranceUnit tolerance_unit,
                     double tolerance) except + nogil 
            # wrap-doc:
            #  Computes PSMExplainedIonCurrent with PeptideIdentification input
            #  
            #  Same as above, but with PeptideIdentification + SearchParameter input 
            #  instead of FeatureMap.

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

        libcpp_vector[Statistics] getResults() except + nogil  # wrap-doc:Returns results

cdef extern from "<OpenMS/QC/PSMExplainedIonCurrent.h>" namespace "OpenMS::PSMExplainedIonCurrent":
    
    cdef cppclass Statistics "OpenMS::PSMExplainedIonCurrent::Statistics":
        # wrap-doc:
        #  Structure for storing results: average and variance over all PSMs
        
        Statistics() except + nogil 
        Statistics(Statistics &) except + nogil 
        
        double average_correctness  # wrap-doc:Average PSM explained ion current
        double variance_correctness  # wrap-doc:Variance of PSM explained ion current
