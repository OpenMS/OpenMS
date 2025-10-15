from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *
from FASTAFile cimport *

cdef extern from "<OpenMS/QC/Contaminants.h>" namespace "OpenMS":
    
    cdef cppclass Contaminants(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  This class is a metric for the QualityControl TOPP tool.
        #  
        #  This class checks whether a peptide is a contaminant (given a protein DB) 
        #  and adds that result as metavalue "is_contaminant" to the first hit of 
        #  each PeptideIdentification.

        Contaminants() except + nogil 
        Contaminants(Contaminants &) except + nogil 

        void compute(FeatureMap & features, libcpp_vector[FASTAEntry] & contaminants) except + nogil 
            # wrap-doc:
            #  Checks if the peptides are in the contaminant database.
            #  
            #  "is_contaminant" metavalue is added to the first hit of each 
            #  PeptideIdentification of each feature and to the first hit of all 
            #  unassigned PeptideIdentifications.

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

        libcpp_vector[ContaminantsSummary] getResults() except + nogil  # wrap-doc:Returns results

cdef extern from "<OpenMS/QC/Contaminants.h>" namespace "OpenMS::Contaminants":
    
    cdef cppclass ContaminantsSummary "OpenMS::Contaminants::ContaminantsSummary":
        # wrap-doc:
        #  Structure for storing results
        
        ContaminantsSummary() except + nogil 
        ContaminantsSummary(ContaminantsSummary &) except + nogil 
        
        double assigned_contaminants_ratio  # wrap-doc:(#contaminants in assigned / #peptides in assigned)
        double unassigned_contaminants_ratio  # wrap-doc:(#contaminants in unassigned / #peptides in unassigned)
        double all_contaminants_ratio  # wrap-doc:(#all contaminants / #peptides in all)
        double assigned_contaminants_intensity_ratio  # wrap-doc:(intensity of contaminants in assigned / intensity of peptides in assigned)
