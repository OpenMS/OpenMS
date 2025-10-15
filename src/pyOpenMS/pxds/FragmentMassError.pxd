from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from QCBase cimport *
from MSExperiment cimport *
from ProteinIdentification cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/QC/FragmentMassError.h>" namespace "OpenMS":
    
    cdef cppclass FragmentMassError(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  QC metric for fragment mass error (FME)
        #  
        #  Calculates the fragment mass error in ppm or Da for each match of a 
        #  theoretical spectrum to an experimental spectrum.

        FragmentMassError() except + nogil 
        FragmentMassError(FragmentMassError &) except + nogil 

        void compute(PeptideIdentificationList & pep_ids, 
                     SearchParameters & search_params, 
                     MSExperiment & exp, 
                     SpectraMap & map_to_spectrum,
                     ToleranceUnit tolerance_unit,
                     double tolerance) except + nogil 
            # wrap-doc:
            #  Calculates fragment mass error and annotates PeptideIdentifications
            #  
            #  For each match of a theoretical spectrum to an experimental spectrum, 
            #  the fragment mass error (FME) in ppm is calculated and stored as 
            #  metavalue in the PeptideIdentification.

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

        libcpp_vector[Statistics] getResults() except + nogil  # wrap-doc:Returns results

cdef extern from "<OpenMS/QC/FragmentMassError.h>" namespace "OpenMS::FragmentMassError":
    
    cdef cppclass Statistics "OpenMS::FragmentMassError::Statistics":
        # wrap-doc:
        #  Structure for storing fragment mass error statistics
        
        Statistics() except + nogil 
        Statistics(Statistics &) except + nogil 
        
        double average_ppm  # wrap-doc:Average fragment mass error in ppm
        double variance_ppm  # wrap-doc:Variance of fragment mass error in ppm
