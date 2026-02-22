from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport uint32_t, uint64_t, int64_t
from OpenMS cimport *
from QCBase cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from WindowMower cimport *

cdef extern from "<OpenMS/QC/PSMExplainedIonCurrent.h>" namespace "OpenMS":

    cdef cppclass PSMExplainedIonCurrent(QCBase):
        # wrap-inherits:
        #    QCBase

        # Nested struct for results
        cppclass Statistics:
            double average_correctness
            double variance_correctness

        PSMExplainedIonCurrent() except +
        PSMExplainedIonCurrent(PSMExplainedIonCurrent &) except +

        void compute(FeatureMap& fmap, 
                     const MSExperiment& exp, 
                     const QCBase.SpectraMap& map_to_spectrum,
                     ToleranceUnit tolerance_unit,
                     double tolerance) except +
        
        void compute(PeptideIdentificationList& pep_ids, 
                     const ProteinIdentification.SearchParameters& search_params,
                     const MSExperiment& exp, 
                     const QCBase.SpectraMap& map_to_spectrum,
                     ToleranceUnit tolerance_unit,
                     double tolerance) except +

        const vector[Statistics]& getResults() const

        String getName() const

        QCBase.Status requirements() const
EOF