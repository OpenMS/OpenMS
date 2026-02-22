cat > /workspaces/OpenMS/pyOpenMS/pxds/Ms2IdentificationRate.pxd << 'EOF'
# distutils: language = c++

from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport uint32_t, uint64_t, int64_t
from OpenMS cimport *
from QCBase cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from MzTabMetaData cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/QC/Ms2IdentificationRate.h>" namespace "OpenMS":

    cdef cppclass Ms2IdentificationRate(QCBase):
        # wrap-inherits:
        #    QCBase

        # Nested struct for results
        cppclass IdentificationRateData:
            uint64_t num_peptide_identification
            uint64_t num_ms2_spectra
            double identification_rate

        Ms2IdentificationRate() except +
        Ms2IdentificationRate(Ms2IdentificationRate &) except +

        void compute(const FeatureMap& feature_map, 
                     const MSExperiment& exp, 
                     bool assume_all_target) except +
        
        void compute(const PeptideIdentificationList& pep_ids, 
                     const MSExperiment& exp, 
                     bool assume_all_target) except +

        const vector[IdentificationRateData]& getResults() const

        String getName() const

        QCBase.Status requirements() const

        void addMetaDataMetricsToMzTab(MzTabMetaData& meta) const
EOF