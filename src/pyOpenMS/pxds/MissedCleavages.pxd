from libcpp.vector cimport vector
from libcpp.map cimport map as cppmap
from libc.stdint cimport uint32_t
from OpenMS cimport *
from QCBase cimport *
from ProteaseDigestion cimport *
from FeatureMap cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/QC/MissedCleavages.h>" namespace "OpenMS":

    cdef cppclass MissedCleavages(QCBase):
        # wrap-inherits:
        #    QCBase

        MissedCleavages() except +
        MissedCleavages(MissedCleavages &) except +

        void compute(FeatureMap& fmap) except +
        void compute(vector[ProteinIdentification]& prot_ids, 
                     PeptideIdentificationList& pep_ids) except +

        const vector[cppmap[uint32_t, uint32_t]]& getResults() const

        String getName() const

        QCBase.Status requirements() const
EOF