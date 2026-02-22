from libcpp.vector cimport vector
from OpenMS cimport *
from QCBase cimport *

cdef extern from "<OpenMS/QC/IdentificationSummary.h>" namespace "OpenMS":

    cdef cppclass IdentificationSummary(QCBase):

        IdentificationSummary() except +

        # Nested classes
        cppclass UniqueID:
            unsigned int count
            float fdr_threshold

        cppclass Result:
            unsigned int peptide_spectrum_matches
            UniqueID unique_peptides
            UniqueID unique_proteins
            float missed_cleavages_mean
            double protein_hit_scores_mean
            double peptide_length_mean

        # Methods
        Result compute(vector[ProteinIdentification]& prot_ids,
                       PeptideIdentificationList& pep_ids) except +

        String getName() const

        QCBase.Status requirements() const