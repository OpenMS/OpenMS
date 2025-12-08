from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *
from ProteinIdentification cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/QC/MissedCleavages.h>" namespace "OpenMS":

    cdef cppclass MissedCleavages(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  This class is a metric for the QualityControl TOPP Tool.
        #
        #  This class counts the number of missed cleavages per PeptideIdentification
        #  given a FeatureMap and returns an agglomeration statistic (observed counts).
        #  Additionally the PeptideHits in the FeatureMap are augmented with meta information.

        MissedCleavages() except + nogil
        MissedCleavages(MissedCleavages &) except + nogil

        void compute(FeatureMap & fmap) except + nogil
            # wrap-doc:
            #  Counts the number of missed cleavages per PeptideIdentification.
            #
            #  The result is a key/value map: #missed_cleavages --> counts.
            #  Additionally the first PeptideHit in each PeptideIdentification of the
            #  FeatureMap is annotated with metavalue 'missed_cleavages'.
            #  The protease and digestion parameters are taken from the first
            #  ProteinIdentication within the FeatureMap itself.
            #
            #  :param fmap: FeatureMap with Peptide and ProteinIdentifications

        void compute(libcpp_vector[ProteinIdentification] & prot_ids, PeptideIdentificationList & pep_ids) except + nogil
            # wrap-doc:
            #  Counts the number of missed cleavages per PeptideIdentification.
            #
            #  :param prot_ids: ProteinIdentifications containing search parameters
            #  :param pep_ids: PeptideIdentifications to analyze

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

        libcpp_vector[libcpp_map[UInt32, UInt32]] getResults() except + nogil  # wrap-doc:Returns results as maps of number of missed_cleavages to counts
