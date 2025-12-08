from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
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
            #  The enzyme and number of missed cleavages used to digest the given
            #  protein DB is taken from ProteinIdentification[0].getSearchParameters().
            #
            #  :param features: Input FeatureMap with peptideidentifications of features
            #  :param contaminants: Vector of FASTAEntries to digest and check against

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
        libcpp_pair[Int64, Int64] empty_features  # wrap-doc:(features without peptideidentification or with peptideidentifications but without hits; all features)
