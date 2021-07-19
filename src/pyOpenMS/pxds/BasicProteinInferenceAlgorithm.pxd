from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from ProteinIdentification cimport *
from ProteinIdentification cimport *
from PeptideHit cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>" namespace "OpenMS":

    cdef cppclass BasicProteinInferenceAlgorithm(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        #
        # wrap-doc:
        #     Algorithm class that implements simple protein inference by aggregation of peptide scores.
        #     It has multiple parameter options like the aggregation method, when to distinguish peptidoforms,
        #     and if you want to use shared peptides ("use_shared_peptides").
        #     First, the best PSM per spectrum is used, then only the best PSM per peptidoform is aggregated.
        #     Peptidoforms can optionally be distinguished via the treat_X_separate parameters:
        #     - Modifications (modified sequence string)
        #     - Charge states
        #     The algorithm assumes posteriors or posterior error probabilities and converts to posteriors initially.
        #     Possible aggregation methods that can be set via the parameter "aggregation_method" are:
        #     - "maximum" (default)
        #     - "sum"
        #     - "product" (ignoring zeroes)
        #     Annotation of the number of peptides used for aggregation can be disabled (see parameters).
        #     Supports multiple runs but goes through them one by one iterating over the full PeptideIdentification vector.

        BasicProteinInferenceAlgorithm() nogil except +

        BasicProteinInferenceAlgorithm(BasicProteinInferenceAlgorithm) nogil except + #wrap-ignore

        void run(libcpp_vector[ PeptideIdentification ] & pep_ids,
                 libcpp_vector[ ProteinIdentification ] & prot_ids) nogil except +

        void run(libcpp_vector[ PeptideIdentification ] & pep_ids,
                                ProteinIdentification & prot_id) nogil except +

cdef extern from "<OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>" namespace "OpenMS::BasicProteinInferenceAlgorithm":
    cdef enum AggregationMethod "OpenMS::BasicProteinInferenceAlgorithm::AggregationMethod":
        # wrap-doc:
        #     Aggregation method
        # wrap-attach:
        #     BasicProteinInferenceAlgorithm
        PROD # wrap-doc:Aggregate by product (ignore zeroes)
        SUM # wrap-doc:Aggregate by summing
        MAXIMUM # wrap-doc:Aggregate by maximum

