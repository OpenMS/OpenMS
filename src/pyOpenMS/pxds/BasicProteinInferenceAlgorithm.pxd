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
        #   Algorithm class that implements simple protein inference by aggregation of peptide scores.
        #   -----
        #   It has multiple parameter options like the aggregation method, when to distinguish peptidoforms,
        #   and if you want to use shared peptides ("use_shared_peptides").
        #   First, the best PSM per spectrum is used, then only the best PSM per peptidoform is aggregated.
        #   Peptidoforms can optionally be distinguished via the treat_X_separate parameters:
        #   - Modifications (modified sequence string)
        #   - Charge states
        #   The algorithm assumes posteriors or posterior error probabilities and converts to posteriors initially.
        #   Possible aggregation methods that can be set via the parameter "aggregation_method" are:
        #   - "maximum" (default)
        #   - "sum"
        #   - "product" (ignoring zeroes)
        #   Annotation of the number of peptides used for aggregation can be disabled (see parameters).
        #   Supports multiple runs but goes through them one by one iterating over the full PeptideIdentification vector.
        #   -----
        #   Usage:
        #     from pyopenms import *
        #     from urllib.request import urlretrieve
        #     urlretrieve("https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/tests/class_tests/openms/data/newMergerTest_out.idXML", "BasicProteinInference_test.idXML")
        #     proteins = []
        #     peptides = []
        #     idf = IdXMLFile()
        #     idf.load("BasicProteinInference_test.idXML", proteins, peptides);
        #     bpia = BasicProteinInferenceAlgorithm()
        #     p = bpia.getParameters();
        #     p.setValue("min_peptides_per_protein", 0);
        #     bpia.setParameters(p);
        #     bpia.run(peptides, proteins);
        #     #
        #     print(proteins[0].getHits()[0].getScore()) # 0.6
        #     print(proteins[0].getHits()[5].getScore()) # 0.9
        #     print(proteins[0].getHits()[0].getMetaValue("nr_found_peptides")) # 1
        #     print(proteins[0].getHits()[3].getMetaValue("nr_found_peptides")) # 2
        #   -----

        BasicProteinInferenceAlgorithm() nogil except +

        BasicProteinInferenceAlgorithm(BasicProteinInferenceAlgorithm) nogil except + #wrap-ignore

        void run(libcpp_vector[ PeptideIdentification ] & pep_ids,
                 libcpp_vector[ ProteinIdentification ] & prot_ids) nogil except +
          # wrap-doc:
          #   Performs inference
          #   -----
          #   Annotation of protein groups is currently only possible for a single protein ID run
          #   -----
          #   :param pep_ids: Vector of peptide identifications
          #   :param prot_ids: Vector of protein identifications
          #   :return: Writes its results into prot_ids

        void run(libcpp_vector[ PeptideIdentification ] & pep_ids,
                                ProteinIdentification & prot_id) nogil except +
          # wrap-doc:
          #   Performs inference
          #   -----
          #   Optionally adds indistinguishable protein groups with separate scores, too
          #   Currently only takes first proteinID run and all peptides
          #   -----
          #   :param pep_ids: Vector of peptide identifications
          #   :param prot_id: Peptide identification
          #   :return: Writes its results into prot_ids

cdef extern from "<OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>" namespace "OpenMS::BasicProteinInferenceAlgorithm":
    cdef enum AggregationMethod "OpenMS::BasicProteinInferenceAlgorithm::AggregationMethod":
        # wrap-doc:
        #     Aggregation method
        # wrap-attach:
        #     BasicProteinInferenceAlgorithm
        PROD # wrap-doc:Aggregate by product (ignore zeroes)
        SUM # wrap-doc:Aggregate by summing
        MAXIMUM # wrap-doc:Aggregate by maximum

