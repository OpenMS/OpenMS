from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from TargetedExperiment cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from ProgressLogger cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from FeatureMap cimport * 

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderIdentificationAlgorithm(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        #
        # wrap-doc:
        #     Algorithm class for FeatureFinderIdentification
        #     External IDs (peptides_ext, proteins_ext) may be empty, 
        #     in which case no machine learning or FDR estimation will be performed.
        #     Optional seeds from e.g. untargeted FeatureFinders can be added with
        #     seeds.
        #     Results will be written to features .
        #     Caution: peptide IDs will be shrunk to best hit, FFid metavalues added
        #     and potential seed IDs added.
        #

        FeatureFinderIdentificationAlgorithm() nogil except +

        FeatureFinderIdentificationAlgorithm(FeatureFinderIdentificationAlgorithm) nogil except + #wrap-ignore

        void run(libcpp_vector[ PeptideIdentification ] peptides,
                 libcpp_vector[ ProteinIdentification ] & proteins,
                 libcpp_vector[ PeptideIdentification ] peptides_ext,
                 libcpp_vector[ ProteinIdentification ] proteins_ext,
                 FeatureMap & features, 
                 FeatureMap & seeds) nogil except + 
                 # wrap-doc:
                 #     Run feature detection
                 #     :param peptides: Vector of identified peptides
                 #     :param proteins: Vector of identified proteins
                 #     :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
                 #     :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
                 #     :param features: Feature detection results will be added here
                 #     :param seeds: Optional seeds for feature detection from e.g. untargeted FeatureFinders

        void runOnCandidates(FeatureMap & features) nogil except + # wrap-doc:Run feature detection on identified features (e.g. loaded from an IdXML file)

        MSExperiment getMSData() nogil except + # wrap-doc:Return ms data as MSExperiment

        MSExperiment getChromatograms() nogil except + # wrap-doc:Return chromatogram data as MSExperiment 

        ProgressLogger getProgressLogger() nogil except + # wrap-doc:Return progress logger

        TargetedExperiment getLibrary() nogil except + # wrap-doc:Return constructed assay library

