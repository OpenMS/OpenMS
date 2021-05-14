from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmMetaboIdent(DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        # wrap-doc:
        #   Perform targeted feature extraction of compounds provided as table and stores them in features.
        #   -----
        #   The algorithms detects quantitative features in MS1 data for a list of targets, typically small molecule/metabolite identifications.
        #   Internally, it uses algorithms for targeted data analysis from the OpenSWATH pipeline.
        #   In the simplest case, only CompoundName, SumFormula, Charge and RetentionTime need to be given, all other values may be zero.
        #   Every combination of compound (mass), RT and charge defines one target for feature detection.
        #   Output:
        #   The main output is a feature map of detected features, with annotations in meta data entries.
        #   Additional outputs are the extracted chromatograms/peak groups, the assay in TraML compatible format, and transformations 
        #   that contain the error between provided and observed peaks.

        FeatureFinderAlgorithmMetaboIdent() nogil except +

        void setMSData(MSExperiment & input) nogil except + #wrap-doc:Set spectra
        const MSExperiment& getMSData() nogil except + #wrap-doc:Get spectra

        void run(const libcpp_vector[ Row ] metaboIdentTable, FeatureMap& features) nogil except + #wrap-doc:Run the experiment

        MSExperiment& getChromatograms() nogil except + #wrap-doc:Retrieve chromatograms (empty if run was not executed)

        const TargetedExperiment& getLibrary () nogil except + #wrap-doc:Retrieve the assay library (e.g., to store as TraML, empty if run was not executed)
        
        const TransformationDescription& getTransformations() nogil except + #wrap-doc:Retrieve deviations between provided coordinates and extacted ones (e.g., to store as TrafoXML or for plotting)

        size_t getNShared() nogil except + #wrap-doc:Retrieve number of features with shared identifications

    cdef cppclass Row "OpenMS::FeatureFinderAlgorithmMetaboIdent::Row":
        # wrap-attach:FeatureFinderAlgorithmMetaboIdent
        
        Row(String name, String formula, double mass, libcpp_vector[ int ] charges, libcpp_vector[ double ] rts, libcpp_vector[ double ] rt_ranges, libcpp_vector[ double ] iso_distrib) nogil except + #wrap-doc:Represents a compound in the ID table 

        String name
        String formula
        double mass
        libcpp_vector[ int ] charges
        libcpp_vector[ double ] rts
        libcpp_vector[ double ] rt_ranges
        libcpp_vector[ double ] iso_distrib