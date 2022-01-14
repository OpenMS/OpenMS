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
        #   Perform targeted feature extraction of compounds provided as table and stores them in features
        #   -----
        #   The algorithms detects quantitative features in MS1 data for a list of targets, typically small molecule/metabolite identifications
        #   Internally, it uses algorithms for targeted data analysis from the OpenSWATH pipeline
        #   In the simplest case, only CompoundName, SumFormula, Charge and RetentionTime need to be given, all other values may be zero
        #   Every combination of compound (mass), RT and charge defines one target for feature detection
        #   Output:
        #   The main output is a feature map of detected features, with annotations in meta data entries
        #   Additional outputs are the extracted chromatograms/peak groups, the assay in TraML compatible format, and transformations 
        #   that contain the error between provided and observed peaks
        #   -----
        #   Usage:
        #       exp = MSExperiment()
        #       MzMLFile().load(path_to_file, exp)
        #       ff = FeatureFinderAlgorithmMetaboIdent()
        #       ff.setMSData(exp)
        #
        #       fm = FeatureMap() # detected features will be stored here
        #
        #       library = []
        #       # fill library with compounds: FeatureFinderMetaboIdentCompound(name, formula, mass, [charges] [RTs_in_sec], [RT_ranges], [isotope distributions])
        #       # e.g. FeatureFinderMetaboIdentCompound('glucose','C6H12O6', 0.0, [-1], [123.4], [0.0], [0.0])
        #       
        #       params = ff.getParameters() # optional!
        #       params[param_name] = new_value # e.g. params[b'extract:n_isotopes'] = 3
        #       ff.setParameters(params)
        #
        #       ff.run(library, fm, path_to_file)

        FeatureFinderAlgorithmMetaboIdent() nogil except +

        void setMSData(MSExperiment & input) nogil except + #wrap-doc:Sets spectra
        const MSExperiment& getMSData() nogil except + #wrap-doc:Returns spectra

        void run(const libcpp_vector[ FeatureFinderMetaboIdentCompound ] metaboIdentTable, FeatureMap& features, String spectra_path) nogil except +
        # wrap-doc:
        #    Run feature extraction. spectra_path get's annotated as primaryMSRunPath in the resulting feature map.

        MSExperiment& getChromatograms() nogil except + #wrap-doc:Retrieves chromatograms (empty if run was not executed)

        const TargetedExperiment& getLibrary () nogil except + #wrap-doc:Retrieves the assay library (e.g., to store as TraML, empty if run was not executed)
        
        const TransformationDescription& getTransformations() nogil except + #wrap-doc:Retrieves deviations between provided coordinates and extacted ones (e.g., to store as TrafoXML or for plotting)

        size_t getNShared() nogil except + #wrap-doc:Retrieves number of features with shared identifications


cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h>" namespace "OpenMS::FeatureFinderAlgorithmMetaboIdent":

    cdef cppclass FeatureFinderMetaboIdentCompound "OpenMS::FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound":
                
        FeatureFinderMetaboIdentCompound(String name,
            String formula,
            double mass,
            libcpp_vector[ int ] charges,
            libcpp_vector[ double ] rts,
            libcpp_vector[ double ] rt_ranges,
            libcpp_vector[ double ] iso_distrib) nogil except +
        
        # wrap-doc:
        #    Represents a compound in the ID table
