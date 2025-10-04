from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from Feature cimport *
from OpenSwathScoring cimport *

cdef extern from "<OpenMS/KERNEL/MRMFeature.h>" namespace "OpenMS":

    cdef cppclass MRMFeature(Feature):
        #
        # wrap-inherits:
        #   Feature
        #
        # wrap-doc:
        #  A multi-chromatogram MRM (Multiple Reaction Monitoring) feature representing a peak group
        #  
        #  An MRMFeature represents a detected signal across multiple chromatograms in targeted proteomics
        #  experiments (MRM/SRM). It contains corresponding features from individual transitions, where each
        #  transition is represented as a Feature object. This class is essential for analyzing peak groups
        #  in targeted MS experiments.
        #  
        #  The MRMFeature stores:
        #  - Individual transition features (via addFeature/getFeature)
        #  - Precursor features for MS1 data (via addPrecursorFeature/getPrecursorFeature)
        #  - Quality scores for the peak group (via getScores/setScores)
        #  
        #  Example usage:
        #    >>> mrm_feature = oms.MRMFeature()
        #    >>> # Add a transition feature with its native ID
        #    >>> feature = oms.Feature()
        #    >>> feature.setRT(100.5)
        #    >>> feature.setMZ(500.25)
        #    >>> feature.setIntensity(10000.0)
        #    >>> mrm_feature.addFeature(feature, "transition_1")
        #    >>> # Retrieve the feature by its ID
        #    >>> retrieved_feature = mrm_feature.getFeature("transition_1")
        #    >>> print(retrieved_feature.getRT())  # Should print 100.5

        MRMFeature() except + nogil 
        MRMFeature(MRMFeature &) except + nogil 

        OpenSwath_Scores getScores() except + nogil 
            # wrap-doc:
            #  Returns the peak group quality scores
            #  
            #  :return: An object containing various quality metrics for the peak group, such as library correlation, signal-to-noise ratio, and other OpenSWATH scoring metrics

        void setScores(OpenSwath_Scores s) except + nogil 
            # wrap-doc:
            #  Sets the peak group quality scores
            #  
            #  :param s: An OpenSwath_Scores object containing quality metrics for this peak group

        Feature getFeature(String key) except + nogil 
            # wrap-doc:
            #  Retrieves a transition feature by its native ID
            #  
            #  :param key: The native ID of the transition (e.g., "transition_1" or a TRAML identifier)
            #  :return: The Feature object corresponding to this transition
            #  
            #  Raises an exception if the key is not found

        void addFeature(Feature & f, String key) except + nogil 
            # wrap-doc:
            #  Adds a transition feature to this MRM feature
            #  
            #  :param f: A Feature object representing the signal from a single transition chromatogram
            #  :param key: The native ID for this transition (should match the transition ID from your method)
            #  
            #  Each transition in an MRM experiment produces one chromatogram, which is represented as a Feature

        libcpp_vector[Feature] getFeatures() except + nogil 
            # wrap-doc:
            #  Returns all transition features in this MRM feature
            #  
            #  :return: A list of all transition features that have been added to this peak group

        void getFeatureIDs(libcpp_vector[String] & result) except + nogil 
            # wrap-doc:
            #  Retrieves the native IDs of all transition features
            #  
            #  :param result: Output parameter that will be filled with the native IDs of all transitions
            #  
            #  This is an output parameter. Pass an empty list and it will be populated with IDs

        Feature getPrecursorFeature(String key) except + nogil 
            # wrap-doc:
            #  Retrieves a precursor feature by its native ID
            #  
            #  :param key: The native ID of the precursor
            #  :return: The Feature object for the precursor (MS1 signal)
            #  
            #  Precursor features represent the MS1 signal for the peptide, if available

        void addPrecursorFeature(Feature & f, String key) except + nogil 
            # wrap-doc:
            #  Adds a precursor feature to this MRM feature
            #  
            #  :param f: A Feature object representing the MS1 precursor signal
            #  :param key: The native ID for this precursor
            #  
            #  Precursor features are optional and represent MS1-level information

        void getPrecursorFeatureIDs(libcpp_vector[String] & result) except + nogil 
            # wrap-doc:
            #  Retrieves the native IDs of all precursor features
            #  
            #  :param result: Output parameter that will be filled with precursor IDs
            #  
            #  This is an output parameter. Pass an empty list and it will be populated with IDs

        bool operator==(MRMFeature) except + nogil 
        bool operator!=(MRMFeature) except + nogil 
