from libcpp cimport bool
from Feature cimport *
from MRMFeature cimport *
from UniqueIdInterface cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from PeptideIdentificationList cimport *
from DataProcessing cimport *
from MetaInfoInterface cimport *
from DocumentIdentifier cimport *
from RangeManager cimport *
from MSExperiment cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/FeatureMap.h>" namespace "OpenMS":

    cdef cppclass FeatureMap(UniqueIdInterface, DocumentIdentifier, RangeManagerRtMzInt, MetaInfoInterface):

        # wrap-inherits:
        #  UniqueIdInterface
        #  DocumentIdentifier
        #  RangeManagerRtMzInt
        #  MetaInfoInterface
        #
        # wrap-instances:
        #  FeatureMap := FeatureMap
        #
        # wrap-doc:
        #  A container for LC-MS features with metadata and identification information
        #  
        #  FeatureMap is one of the core data structures in OpenMS for storing detected features
        #  from LC-MS experiments. A feature represents a detected chemical entity (peptide, protein,
        #  metabolite, etc.) with its elution profile and mass information.
        #  
        #  Key capabilities:
        #  - Store and manage Feature objects (detected analytes)
        #  - Associate protein and peptide identifications with features
        #  - Sort features by various criteria (RT, m/z, intensity, quality)
        #  - Store experimental metadata and data processing information
        #  - Support direct iteration and indexing in Python
        #  
        #  Example usage:
        #    >>> feature_map = oms.FeatureMap()
        #    >>> # Add a feature
        #    >>> feature = oms.Feature()
        #    >>> feature.setRT(1234.5)
        #    >>> feature.setMZ(445.678)
        #    >>> feature.setIntensity(100000.0)
        #    >>> feature_map.push_back(feature)
        #    >>> # Access features
        #    >>> print(f"Number of features: {feature_map.size()}")
        #    >>> first_feature = feature_map[0]
        #    >>> # Sort by RT
        #    >>> feature_map.sortByRT()
        #    >>> # Iterate over features
        #    >>> for feat in feature_map:
        #    ...     print(f"RT: {feat.getRT()}, m/z: {feat.getMZ()}")

        FeatureMap() except + nogil 
        FeatureMap(FeatureMap &) except + nogil 

        bool operator==(FeatureMap) except + nogil 
        bool operator!=(FeatureMap) except + nogil 

        int size()  except + nogil 
            # wrap-doc:
            #  Returns the number of features in the map
            #  
            #  Returns:
            #    int: Number of features stored in this container

        Feature & operator[](size_t)      except + nogil  #wrap-upper-limit:size()
        
        void push_back(Feature spec) except + nogil 
            # wrap-doc:
            #  Adds a Feature to the map
            #  
            #  Args:
            #    spec (Feature): The feature to add to the map

        void push_back(MRMFeature spec) except + nogil 
            # wrap-doc:
            #  Adds an MRMFeature to the map
            #  
            #  Args:
            #    spec (MRMFeature): The MRM feature to add to the map 

        void sortByIntensity() except + nogil 
            # wrap-doc:
            #  Sorts features by ascending intensity
            #  
            #  Note:
            #    After sorting, features can be accessed in order from lowest to highest intensity

        void sortByIntensity(bool reverse) except + nogil 
            # wrap-doc:
            #  Sorts features by intensity with optional reverse order
            #  
            #  Args:
            #    reverse (bool): If True, sorts in descending order (highest to lowest intensity)

        void sortByPosition() except + nogil 
            # wrap-doc:
            #  Sorts features by position using lexicographical comparison
            #  
            #  Note:
            #    Compares RT first, then m/z for features with the same RT

        void sortByRT() except + nogil 
            # wrap-doc:
            #  Sorts features by retention time (RT) in ascending order
            #  
            #  Note:
            #    This is useful for time-based analysis or visualization

        void sortByMZ() except + nogil 
            # wrap-doc:
            #  Sorts features by mass-to-charge ratio (m/z) in ascending order
            #  
            #  Note:
            #    Useful for mass-based grouping or analysis

        void sortByOverallQuality() except + nogil 
            # wrap-doc:
            #  Sorts features by overall quality score in ascending order
            #  
            #  Note:
            #    Higher quality scores indicate better feature detection confidence

        void swap(FeatureMap &) except + nogil  
        void swapFeaturesOnly(FeatureMap swapfrom) except + nogil  # wrap-doc:Swaps the feature content (plus its range information) of this map 
        void clear() except + nogil 
            # wrap-doc:
            #  Clears all feature data and metadata
            #  
            #  Note:
            #    After calling this, the map will be empty (size() returns 0)

        void clear(bool clear_meta_data) except + nogil 
            # wrap-doc:
            #  Clears feature data and optionally metadata
            #  
            #  Args:
            #    clear_meta_data (bool): If True, also clears all metadata; if False, keeps metadata

        FeatureMap operator+(FeatureMap) except + nogil 
        FeatureMap iadd(FeatureMap) except + nogil  # wrap-as:operator+=

        void updateRanges() except + nogil  # TODO

        libcpp_vector[ProteinIdentification] getProteinIdentifications() except + nogil
            # wrap-doc:
            #  Returns the protein identification runs stored in this map
            #  
            #  Returns:
            #    list of ProteinIdentification: Protein identification data from database searches
            #  
            #  Note:
            #    Protein identifications contain metadata about search parameters and protein hits

        void setProteinIdentifications(libcpp_vector[ProteinIdentification]) except + nogil 
            # wrap-doc:
            #  Sets the protein identifications for this map
            #  
            #  Args:
            #    protein_ids (list of ProteinIdentification): Protein identification results to associate with this map

        PeptideIdentificationList getUnassignedPeptideIdentifications() except + nogil
            # wrap-doc:
            #  Returns peptide identifications that are not assigned to any feature
            #  
            #  Returns:
            #    list of PeptideIdentification: Unassigned peptide identification results
            #  
            #  Note:
            #    These are peptide IDs that could not be matched to features, possibly due to
            #    feature detection issues or filtering

        void setUnassignedPeptideIdentifications(PeptideIdentificationList) except + nogil 
            # wrap-doc:
            #  Sets the unassigned peptide identifications
            #  
            #  Args:
            #    peptide_ids (list of PeptideIdentification): Peptide IDs not assigned to features

        Size applyMemberFunction(Size(* fun)()) except + nogil # wrap-ignore

        libcpp_vector[DataProcessing] getDataProcessing() except + nogil 
        void setDataProcessing(libcpp_vector[DataProcessing])   except + nogil  # wrap-doc:Sets the description of the applied data processing

        void setPrimaryMSRunPath(StringList& s) except + nogil  # wrap-doc:Sets the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
        void setPrimaryMSRunPath(StringList& s, MSExperiment& e) except + nogil  # wrap-doc:Sets the file path to the primary MS run using the mzML annotated in the MSExperiment argument `e`
        void getPrimaryMSRunPath(StringList& toFill) except + nogil  # wrap-doc:Returns the file path to the first MS run

        libcpp_vector[Feature].iterator begin() except + nogil     # wrap-iter-begin:__iter__(Feature)
        libcpp_vector[Feature].iterator end()   except + nogil     # wrap-iter-end:__iter__(Feature)

