from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from UniqueIdInterface cimport *
from ConsensusFeature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Feature cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from PeptideIdentificationList cimport *
from DataProcessing cimport *
from Types cimport *
from DocumentIdentifier cimport *
from RangeManager cimport *
from MetaInfoInterface cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/ConsensusMap.h>" namespace "OpenMS::ConsensusMap":

    cdef cppclass ColumnHeader(MetaInfoInterface):
        # wrap-inherits:
        #  MetaInfoInterface
        #
        # wrap-doc:Describes a column (input file) in a ConsensusMap with filename, label, size, and unique ID

        String filename
        String label
        Size size
        UInt64 unique_id

        ColumnHeader() except + nogil  # wrap-doc:Default constructor
        ColumnHeader(ColumnHeader &) except + nogil  # wrap-doc:Copy constructor

    # for msvc++ compiler, see addons/ConsensusMap.pyx
    # ... forgot why Map[..] did not work
    ctypedef libcpp_map[UInt64, ColumnHeader] ColumnHeaders "OpenMS::ConsensusMap::ColumnHeaders"
    ctypedef libcpp_map[UInt64, ColumnHeader].iterator ColumnHeaders_iterator "OpenMS::ConsensusMap::ColumnHeaders::iterator"

cdef extern from "<OpenMS/KERNEL/ConsensusMap.h>" namespace "OpenMS":

    cdef cppclass ConsensusMap(UniqueIdInterface, DocumentIdentifier, RangeManagerRtMzInt, MetaInfoInterface):

        # wrap-inherits:
        #  UniqueIdInterface
        #  DocumentIdentifier
        #  RangeManagerRtMzInt
        #  MetaInfoInterface
        #
        # wrap-doc:
        #  A container for consensus elements.
        #
        #  A ConsensusMap is a container holding 2-dimensional consensus elements
        #  (ConsensusFeature) which in turn represent analytes that have been
        #  quantified across multiple LC-MS/MS experiments. Each analyte in a
        #  ConsensusFeature is linked to its original LC-MS/MS run, the links are
        #  maintained by the ConsensusMap class.
        #  The map is implemented as a vector of elements of type ConsensusFeature.
        #
        #  To be consistent, all maps who are referenced by ConsensusFeature objects
        #  (through a unique id) need to be registered in this class.
        #
        #  This class supports direct iteration in Python.

        ConsensusMap() except + nogil  # wrap-doc:Default constructor creating an empty consensus map
        ConsensusMap(ConsensusMap &) except + nogil  # wrap-doc:Copy constructor

        bool operator==(ConsensusMap) except + nogil
        bool operator!=(ConsensusMap) except + nogil

        int size() except + nogil  # wrap-doc:Returns the number of consensus features in the map
        bool empty() except + nogil  # wrap-doc:Returns True if the map contains no consensus features
        void reserve(Size s) except + nogil  # wrap-doc:Reserves space for the given number of consensus features to avoid repeated memory allocations
        ConsensusFeature & operator[](size_t) except + nogil  #wrap-upper-limit:size()
        void push_back(ConsensusFeature spec) except + nogil  # wrap-doc:Adds a ConsensusFeature to the map

        ConsensusMap appendRows(ConsensusMap) except + nogil  # wrap-doc:Add consensus map entries as new rows
        ConsensusMap appendColumns(ConsensusMap) except + nogil  # wrap-doc:Add consensus map entries as new columns

        void clear(bool clear_meta_data) except + nogil  # wrap-doc:Clears all data and meta data
        void clear() except + nogil  # wrap-doc:Clears all consensus features and meta data from the map

        void updateRanges() except + nogil  # wrap-doc:Updates the RT, m/z, and intensity ranges based on contained consensus features

        libcpp_vector[ProteinIdentification] getProteinIdentifications(
                ) except + nogil  # wrap-doc:Returns the protein identification runs stored in this map

        void setProteinIdentifications(
                libcpp_vector[ProteinIdentification]
                ) except + nogil  # wrap-doc:Sets the protein identifications

        # PeptideIdentificationList methods
        PeptideIdentificationList getUnassignedPeptideIdentifications() except + nogil  # wrap-doc:Returns peptide identifications that are not assigned to any consensus feature
        void setUnassignedPeptideIdentifications(PeptideIdentificationList unassigned_peptide_identifications) except + nogil # wrap-doc:Sets the unassigned PeptideIdentificationList

        libcpp_vector[DataProcessing] getDataProcessing() except + nogil  # wrap-doc:Returns a const reference to the description of the applied data processing
        void setDataProcessing(libcpp_vector[DataProcessing])   except + nogil  # wrap-doc:Sets the description of the applied data processing

        void setPrimaryMSRunPath(StringList& s) except + nogil  # wrap-doc:Sets the file paths to the primary MS run (stored in ColumnHeaders)
        void setPrimaryMSRunPath(StringList& s, MSExperiment& e) except + nogil
        void getPrimaryMSRunPath(StringList& toFill) except + nogil  # wrap-doc:Returns the MS run path (stored in ColumnHeaders)

        libcpp_vector[ConsensusFeature].iterator begin(
                ) except + nogil  # wrap-iter-begin:__iter__(ConsensusFeature)
            # wrap-doc:Returns an iterator to the beginning of the consensus feature list
        libcpp_vector[ConsensusFeature].iterator end(
                ) except + nogil  # wrap-iter-end:__iter__(ConsensusFeature)
            # wrap-doc:Returns an iterator to the end of the consensus feature list

        # wrapped in ../addons/ConsensusMap.pyx:
        void applyMemberFunction(Size(* fun)()) except + nogil  # wrap-ignore

        void sortByIntensity(bool reverse) except + nogil  # wrap-doc:Sorts the peaks according to ascending intensity.
        void sortByIntensity() except + nogil  # wrap-doc:Sorts the consensus features by ascending intensity
        void sortByRT() except + nogil  # wrap-doc:Sorts the peaks according to RT position
        void sortByMZ() except + nogil  # wrap-doc:Sorts the peaks according to m/z position
        void sortByPosition() except + nogil  # wrap-doc:Lexicographically sorts the peaks by their position (First RT then m/z)
        void sortByQuality(bool reverse) except + nogil  # wrap-doc:Sorts the peaks according to ascending quality.
        void sortByQuality() except + nogil  # wrap-doc:Sorts the consensus features by ascending quality score
        void sortBySize() except + nogil  # wrap-doc:Sorts with respect to the size (number of elements)
        void sortByMaps() except + nogil  # wrap-doc:Sorts with respect to the sets of maps covered by the consensus features (lexicographically)

        # wrapped in ../addons/ConsensusMap.pyx:
        void setColumnHeaders(ColumnHeaders &)   #wrap-ignore
        ColumnHeaders & getColumnHeaders()       #wrap-ignore

        String getExperimentType() except + nogil  # wrap-doc:Non-mutable access to the experiment type
        void setExperimentType(String experiment_type) except + nogil  # wrap-doc:Mutable access to the experiment type

        void sortPeptideIdentificationsByMapIndex() except + nogil  # wrap-doc:Sorts PeptideIdentifications of consensus features with respect to their map index.

