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
from DataProcessing cimport *
from Types cimport *
from DocumentIdentifier cimport *
from RangeManager cimport *
from MetaInfoInterface cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/ConsensusMap.h>" namespace "OpenMS::ConsensusMap":

    cdef cppclass ColumnHeader(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        String filename
        String label
        Size size
        UInt64 unique_id

        ColumnHeader() nogil except +
        ColumnHeader(ColumnHeader &) nogil except +

    # for msvc++ compiler, see addons/ConsensusMap.pyx
    # ... forgot why Map[..] did not work
    ctypedef libcpp_map[UInt64, ColumnHeader] ColumnHeaders "OpenMS::ConsensusMap::ColumnHeaders"
    ctypedef libcpp_map[UInt64, ColumnHeader].iterator ColumnHeaders_iterator "OpenMS::ConsensusMap::ColumnHeaders::iterator"

cdef extern from "<OpenMS/KERNEL/ConsensusMap.h>" namespace "OpenMS":

    cdef cppclass ConsensusMap(UniqueIdInterface, DocumentIdentifier, RangeManagerRtMzInt, MetaInfoInterface):

        # wrap-inherits:
        #   UniqueIdInterface
        #   DocumentIdentifier
        #   RangeManagerRtMzInt
        #   MetaInfoInterface
        #
        # wrap-doc:
        #   A container for consensus elements.
        #   -----
        #   A ConsensusMap is a container holding 2-dimensional consensus elements
        #   (ConsensusFeature) which in turn represent analytes that have been
        #   quantified across multiple LC-MS/MS experiments. Each analyte in a
        #   ConsensusFeature is linked to its original LC-MS/MS run, the links are
        #   maintained by the ConsensusMap class.
        #   The map is implemented as a vector of elements of type ConsensusFeature.
        #   -----
        #   To be consistent, all maps who are referenced by ConsensusFeature objects
        #   (through a unique id) need to be registered in this class. 
        #   -----
        #   This class supports direct iteration in Python.

        ConsensusMap() nogil except +
        ConsensusMap(ConsensusMap &) nogil except +

        bool operator==(ConsensusMap) nogil except +
        bool operator!=(ConsensusMap) nogil except +

        int size() nogil except +
        bool empty() nogil except +
        void reserve(Size s) nogil except +
        ConsensusFeature & operator[](int) nogil except + #wrap-upper-limit:size()
        void push_back(ConsensusFeature spec) nogil except +

        ConsensusMap appendRows(ConsensusMap) nogil except + # wrap-doc:Add consensus map entries as new rows
        ConsensusMap appendColumns(ConsensusMap) nogil except + # wrap-doc:Add consensus map entries as new columns

        void clear(bool clear_meta_data) nogil except + # wrap-doc:Clears all data and meta data
        void clear() nogil except +

        void updateRanges() nogil except +

        libcpp_vector[ProteinIdentification] getProteinIdentifications(
                ) nogil except + # TODO

        void setProteinIdentifications(
                libcpp_vector[ProteinIdentification]
                ) nogil except + # wrap-doc:Sets the protein identifications

        libcpp_vector[PeptideIdentification]\
                getUnassignedPeptideIdentifications() nogil except + # TODO

        void setUnassignedPeptideIdentifications(
                libcpp_vector[PeptideIdentification]
                ) nogil except + # wrap-doc:Sets the unassigned peptide identifications

        libcpp_vector[DataProcessing] getDataProcessing() nogil except + # wrap-doc:Returns a const reference to the description of the applied data processing
        void setDataProcessing(libcpp_vector[DataProcessing])   nogil except + # wrap-doc:Sets the description of the applied data processing

        void setPrimaryMSRunPath(StringList& s) nogil except + # wrap-doc:Sets the file paths to the primary MS run (stored in ColumnHeaders)
        void setPrimaryMSRunPath(StringList& s, MSExperiment& e) nogil except +
        void getPrimaryMSRunPath(StringList& toFill) nogil except + # wrap-doc:Returns the MS run path (stored in ColumnHeaders)

        libcpp_vector[ConsensusFeature].iterator begin(
                ) nogil except + # wrap-iter-begin:__iter__(ConsensusFeature)
        libcpp_vector[ConsensusFeature].iterator end(
                ) nogil except + # wrap-iter-end:__iter__(ConsensusFeature)

        # wrapped in ../addons/ConsensusMap.pyx:
        void applyMemberFunction(Size(* fun)()) nogil except + # wrap-ignore

        void sortByIntensity(bool reverse) nogil except + # wrap-doc:Sorts the peaks according to ascending intensity.
        void sortByIntensity() nogil except +
        void sortByRT() nogil except + # wrap-doc:Sorts the peaks according to RT position
        void sortByMZ() nogil except + # wrap-doc:Sorts the peaks according to m/z position
        void sortByPosition() nogil except + # wrap-doc:Lexicographically sorts the peaks by their position (First RT then m/z)
        void sortByQuality(bool reverse) nogil except + # wrap-doc:Sorts the peaks according to ascending quality.
        void sortByQuality() nogil except +
        void sortBySize() nogil except + # wrap-doc:Sorts with respect to the size (number of elements)
        void sortByMaps() nogil except + # wrap-doc:Sorts with respect to the sets of maps covered by the consensus features (lexicographically)

        # wrapped in ../addons/ConsensusMap.pyx:
        void setColumnHeaders(ColumnHeaders &)   #wrap-ignore
        ColumnHeaders & getColumnHeaders()       #wrap-ignore

        String getExperimentType() nogil except + # wrap-doc:Non-mutable access to the experiment type
        void setExperimentType(String experiment_type) nogil except + # wrap-doc:Mutable access to the experiment type

        void sortPeptideIdentificationsByMapIndex() nogil except + # wrap-doc:Sorts PeptideIdentifications of consensus features with respect to their map index.

