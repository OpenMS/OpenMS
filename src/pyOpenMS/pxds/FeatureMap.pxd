from libcpp cimport bool
from Feature cimport *
from MRMFeature cimport *
from UniqueIdInterface cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from DataProcessing cimport *
from MetaInfoInterface cimport *
from DocumentIdentifier cimport *
from RangeManager cimport *
from MSExperiment cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/FeatureMap.h>" namespace "OpenMS":

    cdef cppclass FeatureMap(UniqueIdInterface, DocumentIdentifier, RangeManagerRtMzInt, MetaInfoInterface):

        # wrap-inherits:
        #   UniqueIdInterface
        #   DocumentIdentifier
        #   RangeManagerRtMzInt
        #   MetaInfoInterface
        #
        # wrap-instances:
        #   FeatureMap := FeatureMap
        #
        # wrap-doc:
        #   A container for features.
        #   -----
        #   A feature map is a container holding features, which represent
        #   chemical entities (peptides, proteins, small molecules etc.) found
        #   in an LC-MS/MS experiment.
        #   -----
        #   This class supports direct iteration in Python.

        FeatureMap() nogil except +
        FeatureMap(FeatureMap &) nogil except +

        bool operator==(FeatureMap) nogil except +
        bool operator!=(FeatureMap) nogil except +

        int size()  nogil except +
        Feature & operator[](int)      nogil except + #wrap-upper-limit:size()
        void push_back(Feature spec) nogil except +
        void push_back(MRMFeature spec) nogil except +

        void sortByIntensity() nogil except + # wrap-doc:Sorts the peaks according to ascending intensity
        void sortByIntensity(bool reverse) nogil except + # wrap-doc:Sorts the peaks according to ascending intensity. Order is reversed if argument is `true` ( reverse = true )
        void sortByPosition() nogil except + # wrap-doc:Sorts features by position. Lexicographical comparison (first RT then m/z) is done
        void sortByRT() nogil except + # wrap-doc:Sorts features by RT position
        void sortByMZ() nogil except + # wrap-doc:Sorts features by m/z position
        void sortByOverallQuality() nogil except + # wrap-doc:Sorts features by ascending overall quality. Order is reversed if argument is `true` ( reverse = true )

        void swap(FeatureMap &) nogil except + 
        void swapFeaturesOnly(FeatureMap swapfrom) nogil except + # wrap-doc:Swaps the feature content (plus its range information) of this map 
        void clear() nogil except + # wrap-doc:Clears all data and meta data
        void clear(bool clear_meta_data) nogil except + # wrap-doc:Clears all data and meta data. If 'true' is passed as an argument, all meta data is cleared in addition to the data

        FeatureMap operator+(FeatureMap) nogil except +
        FeatureMap iadd(FeatureMap) nogil except + # wrap-as:operator+=

        void updateRanges() nogil except + # TODO

        libcpp_vector[ProteinIdentification] getProteinIdentifications() nogil except+
        void setProteinIdentifications(libcpp_vector[ProteinIdentification]) nogil except+ # wrap-doc:Sets the protein identifications

        libcpp_vector[PeptideIdentification] getUnassignedPeptideIdentifications() nogil except+
        void setUnassignedPeptideIdentifications(libcpp_vector[PeptideIdentification]) nogil except+ # wrap-doc:Sets the unassigned peptide identifications

        Size applyMemberFunction(Size(* fun)()) nogil except +# wrap-ignore

        libcpp_vector[DataProcessing] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[DataProcessing])   nogil except + # wrap-doc:Sets the description of the applied data processing

        void setPrimaryMSRunPath(StringList& s) nogil except + # wrap-doc:Sets the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
        void setPrimaryMSRunPath(StringList& s, MSExperiment& e) nogil except + # wrap-doc:Sets the file path to the primary MS run using the mzML annotated in the MSExperiment argument `e`
        void getPrimaryMSRunPath(StringList& toFill) nogil except + # wrap-doc:Returns the file path to the first MS run

        libcpp_vector[Feature].iterator begin() nogil except +    # wrap-iter-begin:__iter__(Feature)
        libcpp_vector[Feature].iterator end()   nogil except +    # wrap-iter-end:__iter__(Feature)

