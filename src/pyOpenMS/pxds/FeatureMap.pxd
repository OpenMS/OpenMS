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

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/FeatureMap.h>" namespace "OpenMS":

    cdef cppclass FeatureMap(UniqueIdInterface, DocumentIdentifier, RangeManager2, MetaInfoInterface):

        # wrap-inherits:
        #   UniqueIdInterface
        #   DocumentIdentifier
        #   RangeManager2
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
        Feature operator[](int)      nogil except + #wrap-upper-limit:size()
        void push_back(Feature spec) nogil except +
        void push_back(MRMFeature spec) nogil except +

        void sortByIntensity() nogil except +
        void sortByIntensity(bool reverse) nogil except +
        void sortByPosition() nogil except +
        void sortByRT() nogil except +
        void sortByMZ() nogil except +
        void sortByOverallQuality() nogil except +

        void swap(FeatureMap &) nogil except +
        void swapFeaturesOnly(FeatureMap swapfrom) nogil except +
        void clear() nogil except +
        void clear(bool clear_meta_data) nogil except +

        FeatureMap operator+(FeatureMap) nogil except +
        FeatureMap iadd(FeatureMap) nogil except + # wrap-as:operator+=

        void updateRanges() nogil except +

        libcpp_vector[ProteinIdentification] getProteinIdentifications() nogil except+
        void setProteinIdentifications(libcpp_vector[ProteinIdentification]) nogil except+

        libcpp_vector[PeptideIdentification] getUnassignedPeptideIdentifications() nogil except+
        void setUnassignedPeptideIdentifications(libcpp_vector[PeptideIdentification]) nogil except+

        Size applyMemberFunction(Size(* fun)()) nogil except +# wrap-ignore

        libcpp_vector[DataProcessing] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[DataProcessing])   nogil except +

        void setPrimaryMSRunPath(StringList& s) nogil except +
        void getPrimaryMSRunPath(StringList& toFill) nogil except +

        libcpp_vector[Feature].iterator begin() nogil except +    # wrap-iter-begin:__iter__(Feature)
        libcpp_vector[Feature].iterator end()   nogil except +    # wrap-iter-end:__iter__(Feature)

