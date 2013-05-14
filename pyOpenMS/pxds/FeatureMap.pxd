from libcpp cimport bool
from Feature cimport *
from UniqueIdInterface cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from DataProcessing cimport *

cdef extern from "<OpenMS/KERNEL/FeatureMap.h>" namespace "OpenMS":

    cdef cppclass FeatureMap[FeatureT](UniqueIdInterface):

        # wrap-inherits:
        #   UniqueIdInterface
        #
        # wrap-instances:
        #   FeatureMap := FeatureMap[Feature]

        FeatureMap() nogil except +
        FeatureMap(FeatureMap[FeatureT]) nogil except + # wrap-ignore
        int   size()  nogil except +
        Feature operator[](int)      nogil except +
        void push_back(FeatureT spec) nogil except +

        void sortByIntensity() nogil except +
        void sortByIntensity(bool reverse) nogil except +
        void sortByPosition() nogil except +
        void sortByRT() nogil except +
        void sortByMZ() nogil except +
        void sortByOverallQuality() nogil except +

        void swap(FeatureMap[FeatureT] &) nogil except +
        void clear() nogil except +
        void clear(bool clear_meta_data) nogil except +

        FeatureMap[FeatureT] operator+(FeatureMap[FeatureT]) nogil except +
        FeatureMap[FeatureT] iadd(FeatureMap[FeatureT]) nogil except + # wrap-as:operator+=

        void   updateRanges() nogil except +

        libcpp_vector[ProteinIdentification] getProteinIdentifications() nogil except+
        void setProteinIdentifications(libcpp_vector[ProteinIdentification]) nogil except+

        libcpp_vector[PeptideIdentification] getUnassignedPeptideIdentifications() nogil except+
        void setUnassignedPeptideIdentifications(libcpp_vector[PeptideIdentification]) nogil except+

        void applyMemberFunction(Size(* fun)())  nogil except + # wrap-ignore
 
        libcpp_vector[DataProcessing] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[DataProcessing])   nogil except +







