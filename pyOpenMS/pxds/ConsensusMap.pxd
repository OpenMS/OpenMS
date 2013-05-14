from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from UniqueIdInterface cimport *
from ConsensusFeature cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from DataProcessing cimport *
from Types cimport *
from Map cimport *

cdef extern from "<OpenMS/KERNEL/ConsensusMap.h>" namespace "OpenMS::ConsensusMap":

    cdef cppclass FileDescription:
        String filename
        String label
        Size size
        UInt64 unique_id

        FileDescription()
        FileDescription(FileDescription &)  # wrap-ignore

    # for msvc++ compiler, see addons/ConsensusMap.pyx
    ctypedef Map[unsigned long int, FileDescription] FileDescriptions "OpenMS::ConsensusMap::FileDescriptions"
    ctypedef Map[unsigned long int, FileDescription].iterator FileDescriptions_iterator "OpenMS::ConsensusMap::FileDescriptions::iterator"



cdef extern from "<OpenMS/KERNEL/ConsensusMap.h>" namespace "OpenMS":

    cdef cppclass ConsensusMap(UniqueIdInterface):

        # wrap-inherits:
        #   UniqueIdInterface

        #void setFileDescriptions(Map[long unsigned int,FileDescription] fd)
        #ctypedef Map[long unsigned int, FileDescription] FileDescriptions

        ConsensusMap() nogil except +
        ConsensusMap(ConsensusMap) nogil except + # wrap-ignore

        bool operator==(ConsensusMap) nogil except +
        bool operator!=(ConsensusMap) nogil except +

        void clear(bool clear_meta_data)
        void clear()

        void updateRanges()

        libcpp_vector[ProteinIdentification] getProteinIdentifications() nogil except+
        void setProteinIdentifications(libcpp_vector[ProteinIdentification]) nogil except+

        libcpp_vector[PeptideIdentification] getUnassignedPeptideIdentifications() nogil except+
        void setUnassignedPeptideIdentifications(libcpp_vector[PeptideIdentification]) nogil except+

        libcpp_vector[DataProcessing] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[DataProcessing])   nogil except +


        libcpp_vector[ConsensusFeature].iterator begin() nogil except +    # wrap-iter-begin:__iter__(ConsensusFeature)
        libcpp_vector[ConsensusFeature].iterator end()   nogil except +    # wrap-iter-end:__iter__(ConsensusFeature)

        void applyMemberFunction(Size(* fun)())  nogil except + # wrap-ignore

        void sortByIntensity(bool reverse) nogil except +
        void sortByIntensity() nogil except +
        void sortByRT() nogil except +
        void sortByMZ() nogil except +
        void sortByPosition() nogil except +
        void sortByQuality(bool reverse) nogil except +
        void sortByQuality() nogil except +
        void sortBySize() nogil except +
        void sortByMaps() nogil except +

        int size() nogil except +

        #Map[unsigned long int,FileDescription] getFileDescriptions()
        #void setFileDescriptions(Map[long unsigned int,FileDescription] fd)
        #void setFileDescriptions(FileDescriptions &)

        FileDescriptions & getFileDescriptions()       #wrap-ignore
        #void setFileDescriptions(FileDescriptions &)       #wrap-ignore




