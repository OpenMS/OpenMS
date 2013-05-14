from libcpp cimport bool
from Types cimport *
from BaseFeature cimport *
from Peak2D cimport *
from RichPeak2D cimport *
from UniqueIdInterface cimport *
from FeatureMap cimport *
from BaseFeature cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/KERNEL/ConsensusFeature.h>" namespace "OpenMS":

    # do not wrap BaseFeature, due to overloaded base methods
    # -> see Precursor.pxd

    cdef cppclass ConsensusFeature(UniqueIdInterface,Peak2D):
        # wrap-inherits:
        #    UniqueIdInterface
        #    Peak2D

        ConsensusFeature() nogil except +
        ConsensusFeature(ConsensusFeature) nogil except +  #wrap-ignore
        ConsensusFeature(UInt64, Peak2D, UInt64) nogil except +
        ConsensusFeature(UInt64, BaseFeature) nogil except +
        ConsensusFeature(UInt64, ConsensusFeature) nogil except +

        void computeConsensus()    nogil except +
        void computeMonoisotopicConsensus()    nogil except +
        void computeDechargeConsensus(FeatureMap[Feature], bool)    nogil except +

        void insert(UInt64, Peak2D, UInt64) nogil except +
        void insert(UInt64, BaseFeature) nogil except +
        void insert(UInt64, ConsensusFeature) nogil except +

        Real getQuality()  nogil except +
        void setQuality(Real q) nogil except +

        Real getWidth() nogil except +
        void setWidth(Real q) nogil except +

        Int getCharge() nogil except +
        void setCharge(Int q) nogil except +

        Size size() nogil except +

        # returns a mutable reference to the PeptideIdentification vector
        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except +
        # sets the PeptideIdentification vector
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification] & peptides) nogil except +

        bool operator==(ConsensusFeature) nogil except +
        bool operator!=(ConsensusFeature) nogil except +

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +
        void clearMetaInfo() nogil except +

        void addRatio(Ratio r) nogil except +
        void setRatios(libcpp_vector[Ratio] rs) nogil except +
        libcpp_vector[Ratio] getRatios() nogil except +

        # # Returns the position range of the contained elements
        # DRange2 getPositionRange() nogil except +
        # # Returns the intensity range of the contained elements
        # DRange1 getIntensityRange() nogil except +

cdef extern from "<OpenMS/KERNEL/ConsensusFeature.h>" namespace "OpenMS::ConsensusFeature":

    # slim struct to feed the need for systematically storing of ratios ( @see MSQuantifications ).
    cdef cppclass Ratio:

      Ratio()
      Ratio(Ratio rhs)

      DoubleReal ratio_value_
      String denominator_ref_
      String numerator_ref_
      libcpp_vector[String] description_


