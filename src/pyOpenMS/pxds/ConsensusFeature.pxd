from libcpp cimport bool
from Types cimport *
from BaseFeature cimport *
from Peak2D cimport *
from UniqueIdInterface cimport *
from FeatureMap cimport *
from BaseFeature cimport *
from FeatureHandle cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/KERNEL/ConsensusFeature.h>" namespace "OpenMS":

    # do not wrap BaseFeature, due to overloaded base methods
    # -> see Precursor.pxd

    cdef cppclass ConsensusFeature(UniqueIdInterface,Peak2D):
        # wrap-inherits:
        #    UniqueIdInterface
        #    Peak2D
        #
        # wrap-doc:
        #   A consensus feature spanning multiple LC-MS/MS experiments.
        #   -----
        #   A ConsensusFeature represents analytes that have been
        #   quantified across multiple LC-MS/MS experiments. Each analyte in a
        #   ConsensusFeature is linked to its original LC-MS/MS run through a
        #   unique identifier.
        #   -----
        #   Get access to the underlying features through getFeatureList()

        ConsensusFeature() nogil except +
        ConsensusFeature(ConsensusFeature &) nogil except +
        ConsensusFeature(UInt64, Peak2D, UInt64) nogil except +
        ConsensusFeature(UInt64, BaseFeature) nogil except +
        ConsensusFeature(UInt64, ConsensusFeature) nogil except +

        void computeConsensus()    nogil except +
        void computeMonoisotopicConsensus()    nogil except +
        void computeDechargeConsensus(FeatureMap, bool)    nogil except +

        void insert(UInt64 map_idx, Peak2D, UInt64 element_idx) nogil except +
        void insert(UInt64 map_idx, BaseFeature) nogil except +
        void insert(UInt64 map_idx, ConsensusFeature) nogil except +

        float getQuality()  nogil except +
        void setQuality(float q) nogil except +

        float getWidth() nogil except +
        void setWidth(float q) nogil except +

        Int getCharge() nogil except +
        void setCharge(Int q) nogil except +

        libcpp_vector[FeatureHandle] getFeatureList() nogil except +

        Size size() nogil except +

        # returns a mutable reference to the PeptideIdentification vector
        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except +
        # sets the PeptideIdentification vector
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification] & peptides) nogil except +

        bool operator==(ConsensusFeature) nogil except +
        bool operator!=(ConsensusFeature) nogil except +

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
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

        void clear() nogil except +
        bool empty() nogil except +

        # # Returns the position range of the contained elements
        # DRange2 getPositionRange() nogil except +
        # # Returns the intensity range of the contained elements
        # DRange1 getIntensityRange() nogil except +

cdef extern from "<OpenMS/KERNEL/ConsensusFeature.h>" namespace "OpenMS::ConsensusFeature":

    # slim struct to feed the need for systematically storing of ratios ( @see MSQuantifications ).
    cdef cppclass Ratio:

      Ratio() nogil except +
      Ratio(Ratio rhs) nogil except +

      double ratio_value_
      String denominator_ref_
      String numerator_ref_
      libcpp_vector[String] description_
