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

    cdef cppclass ConsensusFeature(UniqueIdInterface, BaseFeature):
        # wrap-inherits:
        #   UniqueIdInterface
        #   BaseFeature
        #
        # wrap-doc:
        #  A consensus feature spanning multiple LC-MS/MS experiments.
        #  
        #  A ConsensusFeature represents analytes that have been
        #  quantified across multiple LC-MS/MS experiments. Each analyte in a
        #  ConsensusFeature is linked to its original LC-MS/MS run through a
        #  unique identifier.
        #  
        #  Get access to the underlying features through getFeatureList()

        ConsensusFeature() except + nogil 
        ConsensusFeature(ConsensusFeature &) except + nogil 
        ConsensusFeature(UInt64, Peak2D, UInt64) except + nogil 
        ConsensusFeature(UInt64, BaseFeature) except + nogil 
        ConsensusFeature(UInt64, ConsensusFeature) except + nogil 

        void computeConsensus()    except + nogil  # wrap-doc:Computes and updates the consensus position, intensity, and charge
        void computeMonoisotopicConsensus()    except + nogil  # wrap-doc:Computes and updates the consensus position, intensity, and charge
        void computeDechargeConsensus(FeatureMap, bool)    except + nogil  # wrap-doc:Computes the uncharged parent RT & mass, assuming the handles are charge variants

        void insert(UInt64 map_idx, Peak2D, UInt64 element_idx) except + nogil 
        void insert(UInt64 map_idx, BaseFeature) except + nogil 
        void insert(UInt64 map_idx, ConsensusFeature) except + nogil 

        libcpp_vector[FeatureHandle] getFeatureList() except + nogil 

        Size size() except + nogil 

        bool operator==(ConsensusFeature) except + nogil 
        bool operator!=(ConsensusFeature) except + nogil 

        void addRatio(Ratio r) except + nogil  # wrap-doc:Connects a ratio to the ConsensusFeature.
        void setRatios(libcpp_vector[Ratio] rs) except + nogil  # wrap-doc:Connects the ratios to the ConsensusFeature.
        libcpp_vector[Ratio] getRatios() except + nogil  # wrap-doc:Get the ratio vector.

        void clear() except + nogil 
        bool empty() except + nogil 

        # # Returns the position range of the contained elements
        # DRange2 getPositionRange() except + nogil 
        # # Returns the intensity range of the contained elements
        # DRange1 getIntensityRange() except + nogil 

cdef extern from "<OpenMS/KERNEL/ConsensusFeature.h>" namespace "OpenMS::ConsensusFeature":

    # slim struct to feed the need for systematically storing of ratios .
    cdef cppclass Ratio:

      Ratio() except + nogil 
      Ratio(Ratio rhs) except + nogil 

      double ratio_value_
      String denominator_ref_
      String numerator_ref_
      libcpp_vector[String] description_
