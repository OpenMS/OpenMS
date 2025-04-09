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
        #  UniqueIdInterface
        #  DocumentIdentifier
        #  RangeManagerRtMzInt
        #  MetaInfoInterface
        #
        # wrap-instances:
        #  FeatureMap := FeatureMap
        #
        # wrap-doc:
        #  A container for features.
        #  
        #  A feature map is a container holding features, which represent
        #  chemical entities (peptides, proteins, small molecules etc.) found
        #  in an LC-MS/MS experiment.
        #  
        #  This class supports direct iteration in Python.

        FeatureMap() except + nogil 
        FeatureMap(FeatureMap &) except + nogil 

        bool operator==(FeatureMap) except + nogil 
        bool operator!=(FeatureMap) except + nogil 

        int size()  except + nogil 
        Feature & operator[](size_t)      except + nogil  #wrap-upper-limit:size()
        void push_back(Feature spec) except + nogil 
        void push_back(MRMFeature spec) except + nogil 

        void sortByIntensity() except + nogil  # wrap-doc:Sorts the peaks according to ascending intensity
        void sortByIntensity(bool reverse) except + nogil  # wrap-doc:Sorts the peaks according to ascending intensity. Order is reversed if argument is `true` ( reverse = true )
        void sortByPosition() except + nogil  # wrap-doc:Sorts features by position. Lexicographical comparison (first RT then m/z) is done
        void sortByRT() except + nogil  # wrap-doc:Sorts features by RT position
        void sortByMZ() except + nogil  # wrap-doc:Sorts features by m/z position
        void sortByOverallQuality() except + nogil  # wrap-doc:Sorts features by ascending overall quality. Order is reversed if argument is `true` ( reverse = true )

        void swap(FeatureMap &) except + nogil  
        void swapFeaturesOnly(FeatureMap swapfrom) except + nogil  # wrap-doc:Swaps the feature content (plus its range information) of this map 
        void clear() except + nogil  # wrap-doc:Clears all data and meta data
        void clear(bool clear_meta_data) except + nogil  # wrap-doc:Clears all data and meta data. If 'true' is passed as an argument, all meta data is cleared in addition to the data

        FeatureMap operator+(FeatureMap) except + nogil 
        FeatureMap iadd(FeatureMap) except + nogil  # wrap-as:operator+=

        void updateRanges() except + nogil  # TODO

        libcpp_vector[ProteinIdentification] getProteinIdentifications() except + nogil
        void setProteinIdentifications(libcpp_vector[ProteinIdentification]) except + nogil # wrap-doc:Sets the protein identifications

        libcpp_vector[PeptideIdentification] getUnassignedPeptideIdentifications() except + nogil
        void setUnassignedPeptideIdentifications(libcpp_vector[PeptideIdentification]) except + nogil # wrap-doc:Sets the unassigned peptide identifications

        Size applyMemberFunction(Size(* fun)()) except + nogil # wrap-ignore

        libcpp_vector[DataProcessing] getDataProcessing() except + nogil 
        void setDataProcessing(libcpp_vector[DataProcessing])   except + nogil  # wrap-doc:Sets the description of the applied data processing

        void setPrimaryMSRunPath(StringList& s) except + nogil  # wrap-doc:Sets the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
        void setPrimaryMSRunPath(StringList& s, MSExperiment& e) except + nogil  # wrap-doc:Sets the file path to the primary MS run using the mzML annotated in the MSExperiment argument `e`
        void getPrimaryMSRunPath(StringList& toFill) except + nogil  # wrap-doc:Returns the file path to the first MS run

        libcpp_vector[Feature].iterator begin() except + nogil     # wrap-iter-begin:__iter__(Feature)
        libcpp_vector[Feature].iterator end()   except + nogil     # wrap-iter-end:__iter__(Feature)

