from libcpp cimport *
from Types cimport *
from String cimport *
from ProteinIdentification cimport *
from ConsensusMap cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS":
    
    cdef cppclass ExperimentalDesign "OpenMS::ExperimentalDesign":
        ExperimentalDesign() nogil except +
        ExperimentalDesign(ExperimentalDesign) nogil except + #wrap-ignore

        libcpp_vector[ ExperimentalDesign_MSFileSectionEntry ] getMSFileSection() nogil except +
        void setMSFileSection(libcpp_vector[ ExperimentalDesign_MSFileSectionEntry ] msfile_section) nogil except +

        # Returns the Sample Section of the experimental design file
        ExperimentalDesign_SampleSection getSampleSection() nogil except +
        void setSampleSection(ExperimentalDesign_SampleSection sample_section) nogil except +

        # fraction index to file paths (ordered by fraction_group)
        #lib_map[unsigned int, lib_cpp[ String ] ] getFractionToMSFilesMapping() nogil except +
     
        # return <file_path, label> to sample mapping
        #lib_map[ std::pair< String, unsigned >, unsigned] getPathLabelToSampleMapping(bool) nogil except +

        # return <file_path, label> to fraction mapping
        #std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionMapping(bool) const;

        # return <file_path, label> to fraction_group mapping
        #std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionGroupMapping(bool) const;

        # @return the number of samples measured (= highest sample index)
        unsigned int getNumberOfSamples() nogil except +

        # @return the number of fractions (= highest fraction index)
        unsigned int getNumberOfFractions() nogil except +

        # @return the number of labels per file
        unsigned int getNumberOfLabels() nogil except +

        # @return the number of MS files (= fractions * fraction_groups)
        unsigned int getNumberOfMSFiles() nogil except +

        # @return the number of fraction_groups
        # Allows to group fraction ids and source files
        unsigned int getNumberOfFractionGroups() nogil except +

        # @return sample index (depends on fraction_group and label)
        unsigned int getSample(unsigned int fraction_group, unsigned int label) nogil except +

        # @return whether at least one fraction_group in this experimental design is fractionated
        bool isFractionated() nogil except +

        # return if each fraction number is associated with the same number of fraction_group
        bool sameNrOfMSFilesPerFraction() nogil except +
                
# COMMENT: wrap static methods
cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":

        # Extract experimental design from consensus map
        ExperimentalDesign fromConsensusMap(ConsensusMap c) nogil except + #wrap-attach:ExperimentalDesign

        # Extract experimental design from feature map
        ExperimentalDesign fromFeatureMap(FeatureMap f) nogil except + #wrap-attach:ExperimentalDesign

        # Extract experimental design from identifications
        ExperimentalDesign fromIdentifications(libcpp_vector[ProteinIdentification] & proteins) nogil except + #wrap-attach:ExperimentalDesign

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":
    
    cdef cppclass ExperimentalDesign_MSFileSectionEntry "OpenMS::ExperimentalDesign::MSFileSectionEntry":

        ExperimentalDesign_MSFileSectionEntry() nogil except +
        ExperimentalDesign_MSFileSectionEntry(ExperimentalDesign_MSFileSectionEntry) nogil except + #wrap-ignore

        libcpp_string path
        unsigned int fraction_group
        unsigned int fraction
        unsigned int label
        unsigned int sample
        
cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":

    cdef cppclass ExperimentalDesign_SampleSection "OpenMS::ExperimentalDesign::SampleSection":

          ExperimentalDesign_SampleSection() nogil except +
          ExperimentalDesign_SampleSection(ExperimentalDesign_SampleSection) nogil except +

          # ExperimentalDesign_SampleSection(const libcpp_vector[ libcpp_vector[ String ] ] & content,
          #                                  libcpp_map[ unsigned, size_t ] sample_to_rowindex, 
          #                                  libcpp_map[ String, size_t ] columnname_to_columnindex) nogil except +

          # Get set of all samples that are present in the sample section
          libcpp_set[ unsigned int ] getSamples() nogil except +

          # Get set of all factors (column names) that were defined for the sample section
          libcpp_set[ String ] getFactors() nogil except +

          # Checks whether sample section has row for a sample number
          bool hasSample(unsigned int sample) nogil except +

          # Checks whether Sample Section has a specific factor (i.e. column name)
          bool hasFactor(String& factor) nogil except +

          # Returns value of factor for given sample and factor name
          String getFactorValue(unsigned int sample, String &factor) nogil except +
                   

