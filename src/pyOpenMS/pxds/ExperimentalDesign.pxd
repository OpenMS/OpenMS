from libcpp cimport *
from Types cimport *
from String cimport *
from ProteinIdentification cimport *
from ConsensusMap cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS":
    
    cdef cppclass ExperimentalDesign "OpenMS::ExperimentalDesign":
    # wrap-doc:
    #  Representation of an experimental design in OpenMS. Instances can be loaded with the ExperimentalDesignFile class

        ExperimentalDesign() except + nogil  # compiler 
        ExperimentalDesign(ExperimentalDesign &) except + nogil  # compiler

        libcpp_vector[ ExperimentalDesign_MSFileSectionEntry ] getMSFileSection() except + nogil  # TODO
        void setMSFileSection(libcpp_vector[ ExperimentalDesign_MSFileSectionEntry ] msfile_section) except + nogil  # TODO

        # Returns the Sample Section of the experimental design file
        ExperimentalDesign_SampleSection getSampleSection() except + nogil  # wrap-doc:Returns the Sample Section of the experimental design file
        void setSampleSection(ExperimentalDesign_SampleSection sample_section) except + nogil  # wrap-doc:Sets the Sample Section of the experimental design file

        # fraction index to file paths (ordered by fraction_group)
        #lib_map[unsigned int, lib_cpp[ String ] ] getFractionToMSFilesMapping() except + nogil 
     
        # return <file_path, label> to sample mapping
        #lib_map[ std::pair< String, unsigned >, unsigned] getPathLabelToSampleMapping(bool) except + nogil 

        # return <file_path, label> to fraction mapping
        #std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionMapping(bool) const;

        # return <file_path, label> to fraction_group mapping
        #std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionGroupMapping(bool) const;

        # @return the number of samples measured (= highest sample index)
        unsigned int getNumberOfSamples() except + nogil  # wrap-doc:Returns the number of samples measured (= highest sample index)

        # @return the number of fractions (= highest fraction index)
        unsigned int getNumberOfFractions() except + nogil  # wrap-doc:Returns the number of fractions (= highest fraction index)

        # @return the number of labels per file
        unsigned int getNumberOfLabels() except + nogil  # wrap-doc:Returns the number of labels per file

        # @return the number of MS files (= fractions * fraction_groups)
        unsigned int getNumberOfMSFiles() except + nogil  # wrap-doc:Returns the number of MS files (= fractions * fraction_groups)

        # @return the number of fraction_groups
        # Allows to group fraction ids and source files
        unsigned int getNumberOfFractionGroups() except + nogil  # wrap-doc:Allows to group fraction ids and source files. Return the number of fraction_groups

        # @return sample index (depends on fraction_group and label)
        unsigned int getSample(unsigned int fraction_group, unsigned int label) except + nogil  # wrap-doc:Returns sample index (depends on fraction_group and label)

        # @return whether at least one fraction_group in this experimental design is fractionated
        bool isFractionated() except + nogil  # wrap-doc:Returns whether at least one fraction_group in this experimental design is fractionated

        # return if each fraction number is associated with the same number of fraction_group
        bool sameNrOfMSFilesPerFraction() except + nogil  # wrap-doc:Returns if each fraction number is associated with the same number of fraction_group
                
# COMMENT: wrap static methods
cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":

        # Extract experimental design from consensus map
        ExperimentalDesign fromConsensusMap(ConsensusMap c) except + nogil  #wrap-attach:ExperimentalDesign

        # Extract experimental design from feature map
        ExperimentalDesign fromFeatureMap(FeatureMap f) except + nogil  #wrap-attach:ExperimentalDesign

        # Extract experimental design from identifications
        ExperimentalDesign fromIdentifications(libcpp_vector[ProteinIdentification] & proteins) except + nogil  #wrap-attach:ExperimentalDesign

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":
    
    cdef cppclass ExperimentalDesign_MSFileSectionEntry "OpenMS::ExperimentalDesign::MSFileSectionEntry":

        ExperimentalDesign_MSFileSectionEntry() except + nogil 
        ExperimentalDesign_MSFileSectionEntry(ExperimentalDesign_MSFileSectionEntry) except + nogil  #wrap-ignore

        libcpp_string path
        unsigned int fraction_group
        unsigned int fraction
        unsigned int label
        unsigned int sample
        
cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":

    cdef cppclass ExperimentalDesign_SampleSection "OpenMS::ExperimentalDesign::SampleSection":

          ExperimentalDesign_SampleSection() except + nogil 
          ExperimentalDesign_SampleSection(ExperimentalDesign_SampleSection) except + nogil 

          # ExperimentalDesign_SampleSection(const libcpp_vector[ libcpp_vector[ String ] ] & content,
          #                                 libcpp_map[ unsigned, size_t ] sample_to_rowindex, 
          #                                 libcpp_map[ String, size_t ] columnname_to_columnindex) except + nogil 

          # Get set of all samples that are present in the sample section
          libcpp_set[ String ] getSamples() except + nogil  # wrap-doc:Returns a set of all samples that are present in the sample section

          # Get set of all factors (column names) that were defined for the sample section
          libcpp_set[ String ] getFactors() except + nogil  # wrap-doc:Returns a set of all factors (column names) that were defined for the sample section

          # Checks whether sample section has row for a sample number
          bool hasSample(unsigned int sample) except + nogil  # wrap-doc:Checks whether sample section has row for a sample number

          # Checks whether Sample Section has a specific factor (i.e. column name)
          bool hasFactor(String& factor) except + nogil  # wrap-doc:Checks whether Sample Section has a specific factor (i.e. column name)

          # Returns value of factor for given sample and factor name
          String getFactorValue(unsigned int sample, String &factor) except + nogil  # wrap-doc:Returns value of factor for given sample and factor name
                   

