from libcpp cimport *
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS":
    
    cdef cppclass ExperimentalDesign "OpenMS::ExperimentalDesign":
        ExperimentalDesign() nogil except +
        ExperimentalDesign(ExperimentalDesign) nogil except + #wrap-ignore

        RunRows& getRunSection() nogil except +
        setRunSection() nogil except +

        # Returns the Sample Section of the experimental design file
        ExperimentalDesign_SampleSection getSampleSection() nogil except +
        setSampleSection(const SampleSection& sample_section) nogil except +

        # Gets vector of Filenames that appears in the run section, optionally trims to basename
        lib_cpp[ String ] getFileNames(bool basename) nogil except +

        # Returns vector of channels of the run section
        lib_cpp[ unsigned ] getChannels() nogil except +

        lib_cpp[ unsigned ] getFractions() nogil except +

        # return fraction index to file paths (ordered by run id)
        #lib_map[unsigned int, lib_cpp[ String ] ] getFractionToMSFilesMapping() nogil except +
     
        # return <file_path, channel> to sample mapping
        #lib_map[ std::pair< String, unsigned >, unsigned] getPathChannelToSampleMapping(bool) nogil except +

        # return <file_path, channel> to fraction mapping
        #std::map< std::pair< String, unsigned >, unsigned> getPathChannelToFractionMapping(bool) const;

        # return <file_path, channel> to run mapping
        #std::map< std::pair< String, unsigned >, unsigned> getPathChannelToRunMapping(bool) const;

        # @return the number of samples measured (= highest sample index)
        unsigned getNumberOfSamples() nogil except +

        # @return the number of fractions (= highest fraction index)
        unsigned getNumberOfFractions() nogil except +

        # @return the number of channels per file
        unsigned getNumberOfChannels() nogil except +

        # @return the number of MS files (= fractions * runs)
        unsigned getNumberOfMSFiles() nogil except +

       # @return the number of runs (before fractionation)
       # Allows to group fraction ids and source files
       unsigned getNumberOfPrefractionationRuns() nogil except +

       # @return sample index (depends on run and channel)
       unsigned getSample(unsigned run, unsigned channel) nogil except +

       # @return whether at least one run in this experimental design is fractionated
       bool isFractionated() nogil except +

       # return if each fraction number is associated with the same number of runs
       bool sameNrOfMSFilesPerFraction() nogil except +

       # Extract experimental design from consensus map
       # static ExperimentalDesign fromConsensusMap(const ConsensusMap& c);

       # Extract experimental design from feature map
       # static ExperimentalDesign fromFeatureMap(const FeatureMap& f);

       # Extract experimental design from identifications
       # static ExperimentalDesign fromIdentifications(const std::vector<ProteinIdentification> & proteins);
                
        
cdef extern from "<OpenMS/METADATA/ExperimentalDesign.h>" namespace "OpenMS::ExperimentalDesign":
    
    cdef cppclass ExperimentalDesign_RunRow "OpenMS::ExperimentalDesign::RunRow":
        ExperimentalDesign_RunRow() nogil except +
        ExperimentalDesign_RunRow(ExperimentalDesign_RunRow) nogil except + #wrap-ignore
        # libcpp_string file
        # unsigned fraction
        # unsigned technical_replicate        
        
cdef cppclass ExperimentalDesign_SampleSection "OpenMS::ExperimentalDesign::SampleSection":
      ExperimentalDesign_SampleSection() nogil except +

      // Get set of all samples that are present in the sample section
      libcpp_set[ unsigned ] getSamples() nogil except +

      // Get set of all factors (column names) that were defined for the sample section
      libcpp_set[ unsigned ] getFactors() nogil except +

      // Checks whether sample section has row for a sample number
      bool hasSample(unsigned sample) nogil except +

      // Checks whether Sample Section has a specific factor (i.e. column name)
      bool hasFactor(const String& factor) nogil except +

      // Returns value of factor for given sample and factor name
      String getFactorValue(unsigned sample, const String &factor) nogil except +
               
