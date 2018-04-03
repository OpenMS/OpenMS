// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>
#include <map>
#include <set>

namespace OpenMS
{
  /**
  @brief Representation of the Experimental Design in OpenMS. Instances can be loaded via
   the ExperimentalDesignFile class.

  @ingroup Format
  */
  class OPENMS_DLLAPI ExperimentalDesign
  {
  public:

    class OPENMS_DLLAPI RunRow
    {
    public:
      RunRow() = default;
      unsigned run = 1; ///< run index (before prefractionation)
      unsigned fraction = 1; ///< fraction 1..m, mandatory, 1 if not set
      std::string path = "UNKNOWN_FILE"; ///< file name, mandatory
      unsigned channel = 1;  ///< if and how many multiplexed channels are in a file
      unsigned sample = 1;  ///< allows grouping by sample
    };

    class OPENMS_DLLAPI SampleSection
    {
    public:

      SampleSection() = default;

      SampleSection(
        std::vector< std::vector < String > > _content,
        std::map< unsigned, Size > _sample_to_rowindex,
        std::map< String, Size > _columnname_to_columnindex
      ) : content_(_content),
          sample_to_rowindex_(_sample_to_rowindex),
          columnname_to_columnindex_(_columnname_to_columnindex) {}

      // Get set of all samples that are present in the sample section
      std::set< unsigned > getSamples() const;

      // Get set of all factors (column names) that were defined for the sample section
      std::set< String > getFactors() const;

      // Checks whether sample section has row for a sample number
      bool hasSample(unsigned sample) const;

      // Checks whether Sample Section has a specific factor (i.e. column name)
      bool hasFactor(const String &factor) const;

      // Returns value of factor for given sample and factor name
      String getFactorValue(unsigned sample, const String &factor);

    private:

      // The entries of the Sample Section, filled while parsing
      // the Experimental Design File
      std::vector< std::vector < String > > content_;

      // Maps the Sample Entry to the row where the sample
      // appears in the Sample section
      std::map< unsigned, Size > sample_to_rowindex_;

      // Maps the column name of the SampleSection to the
      // Index of the column
      std::map< String, Size > columnname_to_columnindex_;
    };

    using RunRows = std::vector<RunRow>;

    // Experimental Design c'tors
    ExperimentalDesign() = default;

    ExperimentalDesign(
      RunRows run_rows, SampleSection sample_section) : run_section_(run_rows), sample_section_(sample_section)
    {
      sort_();
      checkValidRunSection_();
    }

    const RunRows& getRunSection() const;

    void setRunSection(const RunRows& run_section);

    // Returns the Sample Section of the experimental design file
    const ExperimentalDesign::SampleSection& getSampleSection() const;

    void setSampleSection(const SampleSection& sample_section);

    // Gets vector of Filenames that appears in the run section, optionally trims to basename
    std::vector< String > getFileNames(bool basename) const;

    // Returns vector of channels of the run section
    std::vector<unsigned> getChannels() const;

    std::vector<unsigned> getFractions() const;

    /// return fraction index to file paths (ordered by run id)
    std::map<unsigned int, std::vector<String> > getFractionToMSFilesMapping() const;

   /*
    *   The (Path, Channel) tuples in the experimental design have to be unique, so we can map them
    *   uniquely to the sample number, fraction number, and run number
    */
    /// return <file_path, channel> to sample mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToSampleMapping(bool) const;

    /// return <file_path, channel> to fraction mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToFractionMapping(bool) const;

    /// return <file_path, channel> to run mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToRunMapping(bool) const;

    // @return the number of samples measured (= highest sample index)
    unsigned getNumberOfSamples() const;

    // @return the number of fractions (= highest fraction index)
    unsigned getNumberOfFractions() const;

    // @return the number of channels per file
    unsigned getNumberOfChannels() const;

    // @return the number of MS files (= fractions * runs)
    unsigned getNumberOfMSFiles() const;

    // @return the number of runs (before fractionation)
    // Allows to group fraction ids and source files
    unsigned getNumberOfPrefractionationRuns() const;

    // @return sample index (depends on run and channel)
    unsigned getSample(unsigned run, unsigned channel = 1);

    /// @return whether at least one run in this experimental design is fractionated
    bool isFractionated() const;

    /// return if each fraction number is associated with the same number of runs
    bool sameNrOfMSFilesPerFraction() const;

    /// Extract experimental design from consensus map
    static ExperimentalDesign fromConsensusMap(const ConsensusMap& c);

    /// Extract experimental design from feature map
    static ExperimentalDesign fromFeatureMap(const FeatureMap& f);

    /// Extract experimental design from identifications
    static ExperimentalDesign fromIdentifications(const std::vector<ProteinIdentification> & proteins);

    private:

      /// Generic Mapper (Path, Channel) -> f(row)
      std::map< std::pair< String, unsigned >, unsigned> pathChannelMapper(
          bool,
          unsigned (*f)(const ExperimentalDesign::RunRow&)) const;

      // sort to obtain the default order
      void sort_();

      template<typename T>
      static void errorIfAlreadyExists(std::set<T> &container, T &item, const String &message);

      void checkValidRunSection_();

      RunRows run_section_;
      SampleSection sample_section_;
  };
}
