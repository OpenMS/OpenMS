// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_KERNEL_EXPERIMENTALDESIGN_H
#define OPENMS_KERNEL_EXPERIMENTALDESIGN_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <tuple>


namespace OpenMS
{
  /**
  @brief A TSV and user friendly representation of the experimental design.
  Used for loading and storing the experimental design in OpenMS.

  @ingroup Metadata
  */
  class OPENMS_DLLAPI ExperimentalDesign
  {
  public:
    ExperimentalDesign() {};

    /// 1) Mandatory section with run-level information of the experimental design.
    ///    Required to process fractionated data.
/*
 * Run Section Format:
   Format: Single header line
         Run:           Run index (prior fractionation) used to group fractions and source files.
                        Note: For label-free this has same cardinality as sample.
                              For multiplexed experiments, these might differ as multiple samples can be measured in single files
         Fraction:      1st, 2nd, .., fraction. Note: All runs must have the same number of fractions.
         Path:          Path to mzML files
         Channel:       Channel in MS file:
                          label-free: always 1
                          TMT6Plex: 1..6
                          SILAC with light and heavy: 1..2
         Sample:        Index of sample measured in the specified channel X, in fraction Y of run Z

	Run	Fraction	Path (Spectra File)	Channel		Sample (Condition)
	1	1		SPECTRAFILE_F1_TR1.mzML	1		1
	1	2		SPECTRAFILE_F2_TR1.mzML	1		1
	1	3		SPECTRAFILE_F3_TR1.mzML	1		1
	1	1		SPECTRAFILE_F1_TR1.mzML	2		2
	1	2		SPECTRAFILE_F2_TR1.mzML	2		2
	1	3		SPECTRAFILE_F3_TR1.mzML	2		2
	1	1		SPECTRAFILE_F1_TR1.mzML	3		3
	1	2		SPECTRAFILE_F2_TR1.mzML	3		3
	1	3		SPECTRAFILE_F3_TR1.mzML	3		3
	1	1		SPECTRAFILE_F1_TR1.mzML	4		4
	1	2		SPECTRAFILE_F2_TR1.mzML	4		4
	1	3		SPECTRAFILE_F3_TR1.mzML	4		4
	2	1		SPECTRAFILE_F1_TR2.mzML	1		5
	2	2		SPECTRAFILE_F2_TR2.mzML	1		5
	2	3		SPECTRAFILE_F3_TR2.mzML	1		5
	2	1		SPECTRAFILE_F1_TR2.mzML	2		6
	2	2		SPECTRAFILE_F2_TR2.mzML	2		6
	2	3		SPECTRAFILE_F3_TR2.mzML	2		6
	2	1		SPECTRAFILE_F1_TR2.mzML	3		7
	2	2		SPECTRAFILE_F2_TR2.mzML	3		7
	2	3		SPECTRAFILE_F3_TR2.mzML	3		7
	2	1		SPECTRAFILE_F1_TR2.mzML	4		8
	2	2		SPECTRAFILE_F2_TR2.mzML	4		8
	2	3		SPECTRAFILE_F3_TR2.mzML	4		8
*/
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

      // Number of columns of the Sample Section
      Size n_columns_;

      // The entries of the Sample Section, filled while parsing
      // the Experimental Design File
      std::vector< std::vector < String > > content_;

      // Maps the Sample Entry to the row where the sample
      // appears in the Sample section
      std::map< unsigned, Size > sample_to_rowindex_;

      // Maps the column name of the SampleSection to the
      // Index of the column
      std::map< String, Size > columnname_to_columnindex_;

      // Get Sample Attributes
      void getSamples(std::set< unsigned > &samples) const
      {
        samples.clear();
        for (auto it = sample_to_rowindex_.begin();
             it != sample_to_rowindex_.end(); ++it)
        {
          samples.insert(it->first);
        }
      }
      // Get Sample Attributes
      void getFactors(std::set< String > &factors) const
      {
        factors.clear();
        for (auto it = columnname_to_columnindex_.begin();
             it != columnname_to_columnindex_.end(); ++it)
        {
          factors.insert(it->first);
        }
      }

      bool hasSample(const unsigned sample) const
      {
        return sample_to_rowindex_.find(sample) != sample_to_rowindex_.end();
      }

      bool hasFactor(const String &factor) const
      {
        return columnname_to_columnindex_.find(factor) != columnname_to_columnindex_.end();
      }

      String getFactorValue(const unsigned sample, const String &factor)
      {
        if (! hasSample(sample))
        {
          throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Sample " + String(sample) + " is not present in the Experimental Design");
        }
        if (! hasFactor(factor))
        {
          throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Factor " + factor + " is not present in the Experimental Design");
        }
        StringList sample_row = content_[sample_to_rowindex_[sample]];
        const Size col_index = columnname_to_columnindex_[factor];
        return sample_row[col_index];
      }
    };


    using RunRows = std::vector<RunRow>;


    const RunRows& getRunSection() const
    {
      return run_section_;
    }

    void setRunSection(const RunRows& run_section)
    {
      run_section_ = run_section;
      sort_();
      checkValidRunSection_();
    }


    // Returns the Sample Section of the experimental design file
    const SampleSection& getSampleSection() const
    {
      return sample_section_;
    }

    void getFileNames(std::vector< String > &filenames, const bool basename) const
    {
      filenames.clear();
      for (const RunRow &row : run_section_)
      {
        const String path = String(row.path);
        filenames.push_back(basename ? path : File::basename(path));
      }
    }

    void getChannels(std::vector<unsigned> &channels) const
    {
      channels.clear();
      for (const RunRow &row : run_section_)
      {
        channels.push_back(row.channel);
      }
    }


    bool isLabelFree() const
    {
      std::vector<unsigned> channels;
      getChannels(channels);
      std::unique(channels.begin(), channels.end());

      // At most one channel is label-free
      return channels.size() < 2;
    }

    bool hasFractions() const
    {
      std::set<unsigned> fractions;
      for (RunRow const & r : run_section_)
      {
        fractions.insert(r.fraction);
      }
      return fractions.size() > 1;
    }

    /// return fraction index to file paths (ordered by run id)
    std::map<unsigned int, std::vector<String> > getFractionToMSFilesMapping() const;

   /*
    *   The (Path, Channel) tuples in the experimental design have to be unique, so we can map them
    *   uniquely to the sample number, fraction number, and run number
    */
    /// return <file_path, channel> to sample mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToSampleMapping(const bool) const;

    /// return <file_path, channel> to fraction mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToFractionMapping(const bool) const;

    /// return <file_path, channel> to run mapping
    std::map< std::pair< String, unsigned >, unsigned> getPathChannelToRunMapping(const bool) const;


    // @return the number of samples measured (= highest sample index)
    unsigned getNumberOfSamples() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
        [](const RunRow& f1, const RunRow& f2)
        {
          return f1.sample < f2.sample;
        })->sample;
    }

    // @return the number of fractions (= highest fraction index)
    unsigned getNumberOfFractions() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
        [](const RunRow& f1, const RunRow& f2)
        {
          return f1.fraction < f2.fraction;
        })->fraction;
    }

    // @return the number of channels per file
    unsigned getNumberOfChannels() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
        [](const RunRow& f1, const RunRow& f2)
        {
          return f1.fraction < f2.fraction;
        })->channel;
    }

    // @return the number of MS files (= fractions * runs)
    unsigned getNumberOfMSFiles() const
    {
      std::set<std::string> unique_paths;
      for (auto const & r : run_section_) { unique_paths.insert(r.path); }
      return unique_paths.size();
    }

    // @return the number of runs (before fractionation)
    // Allows to group fraction ids and source files
    unsigned getNumberOfPrefractionationRuns() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
        [](const RunRow& f1, const RunRow& f2)
        {
          return f1.run < f2.run;
        })->run;
    }

    // @return sample index (depends on run and channel)
    unsigned getSample(unsigned run, unsigned channel = 1)
    {
      return std::find_if(run_section_.begin(), run_section_.end(),
        [&run, &channel](const RunRow& r)
        {
          return r.run == run && r.channel == channel;
        })->sample;
    }

    /// return if each fraction number is associated with the same number of runs
    bool sameNrOfMSFilesPerFraction() const;

    /// Loads an experimental design from a tabular separated file
    static ExperimentalDesign load(const String & tsv_file);

    /// Extract experimental design from consensus map
    static ExperimentalDesign fromConsensusMap(const ConsensusMap& c);

    /// Extract experimental design from feature map
    static ExperimentalDesign fromFeatureMap(const FeatureMap& f);

    /// Extract experimental design from identifications
    static ExperimentalDesign fromIdentifications(const std::vector<ProteinIdentification> & proteins);

    private:

      /// Generic Mapper (Path, Channel) -> f(row)
      std::map< std::pair< String, unsigned >, unsigned> pathChannelMapper(
          const bool,
          unsigned (*f)(const ExperimentalDesign::RunRow&)) const;

      // sort to obtain the default order
      void sort_()
      {
        std::sort(run_section_.begin(), run_section_.end(),
        [](const RunRow& a, const RunRow& b)
        {
          return std::tie(a.run, a.fraction, a.channel, a.sample, a.path) <
            std::tie(b.run, b.fraction, b.channel, b.sample, b.path);
        });
      }

      template<typename T>
      void errorIfAlreadyExists(std::set<T> &container, T &item, const String &message)
      {
        if (container.find(item) != container.end())
        {
          throw Exception::MissingInformation(
                  __FILE__,
                  __LINE__,
                  OPENMS_PRETTY_FUNCTION, message);
        }
        container.insert(item);
      }

      void checkValidRunSection_()
      {
        if (getNumberOfMSFiles() == 0)
        {
          throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "No MS files provided.");
        }

        std::set< std::tuple< unsigned, unsigned, unsigned > > run_fraction_sample_set;
        std::set< std::tuple< unsigned, unsigned, unsigned > > run_fraction_channel_set;
        std::set< std::tuple< std::string, unsigned > > path_channel_set;
        std::map< std::tuple< unsigned, unsigned >, std::set< unsigned > > run_channel_to_sample;

        for (const RunRow &row : run_section_)
        {
          // Fail if sample section does not contain the run
          if (sample_section_.hasSample(row.sample) == false)
          {
            throw Exception::MissingInformation(
                    __FILE__,
                    __LINE__,
                    OPENMS_PRETTY_FUNCTION,
                    "Sample Section does not contain sample for run " + String(row.run));
          }

          // RUN_FRACTION_SAMPLE TUPLE
          std::tuple<unsigned, unsigned, unsigned> run_fraction_sample = std::make_tuple(row.run, row.fraction, row.sample);
          errorIfAlreadyExists(
                  run_fraction_sample_set,
                  run_fraction_sample,
                  "(Run, Fraction, Sample) combination can only appear once");

          // RUN_FRACTION_CHANNEL TUPLE
          std::tuple<unsigned, unsigned, unsigned> run_fraction_channel = std::make_tuple(row.run, row.fraction, row.channel);
          errorIfAlreadyExists(
                  run_fraction_channel_set,
                  run_fraction_channel,
                  "(Run, Fraction, Channel) combination can only appear once");


          // PATH_CHANNEL_TUPLE
          std::tuple<std::string, unsigned> path_channel = std::make_tuple(row.path, row.channel);
          errorIfAlreadyExists(
                  path_channel_set,
                  path_channel,
                  "(Path, Channel) combination can only appear once");

          // RUN_CHANNEL TUPLE
          std::tuple<unsigned, unsigned> run_channel = std::make_tuple(row.run, row.channel);
          run_channel_to_sample[run_channel].insert(row.sample);

          if (run_channel_to_sample[run_channel].size() > 1)
          {
            throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Multiple Samples encountered for the same Run and the same Channel");
          }
        }
      }

      RunRows run_section_;
      SampleSection sample_section_;
  };
}


#endif // header guard
