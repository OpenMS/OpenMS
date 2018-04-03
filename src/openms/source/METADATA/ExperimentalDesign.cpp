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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <QtCore/QString>
#include <QtCore/QFileInfo>

#include <iostream>

using namespace std;

namespace OpenMS
{
    using RunRows = std::vector<ExperimentalDesign::RunRow>;


    ExperimentalDesign ExperimentalDesign::fromConsensusMap(const ConsensusMap &cm)
    {
      ExperimentalDesign experimental_design;
      // path of the original MS run (mzML / raw file)
      StringList ms_run_paths;
      cm.getPrimaryMSRunPath(ms_run_paths);

      // no fractionation -> as many runs as samples
      // each consensus element corresponds to one sample abundance
      size_t sample(1);
      ExperimentalDesign::RunRows rows;
      for (const auto &f : ms_run_paths)
      {
        ExperimentalDesign::RunRow r;
        r.path = f;
        r.fraction = 1;
        r.sample = sample;
        r.run = sample;
        r.channel = 1; // TODO MULTIPLEXING: adapt for non-label-free
        rows.push_back(r);
        ++sample;
      }
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (ConsensusMap derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    ExperimentalDesign ExperimentalDesign::fromFeatureMap(const FeatureMap &fm)
    {
      ExperimentalDesign experimental_design;
      // path of the original MS run (mzML / raw file)
      StringList ms_paths;
      fm.getPrimaryMSRunPath(ms_paths);

      if (ms_paths.size() != 1)
      {
        throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "FeatureMap annotated with " + String(ms_paths.size()) + " MS files. Must be exactly one.");
      }

      // Feature map is simple. One file, one fraction, one sample, one run
      ExperimentalDesign::RunRow r;
      r.path = ms_paths[0];
      r.fraction = 1;
      r.sample = 1;
      r.run = 1;
      r.channel = 1;

      ExperimentalDesign::RunRows rows(1, r);
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (FeatureMap derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    ExperimentalDesign ExperimentalDesign::fromIdentifications(const vector <ProteinIdentification> &proteins)
    {
      ExperimentalDesign experimental_design;
      // path of the original MS files (mzML / raw file)
      StringList ms_run_paths;
      for (const auto &protein : proteins)
      {
        StringList tmp_ms_run_paths;
        protein.getPrimaryMSRunPath(tmp_ms_run_paths);
        if (tmp_ms_run_paths.size() != 1)
        {
          throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "ProteinIdentification annotated with " + String(tmp_ms_run_paths.size()) +
            " MS files. Must be exactly one.");
        }
        ms_run_paths.push_back(tmp_ms_run_paths[0]);
      }

      // no fractionation -> as many runs as samples
      // each identification run corresponds to one sample abundance
      unsigned sample(1);
      ExperimentalDesign::RunRows rows;
      for (const auto &f : ms_run_paths)
      {
        ExperimentalDesign::RunRow r;
        r.path = f;
        r.fraction = 1;
        r.sample = sample;
        r.run = sample;
        r.channel = 1;

        rows.push_back(r);
        ++sample;
      }
      experimental_design.setRunSection(rows);
      LOG_INFO << "Experimental design (Identification derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  channels: " << experimental_design.getNumberOfChannels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    map<unsigned, vector<String> > ExperimentalDesign::getFractionToMSFilesMapping() const
    {
      map<unsigned, vector<String> > ret;

      for (RunRow const &r : run_section_)
      {
        ret[r.fraction].emplace_back(r.path);
      }
      return ret;
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesign::pathChannelMapper(
            const bool basename,
            unsigned (*f)(const ExperimentalDesign::RunRow &row)) const
    {
      map<pair<String, unsigned>, unsigned> ret;
      for (RunRow const &r : run_section_)
      {
        const String path = String(r.path);
        pair<String, unsigned> tpl = make_pair((basename ? File::basename(path) : path), r.channel);
        ret[tpl] = f(r);
      }
      return ret;
    }


    map<pair<String, unsigned>, unsigned> ExperimentalDesign::getPathChannelToSampleMapping(
            const bool basename) const
    {
      return pathChannelMapper(basename, [](const RunRow &r)
      { return r.sample; });
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesign::getPathChannelToFractionMapping(
            const bool basename) const
    {
      return pathChannelMapper(basename, [](const RunRow &r)
      { return r.fraction; });
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesign::getPathChannelToRunMapping(
            const bool basename) const
    {
      return pathChannelMapper(basename, [](const RunRow &r)
      { return r.run; });
    }


    bool ExperimentalDesign::sameNrOfMSFilesPerFraction() const
    {
      map<unsigned, vector<String>> frac2files = getFractionToMSFilesMapping();
      if (frac2files.size() <= 1) { return true; }

      Size files_per_fraction(0);
      for (auto const &f : frac2files)
      {
        if (files_per_fraction == 0) // first fraction, initialize
        {
          files_per_fraction = f.second.size();
        }
        else // fraction >= 2
        {
          // different number of associated MS files?
          if (f.second.size() != files_per_fraction)
          {
            return false;
          }
        }
      }
      return true;
    }

    /* Implementation for the Experimental Design and RunSection */

    const RunRows& ExperimentalDesign::getRunSection() const
    {
      return run_section_;
    }


    void ExperimentalDesign::setRunSection(const RunRows& run_section)
    {
      run_section_ = run_section;
      sort_();
      checkValidRunSection_();
    }

    void ExperimentalDesign::setSampleSection(const SampleSection& sample_section)
    {
      sample_section_ = sample_section;
    }

    unsigned ExperimentalDesign::getNumberOfSamples() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
                              [](const RunRow& f1, const RunRow& f2)
                              {
                                return f1.sample < f2.sample;
                              })->sample;
    }

    unsigned ExperimentalDesign::getNumberOfFractions() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
                              [](const RunRow& f1, const RunRow& f2)
                              {
                                  return f1.fraction < f2.fraction;
                              })->fraction;
    }

    // @return the number of channels per file
    unsigned ExperimentalDesign::getNumberOfChannels() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
                              [](const RunRow& f1, const RunRow& f2)
                              {
                                return f1.fraction < f2.fraction;
                              })->channel;
    }

    // @return the number of MS files (= fractions * runs)
    unsigned ExperimentalDesign::getNumberOfMSFiles() const
    {
      std::set<std::string> unique_paths;
      for (auto const & r : run_section_) { unique_paths.insert(r.path); }
      return unique_paths.size();
    }

    bool ExperimentalDesign::isFractionated() const
    {
      std::vector<unsigned> fractions = this->getFractions();
      std::set<unsigned> fractions_set(fractions.begin(), fractions.end());
      return fractions_set.size() < 2;
    }

    // @return the number of runs (before fractionation)
    // Allows to group fraction ids and source files
    unsigned ExperimentalDesign::getNumberOfPrefractionationRuns() const
    {
      if (run_section_.empty()) { return 0; }
      return std::max_element(run_section_.begin(), run_section_.end(),
                              [](const RunRow& f1, const RunRow& f2)
                              {
                                  return f1.run < f2.run;
                              })->run;
    }

    // @return sample index (depends on run and channel)
    unsigned ExperimentalDesign::getSample(unsigned run, unsigned channel)
    {
      return std::find_if(run_section_.begin(), run_section_.end(),
                          [&run, &channel](const RunRow& r)
                          {
                              return r.run == run && r.channel == channel;
                          })->sample;
    }

    const ExperimentalDesign::SampleSection& ExperimentalDesign::getSampleSection() const
    {
      return sample_section_;
    }

    std::vector< String > ExperimentalDesign::getFileNames(const bool basename) const
    {
      std::vector<String> filenames;
      for (const RunRow &row : run_section_)
      {
        const String path = String(row.path);
        filenames.push_back(basename ? path : File::basename(path));
      }
      return filenames;
    }

    template<typename T>
    void ExperimentalDesign::errorIfAlreadyExists(std::set<T> &container, T &item, const String &message)
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

    void ExperimentalDesign::checkValidRunSection_()
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

  std::vector<unsigned> ExperimentalDesign::getChannels() const
  {
    std::vector<unsigned> channels;
    for (const RunRow &row : run_section_)
    {
      channels.push_back(row.channel);
    }
    return channels;
  }

  std::vector<unsigned> ExperimentalDesign::getFractions() const
  {
    std::vector<unsigned> fractions;
    for (const RunRow &row : run_section_)
    {
      fractions.push_back(row.fraction);
    }
    return fractions;
  }

  void ExperimentalDesign::sort_()
  {
    std::sort(run_section_.begin(), run_section_.end(),
                [](const RunRow& a, const RunRow& b)
                {
                  return std::tie(a.run, a.fraction, a.channel, a.sample, a.path) <
                         std::tie(b.run, b.fraction, b.channel, b.sample, b.path);
                });
  }

  /* Implementations of SampleSection */

  std::set<unsigned> ExperimentalDesign::SampleSection::getSamples() const
  {
    std::set<unsigned> samples;
    for (const auto &kv : sample_to_rowindex_)
    {
      samples.insert(kv.first);
    }
    return samples;
  }

  std::set< String > ExperimentalDesign::SampleSection::getFactors() const
  {
    std::set<String> factors;
    for (const auto &kv : columnname_to_columnindex_)
    {
      factors.insert(kv.first);
    }
    return factors;
  }

  bool ExperimentalDesign::SampleSection::hasSample(const unsigned sample) const
  {
    return sample_to_rowindex_.find(sample) != sample_to_rowindex_.end();
  }

  bool ExperimentalDesign::SampleSection::hasFactor(const String &factor) const
  {
    return columnname_to_columnindex_.find(factor) != columnname_to_columnindex_.end();
  }

  String ExperimentalDesign::SampleSection::getFactorValue(const unsigned sample, const String &factor)
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
}
