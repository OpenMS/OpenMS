// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ExperimentalDesign.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <QtCore/QString>
#include <QtCore/QFileInfo>

#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{
    using MSFileSection = std::vector<ExperimentalDesign::MSFileSection>;

    ExperimentalDesign::SampleSection::SampleSection(
        const std::vector< std::vector < String > >& content,
        const std::map< String, Size >& sample_to_rowindex,
        const std::map< String, Size >& columnname_to_columnindex
      ) : 
      content_(content),
      sample_to_rowindex_(sample_to_rowindex),
      columnname_to_columnindex_(columnname_to_columnindex) 
    {    
    }

    ExperimentalDesign::ExperimentalDesign(
      const ExperimentalDesign::MSFileSection& msfile_section, 
      const ExperimentalDesign::SampleSection& sample_section) : 
        msfile_section_(msfile_section), 
        sample_section_(sample_section)
    {
      sort_();
      isValid_();
    }

    ExperimentalDesign ExperimentalDesign::fromConsensusMap(const ConsensusMap &cm)
    {
      ExperimentalDesign experimental_design;

      // one of label-free, labeled_MS1, labeled_MS2
      const String & experiment_type = cm.getExperimentType();

      // path of the original MS run (mzML / raw file)
      StringList ms_run_paths;
      cm.getPrimaryMSRunPath(ms_run_paths);
      
      // Note: consensus elements of the same fraction group corresponds to one sample abundance
      ExperimentalDesign::MSFileSection msfile_section;
      ExperimentalDesign::SampleSection sample_section;

      // determine vector of ms file names (in order of appearance)
      vector<String> msfiles;
      std::map<pair<UInt,UInt>, UInt> fractiongroup_label_to_sample_mapping;
      std::map<String, UInt> samplename_to_sample_mapping;

      for (const auto &f : cm.getColumnHeaders())
      {
        if (std::find(msfiles.begin(), msfiles.end(), f.second.filename) == msfiles.end())
        {
          msfiles.push_back(f.second.filename);
        }
      }

      bool no_fractions = true;
      for (const auto &f : cm.getColumnHeaders())
      {
        ExperimentalDesign::MSFileSectionEntry r;
        r.path = f.second.filename;
        if (f.second.metaValueExists("fraction"))
        {
          no_fractions = false;
          r.fraction = static_cast<unsigned int>(f.second.getMetaValue("fraction"));

          if (f.second.metaValueExists("fraction_group"))
          {
            r.fraction_group = static_cast<unsigned int>(f.second.getMetaValue("fraction_group"));
          }
          else
          {
            // if we have annotated fractions, we need to also know how these are grouped
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Fractions annotated but no grouping provided.");
          }
        }
        else
        { // no fractions and fraction group information annotated, deduce from data
          r.fraction = 1;

          // no fractions -> one fraction group for each MS file
          // to create a unique group identifier [1..n], we take the index of the (unique) filenames  
          r.fraction_group = (std::find(msfiles.begin(), msfiles.end(), f.second.filename) 
                            - msfiles.begin()) + 1;
        }
        if (no_fractions) OPENMS_LOG_INFO << "No fractions annotated in consensusXML. Assuming unfractionated." << endl;

        r.label = f.second.getLabelAsUInt(experiment_type);

        if (!f.second.metaValueExists("sample_name"))
        {
          // We create one sample for each fractiongroup/label combination
          // this assumes that fractionation took place after labelling. Otherwise, a design needs to be given.
          // check fractiongroup_label_to_sample_mapping and add if not present, otherwise use present
          auto key = make_pair(r.fraction_group, r.label);
          auto it = fractiongroup_label_to_sample_mapping.emplace(key, fractiongroup_label_to_sample_mapping.size());
          r.sample = it.first->second;
          r.sample_name = r.sample;
        } else {
          r.sample_name = f.second.getMetaValue("sample_name");
          [[maybe_unused]] const auto& [it, inserted] = samplename_to_sample_mapping.emplace(r.sample_name,samplename_to_sample_mapping.size());
          r.sample = it->second;
        }

        msfile_section.push_back(r);
        if (!sample_section.hasSample(r.sample))
          sample_section.addSample(r.sample);

      }

      experimental_design.setMSFileSection(msfile_section);
      experimental_design.setSampleSection(sample_section);
      OPENMS_LOG_DEBUG << "Experimental design (ConsensusMap derived):\n"
               << "  Files: " << experimental_design.getNumberOfMSFiles()
               << "  Fractions: " << experimental_design.getNumberOfFractions()
               << "  Labels: " << experimental_design.getNumberOfLabels()
               << "  Samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    void ExperimentalDesign::SampleSection::addSample(const String& samplename, const vector<String>& content)
    {
      //TODO warn when already present? Overwrite?
      //TODO check content size
      sample_to_rowindex_.try_emplace(samplename, sample_to_rowindex_.size());
      content_.push_back(content);
    }

    String ExperimentalDesign::SampleSection::getSampleName(unsigned sample_row) const
    {
      return content_.at(sample_row).at(columnname_to_columnindex_.at("Sample"));
    }

    unsigned ExperimentalDesign::SampleSection::getSampleRow(const String& sample) const
    {
      return sample_to_rowindex_.at(sample);
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

      // Feature map is simple. One file, one fraction, one sample, one fraction_group
      ExperimentalDesign::MSFileSectionEntry r;
      r.path = ms_paths[0];
      r.fraction = 1;
      r.sample = 0;
      r.sample_name = "0";
      r.fraction_group = 1;
      r.label = 1;

      ExperimentalDesign::MSFileSection rows(1, r);

      ExperimentalDesign::SampleSection s;
      s.addSample(r.sample_name);

      experimental_design.setMSFileSection(rows);
      experimental_design.setSampleSection(s);

      OPENMS_LOG_INFO << "Experimental design (FeatureMap derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  labels: " << experimental_design.getNumberOfLabels()
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
        // Due to merging it is valid to have multiple origins in a ProteinIDRun
        StringList tmp_ms_run_paths;
        protein.getPrimaryMSRunPath(tmp_ms_run_paths);
        ms_run_paths.insert(ms_run_paths.end(), tmp_ms_run_paths.begin(), tmp_ms_run_paths.end());
        //TODO think about uniquifying or warning if duplicates occur
      }

      // For now and without further info we have to assume a labelfree, unfractionated experiment
      // no fractionation -> as many fraction_groups as samples
      // each identification run corresponds to one sample abundance
      unsigned sample(0);
      ExperimentalDesign::MSFileSection rows;
      ExperimentalDesign::SampleSection srows;
      for (const auto &f : ms_run_paths)
      {
        ExperimentalDesign::MSFileSectionEntry r;
        r.path = f;
        r.fraction = 1;
        r.sample = sample;
        r.sample_name = String(sample);
        r.fraction_group = sample;
        r.label = 1;

        rows.push_back(r);
        srows.addSample(String(sample));
        ++sample;
      }
      experimental_design.setMSFileSection(rows);
      experimental_design.setSampleSection(srows);
      OPENMS_LOG_INFO << "Experimental design (Identification derived):\n"
               << "  files: " << experimental_design.getNumberOfMSFiles()
               << "  fractions: " << experimental_design.getNumberOfFractions()
               << "  labels: " << experimental_design.getNumberOfLabels()
               << "  samples: " << experimental_design.getNumberOfSamples() << "\n"
               << endl;
      return experimental_design;
    }

    map<unsigned, vector<String> > ExperimentalDesign::getFractionToMSFilesMapping() const
    {
      map<unsigned, vector<String> > ret;

      for (MSFileSectionEntry const& r : msfile_section_)
      {
        ret[r.fraction].emplace_back(r.path);
      }
      return ret;
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesign::pathLabelMapper_(
            const bool basename,
            unsigned (*f)(const ExperimentalDesign::MSFileSectionEntry &entry)) const
    {
      map<pair<String, unsigned>, unsigned> ret;
      for (MSFileSectionEntry const& r : msfile_section_)
      {
        const String path = String(r.path);
        pair<String, unsigned> tpl = make_pair((basename ? File::basename(path) : path), r.label);
        ret[tpl] = f(r);
      }
      return ret;
    }

    map<vector<String>, set<String>> ExperimentalDesign::getUniqueSampleRowToSampleMapping() const
    {
      map<vector<String>, set<String> > rowContent2RowIdx;
      auto factors = sample_section_.getFactors();
      assert(!factors.empty());

      factors.erase("Sample"); // we do not care about ID in duplicates

      for (const String& u : sample_section_.getSamples())
      {
        std::vector<String> valuesToHash{};
        valuesToHash.reserve(factors.size());
        for (const String& fac : factors)
        {
          valuesToHash.emplace_back(sample_section_.getFactorValue(u, fac));
        }
        auto emplace_pair = rowContent2RowIdx.emplace(valuesToHash, set<String>{});
        emplace_pair.first->second.insert(u);
      }

      return rowContent2RowIdx;
    }

    map<String, unsigned> ExperimentalDesign::getSampleToPrefractionationMapping() const
    {
      map<String, unsigned> res;

      // could happen when the Experimental Design was loaded from an idXML or consensusXML
      // without additional Experimental Design file
      if (sample_section_.getFactors().empty())
      {
        // no information about the origin of the samples -> assume uniqueness of all
        Size i(0);
        for (const auto& s : sample_section_.getSamples())
        {
          res[s] = i;
          i++;
        }
      }
      else
      {
        const map<vector<String>, set<String>>& rowContent2RowIdx = getUniqueSampleRowToSampleMapping();
        Size s(0);
        for (const auto &condition : rowContent2RowIdx)
        {
          for (auto &sample : condition.second)
          {
            res.emplace(sample, s);
          }
          ++s;
        }
      }
      return res;
    }

    map<vector<String>, set<unsigned>> ExperimentalDesign::getConditionToSampleMapping() const
    {
      const auto& facset = sample_section_.getFactors();
      // assert(!facset.empty()); // not needed: If no factors are given, same condition is assumed for every run
      set<String> nonRepFacs{};

      for (const String& fac : facset)
      {
        if (fac != "Sample" && !fac.hasSubstring("replicate") && !fac.hasSubstring("Replicate"))
        {
          nonRepFacs.insert(fac);
        }
      }

      map<vector<String>, set<unsigned> > rowContent2RowIdx;
      for (const auto& u : sample_section_.getSamples())
      {
        std::vector<String> valuesToHash{};
        valuesToHash.reserve(nonRepFacs.size());
        for (const String& fac : nonRepFacs)
        {
          valuesToHash.emplace_back(sample_section_.getFactorValue(u, fac));
        }
        auto emplace_pair = rowContent2RowIdx.emplace(valuesToHash, set<unsigned>{});
        emplace_pair.first->second.insert(sample_section_.getSampleRow(u));
      }
      return rowContent2RowIdx;
    }

    map<String, unsigned> ExperimentalDesign::getSampleToConditionMapping() const
    {
      map<String, unsigned> res;
      // could happen when the Experimental Design was loaded from an idXML or consensusXML
      // without additional Experimental Design file
      if (sample_section_.getFactors().empty())
      {
        // no information about the origin of the samples -> assume uniqueness of all
        unsigned nr(getNumberOfSamples());
        for (unsigned i(0); i <= nr; ++i)
        {
          res[i] = i;
        }
      }
      else
      {
        const map<vector<String>, set<unsigned>>& rowContent2RowIdx = getConditionToSampleMapping();
        Size s(0);
        for (const auto &condition : rowContent2RowIdx)
        {
          for (auto &sample : condition.second)
          {
            res.emplace(sample, s);
          }
          ++s;
        }
      }
      return res;
    }

    vector<vector<pair<String, unsigned>>> ExperimentalDesign::getConditionToPathLabelVector() const
    {
      const map<vector<String>, set<unsigned>>& rowContent2RowIdx = getConditionToSampleMapping();

      const map<pair<String, unsigned>, unsigned>& pathLab2Sample = getPathLabelToSampleMapping(false);
      vector<vector<pair<String, unsigned>>> res{rowContent2RowIdx.size()};
      Size s(0);
      // ["wt","24h","10mg"] -> sample [1, 3]
      for (const auto& rcri : rowContent2RowIdx)
      {
        // sample 1
        for (const auto& ri : rcri.second)
        {
          // [foo.mzml, ch1] -> sample 1
          for (const auto& pl2Sample : pathLab2Sample)
          {
            // sample 1 == sample 1
            if (pl2Sample.second == ri)
            {
              // res[0] -> [[foo, 1],...]
              res[s].emplace_back(pl2Sample.first);
            }
          }
        }
        ++s;
      }
      return res;
    }

    map<pair< String, unsigned >, unsigned> ExperimentalDesign::getPathLabelToPrefractionationMapping(const bool basename) const
    {
      const auto& sToPreFrac = getSampleToPrefractionationMapping();
      const auto& pToS = getPathLabelToSampleMapping(basename);
      map<pair<String, unsigned>, unsigned> ret;
      for (const auto& entry : pToS)
      {
        ret.emplace(entry.first, sToPreFrac.at(entry.second));
      }
      return ret;
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesign::getPathLabelToConditionMapping(const bool basename) const
    {
      const auto& sToC = getSampleToConditionMapping();
      const auto& pToS = getPathLabelToSampleMapping(basename);
      map<pair<String, unsigned>, unsigned> ret;
      for (const auto& entry : pToS)
      {
        ret.emplace(entry.first, sToC.at(entry.second));
      }
      return ret;
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesign::getPathLabelToSampleMapping(
            const bool basename) const
    {
      return pathLabelMapper_(basename, [](const MSFileSectionEntry &r)
      { return r.sample; });
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesign::getPathLabelToFractionMapping(
            const bool basename) const
    {
      return pathLabelMapper_(basename, [](const MSFileSectionEntry &r)
      { return r.fraction; });
    }

    map<pair<String, unsigned>, unsigned> ExperimentalDesign::getPathLabelToFractionGroupMapping(
            const bool basename) const
    {
      return pathLabelMapper_(basename, [](const MSFileSectionEntry &r)
      { return r.fraction_group; });
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

    const ExperimentalDesign::MSFileSection& ExperimentalDesign::getMSFileSection() const
    {
      return msfile_section_;
    }

    void ExperimentalDesign::setMSFileSection(const MSFileSection& msfile_section)
    {
      msfile_section_ = msfile_section;
      sort_();
    }

    void ExperimentalDesign::setSampleSection(const SampleSection& sample_section)
    {
      sample_section_ = sample_section;
    }

    unsigned ExperimentalDesign::getNumberOfSamples() const
    {
      /*
      if (msfile_section_.empty()) { return 0; }
      return std::max_element(msfile_section_.begin(), msfile_section_.end(),
                              [](const MSFileSectionEntry& f1, const MSFileSectionEntry& f2)
                              {
                                return f1.sample < f2.sample;
                              })->sample;
      */
      return sample_section_.getContentSize();
    }

    // TODO IMHO this is ill-defined and only works if exactly the same identifiers are used
    //   across fraction groups. Also, the number of fractions can change across
    //   fraction groups, e.g., if one is missing.
    unsigned ExperimentalDesign::getNumberOfFractions() const
    {
      /*
      if (msfile_section_.empty()) { return 0; }
      return std::max_element(msfile_section_.begin(), msfile_section_.end(),
                              [](const MSFileSectionEntry& f1, const MSFileSectionEntry& f2)
                              {
                                  return f1.fraction < f2.fraction;
                              })->fraction;
      */
      auto fs = set<Size>();
      for (const auto& row: msfile_section_)
      {
        fs.insert(row.fraction);
      }
      return fs.size();
    }

    // @return the number of labels per file
    unsigned ExperimentalDesign::getNumberOfLabels() const
    {
      if (msfile_section_.empty()) { return 0; }
      return std::max_element(msfile_section_.begin(), msfile_section_.end(),
                              [](const MSFileSectionEntry& f1, const MSFileSectionEntry& f2)
                              {
                                return f1.label < f2.label;
                              })->label;
    }

    // @return the number of MS files (= fractions * fraction_groups)
    unsigned ExperimentalDesign::getNumberOfMSFiles() const
    {
      std::set<std::string> unique_paths;
      for (auto const & r : msfile_section_) 
      { 
        unique_paths.insert(r.path);
      }
      return unique_paths.size();
    }

    bool ExperimentalDesign::isFractionated() const
    {
      // TODO Warning: only works if fractions are always the same in every fraction group!
      std::vector<unsigned> fractions = getFractions_();
      std::set<unsigned> fractions_set(fractions.begin(), fractions.end());
      return fractions_set.size() > 1;
    }

    // TODO IMHO this is ill-defined and only works if fraction group names/IDs change over every sample
    unsigned ExperimentalDesign::getNumberOfFractionGroups() const
    {
      /*
      if (msfile_section_.empty()) { return 0; }
      return std::max_element(msfile_section_.begin(), msfile_section_.end(),
                              [](const MSFileSectionEntry& f1, const MSFileSectionEntry& f2)
                              {
                                  return f1.fraction_group < f2.fraction_group;
                              })->fraction_group;
      */
      auto fgs = set<Size>();
      for (const auto& row: msfile_section_)
      {
        fgs.insert(row.fraction_group);
      }
      return fgs.size();
    }

    unsigned ExperimentalDesign::getSample(unsigned fraction_group, unsigned label)
    {
      return std::find_if(msfile_section_.begin(), msfile_section_.end(),
                          [&fraction_group, &label](const MSFileSectionEntry& r)
                          {
                              return r.fraction_group == fraction_group && r.label == label;
                          })->sample;
    }

    const ExperimentalDesign::SampleSection& ExperimentalDesign::getSampleSection() const
    {
      return sample_section_;
    }

    std::vector< String > ExperimentalDesign::getFileNames_(const bool basename) const
    {
      std::vector<String> filenames;
      for (const MSFileSectionEntry& row : msfile_section_)
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

    void ExperimentalDesign::isValid_()
    {
      std::set< std::tuple< unsigned, unsigned, unsigned > > fractiongroup_fraction_label_set;
      std::set< std::tuple< std::string, unsigned > > path_label_set;
      std::map< std::tuple< unsigned, unsigned >, std::set< unsigned > > fractiongroup_label_to_sample;
      std::set< unsigned > label_set;
      std::set< unsigned > fraction_group_set;
      std::set< unsigned > sample_set;

      if (msfile_section_.empty())
      {
        return;
      }

      bool labelfree = true;
      for (const MSFileSectionEntry& row : msfile_section_)
      {

        // FRACTIONGROUP_FRACTION_LABEL TUPLE
        std::tuple<unsigned, unsigned, unsigned> fractiongroup_fraction_label = std::make_tuple(row.fraction_group, row.fraction, row.label);
        errorIfAlreadyExists(
          fractiongroup_fraction_label_set,
          fractiongroup_fraction_label,
        "(Fraction Group, Fraction, Label) combination can only appear once");

        // PATH_LABEL_TUPLE
        std::tuple<std::string, unsigned> path_label = std::make_tuple(row.path, row.label);
        errorIfAlreadyExists(
          path_label_set,
          path_label,
          "(Path, Label) combination (" + String(get<0>(path_label)) + "," + String(get<1>(path_label)) + ") can only appear once");

        // FRACTIONGROUP_LABEL TUPLE
        std::tuple<unsigned, unsigned> fractiongroup_label = std::make_tuple(row.fraction_group, row.label);
        fractiongroup_label_to_sample[fractiongroup_label].insert(row.sample);

        label_set.insert(row.label);
        fraction_group_set.insert(row.fraction_group);
      }
      if (label_set.size() > 1)
      {
        labelfree = false;
      }

      // TODO remove that restriction. Check if that order is still assumed somewhere.
      //  Except for in the getNrFractionGroups
      if ((*fraction_group_set.begin()) != 1)
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Fraction groups have to be integers and their set needs to be consecutive and start with 1.",
          String(*fraction_group_set.begin()));
      }
      Size s = 0;
      for (const auto& fg : fraction_group_set)
      {
        ++s;
        if (fg != s)
        {
          throw Exception::InvalidValue(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Fraction groups have to be integers and their set needs to be consecutive and start with 1.",
            String(*fraction_group_set.begin()));
        }
      }

      for (const auto& [k,v] : fractiongroup_label_to_sample)
      {
        if (labelfree && v.size() > 1)
        {
          throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION, "Multiple samples encountered for the same fraction group and the same label"
                                    "Please correct your experimental design if this is a label free experiment, otherwise"
                                    "check your labels. Occurred at fraction group " +
                                    String(std::get<0>(k)) + "and label" + String(std::get<1>(k)));
        }
      }
    }

    std::vector<unsigned> ExperimentalDesign::getLabels_() const
    {
      std::vector<unsigned> labels;
      for (const MSFileSectionEntry &row : msfile_section_)
      {
        labels.push_back(row.label);
      }
      return labels;
    }

    std::vector<unsigned> ExperimentalDesign::getFractions_() const
    {
      std::vector<unsigned> fractions;
      for (const MSFileSectionEntry &row : msfile_section_)
      {
        fractions.push_back(row.fraction);
      }
      return fractions;
    }

    void ExperimentalDesign::sort_()
    {
      std::sort(msfile_section_.begin(), msfile_section_.end(),
        [](const MSFileSectionEntry& a, const MSFileSectionEntry& b)
        {
          return std::tie(a.fraction_group, a.fraction, a.label, a.sample, a.path) <
                 std::tie(b.fraction_group, b.fraction, b.label, b.sample, b.path);
        });
    }

    Size ExperimentalDesign::filterByBasenames(const set<String>& bns)
    {
      Size before = msfile_section_.size();
      msfile_section_.erase(std::remove_if(msfile_section_.begin(), msfile_section_.end(),
      [&bns](MSFileSectionEntry& e)
      {
        return bns.find(File::basename(e.path)) == bns.end();
      }), msfile_section_.end());

      int diff = before - msfile_section_.size();
      if (diff > 0)
      {
        OPENMS_LOG_WARN << "Removed " << diff << " files from design to match given mzML/idXML subset." << std::endl;
      }

      if (msfile_section_.empty())
      {
        OPENMS_LOG_FATAL_ERROR << "Given basename set does not overlap with design. Design would be empty." << std::endl;
      }

      return (before - msfile_section_.size());
    }

    /* Implementations of SampleSection */

    std::set<String> ExperimentalDesign::SampleSection::getSamples() const
    {
      std::set<String> samples;
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

    bool ExperimentalDesign::SampleSection::hasSample(const String& sample) const
    {
      return sample_to_rowindex_.find(sample) != sample_to_rowindex_.end();
    }

    bool ExperimentalDesign::SampleSection::hasFactor(const String &factor) const
    {
      return columnname_to_columnindex_.find(factor) != columnname_to_columnindex_.end();
    }

    String ExperimentalDesign::SampleSection::getFactorValue(unsigned sample_idx, const String &factor) const
    {

      if (!hasFactor(factor))
      {
        throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Factor " + factor + " is not present in the Experimental Design");
      }
      const StringList& sample_row = content_.at(sample_idx);
      const Size col_index = columnname_to_columnindex_.at(factor);
      return sample_row[col_index];
    }

    String ExperimentalDesign::SampleSection::getFactorValue(const String& sample_name, const String &factor) const
    {
     if (!hasSample(sample_name))
     {
      throw Exception::MissingInformation(
                  __FILE__,
                  __LINE__,
                  OPENMS_PRETTY_FUNCTION,
                  "Sample " + sample_name + " is not present in the Experimental Design");
     }
     if (!hasFactor(factor))
     {
      throw Exception::MissingInformation(
                  __FILE__,
                  __LINE__,
                  OPENMS_PRETTY_FUNCTION,
                  "Factor " + factor + " is not present in the Experimental Design");
     }
     const StringList& sample_row = content_.at(sample_to_rowindex_.at(sample_name));
     const Size col_index = columnname_to_columnindex_.at(factor);
     return sample_row[col_index];
    }

    Size ExperimentalDesign::SampleSection::getFactorColIdx(const String &factor) const
    {
      if (! hasFactor(factor))
      {
        throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Factor " + factor + " is not present in the Experimental Design");
      }
      return columnname_to_columnindex_.at(factor);
    }

    Size ExperimentalDesign::SampleSection::getContentSize() const
    {
      return content_.size();
    }
}
