// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_TopDownConsensusFeatureGroup TopDownConsensusFeatureGroup
  @brief TopDownConsensusFeatureGroup build ConsensusFeatureGroup from FLASHQuant outputs
**/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TopDownConsensusFeatureGroup:
    public TOPPBase,
    public ProgressLogger
{
public:
  TopDownConsensusFeatureGroup():
      TOPPBase("TopDownConsensusFeatureGroup", "TopDownConsensusFeatureGroup from FLASHQuant", false, {}, false),
      ProgressLogger()
  {
  }

private:
  String QUANT_METHOD_;
  double MASS_TOL_;
  String MASS_TOL_UNIT_;
  int RT_TOL_;
  Size REP_COUNT_;

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files to align (*.tsv)", true);
    setValidFormats_("in", ListUtils::create<String>("tsv"));
    registerOutputFile_("out", "<file>", "", "Output tsv file", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerDoubleOption_("mass_tol", "<value>", 100, "Mass tolerance.", false);
    setMinFloat_("mass_tol", 0.0);
    registerStringOption_("mass_tol_unit", "<choice>", "ppm", "Mass tolerance unit", false);
    setValidStrings_("mass_tol_unit", {"ppm", "dalton"});
    registerIntOption_("rt_tol", "<integer>", 180, "Retention time tolerance for MedianApexRetentionTime in second", false);
    setMinInt_("rt_tol", 0);
    registerStringOption_("quant_method", "<choice>", "FeatureGroupQuantity", "Quantity value to use from FLASHQuant result", false);
    setValidStrings_("quant_method", {"FeatureGroupQuantity", "AllAreaUnderTheCurve", "SumIntensity"});
    registerStringOption_("consensus_as_input", "<choice>", "false", "Set it true when input files are consensus files", false);
    setValidStrings_("consensus_as_input", {"false", "true"});
  }

  struct FeatureGroup
  {
    Size rep_index; /// replicate index
    Size fgroup_index; /// FeatureGroupIndex in FLASHQuant
    double mass; /// MonoisotopicMass from FLASHQuant
    double apex_rt; /// MedianApexRetentionTime from FLASHQuant
    double abundance; /// Abundance from quant_option

    bool operator<(const FeatureGroup& fg) const
    {
      return this->mass < fg.mass;
    }
  };

  struct ConsensusFeatureGroup
  {
    double avg_mass; /// median of MonoisotopicMass from FLASHQuant
    double avg_apex_rt; /// MedianApexRetentionTime from FLASHQuant
    double cv;
    std::vector<Size> fgroup_indices; /// FeatureGroupIndex in FLASHQuant. Sorted by rep_index
    std::vector<double> abundances; /// values of QUANT_METHOD in FLASHQuant. Sorted by rep_index

    void calculate_cv()
    {
      double mean = Math::mean(abundances.begin(), abundances.end());
      double std = Math::sd(abundances.begin(), abundances.end(), mean);
      cv = std/mean;
    }
  };

  void writeConsensusFeatureGroupInTsv(std::vector<ConsensusFeatureGroup> &consensus, String& out_path, Size input_file_count)
  {
    String header = "ConsensusFeatureGroupIndex\tAvgMonoisotopicMass\tAvgApexRetentionTime\tCoefficientOfVariation\t";
    String abundance_col_tag = "AbundanceInFile";
    String fg_index_col_tag = "FeatureGroupIndexInFile";
    String header_abundance = "";
    String header_fg_index = "";
    for (Size i = 0; i < input_file_count; ++i)
    {
      header_abundance += abundance_col_tag + i + "\t";
      header_fg_index += fg_index_col_tag + i + "\t";
    }
    header_fg_index.pop_back();
    header += header_abundance + header_fg_index + "\n";

    // start writing
    std::ofstream os(out_path);
    os << header;
    for (Size i = 0; i < consensus.size(); i++)
    {
      auto &tmp = consensus[i];
      String line = to_string(i) + "\t" + to_string(tmp.avg_mass) + "\t" + to_string(tmp.avg_apex_rt) + "\t" + to_string(tmp.cv) + "\t";

      for (auto &a : tmp.abundances)
      {
        line += to_string(a) + "\t";
      }
      for (auto &f : tmp.fgroup_indices)
      {
        line += to_string(f) + "\t";
      }
      line.pop_back();
      os << line << std::endl;
    }
    os.close();
  }

  void readFLASHQuantResultFile(String &filepath, std::vector<FeatureGroup> &out_fgroups, Size rep_index)
  {
    std::vector<FeatureGroup> fgroups;

    std::ifstream data(filepath);
    std::string line;
    std::string tmp;

    // read header
    std::map<std::string, Size> header_dict;
    TextFile::getLine(data, line);
    std::stringstream tmp_lstream(line);
    Size i = 0;
    while (std::getline(tmp_lstream, tmp, '\t'))
    {
      header_dict[tmp] = i++;
    }

    // read data
    while(TextFile::getLine(data, line)) // iterate over lines
    {
      std::stringstream lstream(line);
      std::vector<String> tmp_line;
      while (std::getline(lstream, tmp, '\t')) // iterate over column
      {
        tmp_line.push_back(tmp);
      }
      FeatureGroup fg;
      fg.rep_index = rep_index;
      fg.fgroup_index = (Size) tmp_line[header_dict.at("FeatureGroupIndex")].toInt();
      fg.mass = tmp_line[header_dict.at("MonoisotopicMass")].toDouble();
      fg.apex_rt = tmp_line[header_dict.at("MedianApexRetentionTime")].toDouble();
      fg.abundance = tmp_line[header_dict.at(QUANT_METHOD_)].toDouble();
      fgroups.push_back(fg);
    }
    OPENMS_LOG_INFO << ", #FeatureGroup " << fgroups.size() << std::endl;

    // update out_fgroups with result
    out_fgroups.reserve(out_fgroups.size() + distance(fgroups.begin(), fgroups.end()));
    out_fgroups.insert(out_fgroups.end(),fgroups.begin(),fgroups.end());
  }

  void readFLASHQuantConsensusFile(String &filepath, std::vector<FeatureGroup> &out_fgroups, Size rep_index)
  {
    std::vector<FeatureGroup> fgroups;

    std::ifstream data(filepath);
    std::string line;
    std::string tmp;

    // read header
    std::map<std::string, Size> header_dict;
    TextFile::getLine(data, line);
    std::stringstream tmp_lstream(line);
    Size i = 0;
    while (std::getline(tmp_lstream, tmp, '\t'))
    {
      header_dict[tmp] = i++;
    }
    Size fg_index = header_dict.at("ConsensusFeatureGroupIndex");
    Size mass_index = header_dict.at("AvgMonoisotopicMass");
    Size rt_index = header_dict.at("AvgApexRetentionTime");
    std::set<Size> abundance_index_set;
    for (auto const& [key, index] : header_dict)
    {
      if (key.find("Abundance") != std::string::npos) {
        abundance_index_set.insert(index);
      }
    }

    // read data
    while(TextFile::getLine(data, line)) // iterate over lines
    {
      std::stringstream lstream(line);
      std::vector<String> tmp_line;
      while (std::getline(lstream, tmp, '\t')) // iterate over column
      {
        tmp_line.push_back(tmp);
      }
      FeatureGroup fg;
      fg.rep_index = rep_index;
      fg.fgroup_index = (Size) tmp_line[fg_index].toInt();
      fg.mass = tmp_line[mass_index].toDouble();
      fg.apex_rt = tmp_line[rt_index].toDouble();
      std::vector<double> masses;
      for (auto const &index : abundance_index_set)
      {
        masses.push_back(tmp_line[index].toDouble());
      }
      // average masses
      fg.abundance = std::accumulate(masses.begin(), masses.end(), 0.0) / masses.size();
      fgroups.push_back(fg);
    }
    OPENMS_LOG_INFO << ", #FeatureGroup " << fgroups.size() << std::endl;

    // update out_fgroups with result
    out_fgroups.reserve(out_fgroups.size() + distance(fgroups.begin(), fgroups.end()));
    out_fgroups.insert(out_fgroups.end(),fgroups.begin(),fgroups.end());
  }

  void computeConsensusFeatureGroup(std::vector<FeatureGroup> &fgroups, std::vector<ConsensusFeatureGroup> &consensus)
  {
    // sort by masses
    std::sort(fgroups.begin(), fgroups.end());

    consensus.clear();
    consensus.reserve(fgroups.size()/REP_COUNT_);

    // vector of pointer for fgroups. To link the leftover fgroups for computing consensus.
    vector<FeatureGroup*> fgroup_pointers;
    fgroup_pointers.reserve(fgroups.size());
    for (auto &f: fgroups)
    {
      fgroup_pointers.push_back(&f);
    }

    while(!fgroup_pointers.empty())
    {
      /// find the FeatureGroup with the maximum abundance
      auto iter = std::max_element(fgroup_pointers.begin(), fgroup_pointers.end(),
                                   [](const FeatureGroup* a, const FeatureGroup* b){return a->abundance < b->abundance;});
      Size reference_index = iter - fgroup_pointers.begin();
      FeatureGroup* reference_fg = *iter;
      double reference_mass = reference_fg->mass;
      double mass_tolerance = MASS_TOL_;
      if (MASS_TOL_UNIT_ == "ppm")
      {
        mass_tolerance *= reference_mass;
      }

      /// collect FeatureGroup within mass and rt tolerance
      std::vector<FeatureGroup*> candidate_fgs;
      std::vector<Size> candidate_indices; // index from fgroup_pointers
      // iter to right side (larger masses than reference FeatureGroup)
      ++iter; // move forward to reference+1 point
      for(; (iter<fgroup_pointers.end()) && ((*iter)->mass - reference_mass <= mass_tolerance) ; ++iter)
      {
        // check if RT tolerance matches
        if (std::abs((*iter)->apex_rt - reference_fg->apex_rt) > RT_TOL_)
        {
          continue;
        }
        // ignore the FeatureGroup from the same RepIndex to the reference
        if ((*iter)->rep_index == reference_fg->rep_index)
        {
          continue;
        }
        candidate_fgs.push_back(*iter);
        candidate_indices.push_back(iter - fgroup_pointers.begin());
      }
      // iter to left side (smaller masses than reference FeatureGroup)
      iter = fgroup_pointers.begin() + reference_index - 1; // back to reference-1 point
      for(; (iter>=fgroup_pointers.begin()) && (reference_mass - (*iter)->mass <= mass_tolerance) ; --iter)
      {
        if (std::abs((*iter)->apex_rt - reference_fg->apex_rt) > RT_TOL_)
        {
          continue;
        }
        // ignore the FeatureGroup from the same RepIndex to the reference
        if ((*iter)->rep_index == reference_fg->rep_index)
        {
          continue;
        }
        candidate_fgs.push_back(*iter);
        candidate_indices.push_back(iter - fgroup_pointers.begin());
      }

      /// check if the collected masses are from multiple replicates (Except for the reference)
      std::set<Size> rep_set;
      for (auto& tmp : candidate_fgs)
      {
        rep_set.insert(tmp->rep_index);
      }
      // remove reference featgroup, if not eligible for consensus
      if (rep_set.size() < REP_COUNT_-1)
      {
        fgroup_pointers.erase(fgroup_pointers.begin() + reference_index);
        continue;
      }

      /// Among the candidates from the same replicate, the one with the highest abundance wins.
      std::vector<FeatureGroup*> collected_fgs(REP_COUNT_);
      std::vector<Size> collected_indices; // index from fgroups (for erase, later)
      for (Size rep_index : rep_set)
      {
        double max_abundance = .0;
        Size chosen_index = 0; // index of candidate_fgs
        for (Size i = 0; i < candidate_fgs.size(); ++i)
        {
          if (candidate_fgs[i]->rep_index != rep_index) // ignore the candidate from different replicate
          {
            continue;
          }
          if (max_abundance > candidate_fgs[i]->abundance)
          {
            continue;
          }
          max_abundance = candidate_fgs[i]->abundance;
          chosen_index = i;
        }
        collected_fgs[rep_index] = candidate_fgs[chosen_index];
        collected_indices.push_back(candidate_indices[chosen_index]);
      }
      // add reference to the consensus fg list
      collected_fgs[reference_fg->rep_index] = reference_fg;
      collected_indices.push_back(reference_index);

      /// save the collected FeatureGroups to output
      ConsensusFeatureGroup cfg;
      std::sort(collected_indices.begin(), collected_indices.end(), greater<Size>()); // to remove indices from the largest number
      double accum_mass = .0;
      double accum_rt = .0;
      for (Size i = 0; i < collected_fgs.size(); ++i)
      {
        // save data
        cfg.fgroup_indices.push_back(collected_fgs[i]->fgroup_index);
        cfg.abundances.push_back(collected_fgs[i]->abundance);
        accum_mass += collected_fgs[i]->mass;
        accum_rt += collected_fgs[i]->apex_rt;

        // remove from fgroup_pointer
        fgroup_pointers.erase(fgroup_pointers.begin() + collected_indices[i]);
      }
      cfg.avg_mass = accum_mass / collected_fgs.size();
      cfg.avg_apex_rt = accum_rt / collected_fgs.size();
      cfg.calculate_cv();
      consensus.push_back(cfg);
    }
    consensus.shrink_to_fit();
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // Parameter handling
    //-------------------------------------------------------------
    StringList ins = getStringList_("in");
    String out = getStringOption_("out");

    MASS_TOL_ = getDoubleOption_("mass_tol");
    MASS_TOL_UNIT_ = getStringOption_("mass_tol_unit");
    RT_TOL_ = getIntOption_("rt_tol");
    QUANT_METHOD_ = getStringOption_("quant_method");
    REP_COUNT_ = ins.size();
    String consensus_input = getStringOption_("consensus_as_input");
    if (MASS_TOL_UNIT_ == "ppm")
    {
      MASS_TOL_ *= 1e-6;
    }

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    std::vector<FeatureGroup> feat_groups;
    if (consensus_input == "false")
    {
      for (Size i = 0; i < ins.size(); ++i)
      {
        OPENMS_LOG_INFO << ins[i] << " as File" << i ;
        readFLASHQuantResultFile(ins[i], feat_groups, i);
      }
    }
    else
    {
      for (Size i = 0; i < ins.size(); ++i)
      {
        OPENMS_LOG_INFO << ins[i] << " as File" << i ;
        readFLASHQuantConsensusFile(ins[i], feat_groups, i);
      }
    }

    //-------------------------------------------------------------
    // calculate consensus feature groups
    //-------------------------------------------------------------
    std::vector<ConsensusFeatureGroup> consensus;
    computeConsensusFeatureGroup(feat_groups, consensus);
    cout << "#consensus=" << consensus.size() << endl;

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    writeConsensusFeatureGroupInTsv(consensus, out, ins.size());

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TopDownConsensusFeatureGroup tool;
  return tool.main(argc, argv);
}
