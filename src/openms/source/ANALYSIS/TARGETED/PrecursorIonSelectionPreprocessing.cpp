// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/SIMULATION/DetectabilitySimulation.h>
#include <OpenMS/SIMULATION/RTSimulation.h>

#include <boost/math/distributions/normal.hpp>

using namespace std;
//#define PISP_DEBUG
#ifdef PISP_DEBUG
#include <OpenMS/SYSTEM/StopWatch.h>
#endif
//#undef PISP_DEBUG
namespace OpenMS
{
  PrecursorIonSelectionPreprocessing::PrecursorIonSelectionPreprocessing() :
    DefaultParamHandler("PrecursorIonSelectionPreprocessing"),
    f_max_(0)
  {
    defaults_.setValue("precursor_mass_tolerance", 10., "Precursor mass tolerance which is used to query the peptide database for peptides");
    defaults_.setMinFloat("precursor_mass_tolerance", 0.);
    defaults_.setValue("rt_settings:min_rt", 960., "Minimal RT in the experiment (in seconds)");
    defaults_.setMinFloat("rt_settings:min_rt", 0.);
    defaults_.setValue("rt_settings:max_rt", 3840., "Maximal RT in the experiment (in seconds)");
    defaults_.setMinFloat("rt_settings:min_rt", 1.);
    defaults_.setValue("rt_settings:rt_step_size", 30., "Time between two consecutive spectra (in seconds)");
    defaults_.setMinFloat("rt_settings:min_rt", 1.);
    defaults_.setValue("rt_settings:gauss_mean", -1.0, "mean of the gauss curve");
    defaults_.setValue("rt_settings:gauss_sigma", 3., "std of the gauss curve");
    defaults_.setValue("precursor_mass_tolerance_unit", "ppm", "Precursor mass tolerance unit.");
    defaults_.setValidStrings("precursor_mass_tolerance_unit", ListUtils::create<String>("ppm,Da"));
    defaults_.setValue("preprocessed_db_path", "", "Path where the preprocessed database should be stored");
    defaults_.setValue("preprocessed_db_pred_rt_path", "", "Path where the predicted rts of the preprocessed database should be stored");
    defaults_.setValue("preprocessed_db_pred_dt_path", "", "Path where the predicted rts of the preprocessed database should be stored");
    defaults_.setValue("max_peptides_per_run", 100000, "Number of peptides for that the pt and rt are parallely predicted.");
    defaults_.setMinInt("max_peptides_per_run", 1);
    defaults_.setValue("missed_cleavages", 1, "Number of allowed missed cleavages.");
    defaults_.setMinInt("missed_cleavages", 0);
    defaults_.setValue("taxonomy", "", "Taxonomy");
    defaults_.setValue("tmp_dir", "", "Absolute path to tmp data directory used to store files needed for rt and dt prediction.");
    defaults_.setValue("store_peptide_sequences", "false", "Flag if peptide sequences should be stored.");
    defaultsToParam_();
    updateMembers_();
  }

  PrecursorIonSelectionPreprocessing::PrecursorIonSelectionPreprocessing(const PrecursorIonSelectionPreprocessing& source) :
    DefaultParamHandler(source),
    sequences_(source.sequences_),
    prot_masses_(source.prot_masses_),
    bin_masses_(source.bin_masses_),
    f_max_(source.f_max_)
  {
    updateMembers_();
  }

  PrecursorIonSelectionPreprocessing::~PrecursorIonSelectionPreprocessing()
  {
    //????
  }

  PrecursorIonSelectionPreprocessing& PrecursorIonSelectionPreprocessing::operator=(const PrecursorIonSelectionPreprocessing& source)
  {
    if (&source != this)
    {
      DefaultParamHandler::operator=(source);
      sequences_ = source.sequences_;
      prot_masses_ = source.prot_masses_;
      bin_masses_ = source.bin_masses_;
      f_max_ = source.f_max_;
    }
    return *this;
  }

  void PrecursorIonSelectionPreprocessing::setFixedModifications(StringList& modifications)
  {
    fixed_modifications_.clear();
    for (Size i = 0; i < modifications.size(); ++i)
    {
#ifdef DEBUG_PISP
      std::cout << modifications[i] << std::endl;
#endif
      // get specification
      String spec = modifications[i].suffix('(');
#ifdef DEBUG_PISP
      std::cerr << "specificity: " << spec[0] << " " << spec << std::endl;
      std::cout << "modname: " << modifications[i].prefix(' ') << std::endl;
#endif
      if (fixed_modifications_.find(spec[0]) != fixed_modifications_.end())
      {
        fixed_modifications_[spec[0]].push_back(modifications[i].prefix(' '));
      }
      else
      {
        std::vector<String> mod_names;
        mod_names.push_back(modifications[i].prefix(' '));
        fixed_modifications_.insert(make_pair(spec[0], mod_names));
      }
    }
    if (!fixed_modifications_.empty())
      fixed_mods_ = true;
  }

  const std::map<String, std::vector<double> >& PrecursorIonSelectionPreprocessing::getProtMasses() const
  {
    return prot_masses_;
  }

  const std::vector<double>& PrecursorIonSelectionPreprocessing::getMasses(String acc) const
  {
    std::map<String, std::vector<double> >::const_iterator iter = prot_masses_.begin();
    while (iter != prot_masses_.end() && acc != iter->first)
      ++iter;
    if (iter != prot_masses_.end())
      return iter->second;
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "PrecursorIonSelectionPreprocessing: protein " + acc + " could not be found.");
    }
  }

  const std::map<String, std::vector<double> >& PrecursorIonSelectionPreprocessing::getProteinRTMap() const
  {
    return rt_prot_map_;
  }

  const std::map<String, std::vector<double> >& PrecursorIonSelectionPreprocessing::getProteinPTMap() const
  {
    return pt_prot_map_;
  }

  const std::map<String, std::vector<String> >& PrecursorIonSelectionPreprocessing::getProteinPeptideSequenceMap() const
  {
    return prot_peptide_seq_map_;
  }

  double PrecursorIonSelectionPreprocessing::getRT(String prot_id, Size peptide_index)
  {
    if (rt_prot_map_.size() > 0)
    {
      if (rt_prot_map_.find(prot_id) != rt_prot_map_.end())
      {
        if (rt_prot_map_[prot_id].size() > peptide_index)
          return rt_prot_map_[prot_id][peptide_index];
        else
          return -1;
      }
      else
        return -1;
    }
    std::cout << "rt_map is empty, no rts predicted!" << std::endl;
    return -1;
  }

  double PrecursorIonSelectionPreprocessing::getPT(String prot_id, Size peptide_index)
  {
    if (pt_prot_map_.size() > 0)
    {
      if (pt_prot_map_.find(prot_id) != pt_prot_map_.end())
      {
        if (pt_prot_map_[prot_id].size() > peptide_index)
          return pt_prot_map_[prot_id][peptide_index];
        else
          return 1;
      }
      else
        return 1;
    }
    std::cout << "pt_map is empty, no detectabilities predicted!" << std::endl;
    return 1;
  }

  double PrecursorIonSelectionPreprocessing::getWeight(double mass)
  {
    if (param_.getValue("precursor_mass_tolerance_unit") == "Da")
    {
      return (double)counter_[(Size) floor((mass - masses_[0]) / (double)param_.getValue("precursor_mass_tolerance") + 0.5)] / (double)f_max_;
    }
    else //
    {
#ifdef PISP_DEBUG
      std::cout << bin_masses_.size() << " " << mass << std::endl;
      std::cout << *(bin_masses_.begin()) << " " << *(bin_masses_.end() - 1) << std::endl;
#endif
      std::vector<double>::iterator tmp_iter = bin_masses_.begin();
      while (tmp_iter != bin_masses_.end() && *tmp_iter < mass)
      {
        ++tmp_iter;
      }
      if (tmp_iter != bin_masses_.begin())
        --tmp_iter;
#ifdef PISP_DEBUG
      if (tmp_iter != bin_masses_.end())
        std::cout << "tmp_iter " << *tmp_iter << std::endl;
      else
        std::cout << "tmp_iter am ende" << std::endl;
      std::cout << "fmax " << f_max_ << std::endl;
#endif
      if ((tmp_iter + 1) == bin_masses_.end()
         || fabs(*tmp_iter - mass) < fabs(*(tmp_iter + 1) - mass))
      {
        return (double)counter_[distance(bin_masses_.begin(), tmp_iter)] / (double)f_max_;
      }
      else
        return (double)counter_[distance(bin_masses_.begin(), tmp_iter + 1)] / (double)f_max_;
    }
  }

  void PrecursorIonSelectionPreprocessing::loadPreprocessing()
  {
    // first check if preprocessed db already exists
    String path = param_.getValue("preprocessed_db_path");

    // check if file exists
    std::ifstream test(path.c_str());
    if (test)
    {
      loadPreprocessedDB_(path);
    }
    else
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path);
    }
  }

  void PrecursorIonSelectionPreprocessing::filterTaxonomyIdentifier_(FASTAFile::FASTAEntry& entry)
  {

#ifdef PISP_DEBUG
    OPENMS_LOG_DEBUG << "Original Entry-Identifier:" << entry.identifier << std::endl;
#endif
    if (entry.identifier.hasPrefix("sp|") || entry.identifier.hasPrefix("tr|") || entry.identifier.hasPrefix("gi|"))
    {
      entry.identifier = entry.identifier.suffix(entry.identifier.size() - 3);
    }
    else if (entry.identifier.hasPrefix("IPI:"))
    {
      entry.identifier = entry.identifier.suffix(entry.identifier.size() - 4);
    }

    if (entry.identifier.has('|'))
    {
      entry.identifier = entry.identifier.prefix('|');
    }

#ifdef PISP_DEBUG
    OPENMS_LOG_DEBUG << "Processed Entry-Identifier:" << entry.identifier << std::endl;
#endif

  }

  void PrecursorIonSelectionPreprocessing::dbPreprocessing(String db_path, String rt_model_path,
                                                           String dt_model_path, bool save)
  {
#ifdef PISP_DEBUG
    std::cout << "Parameters: " << param_.getValue("preprocessed_db_path")
              << "\t" << param_.getValue("precursor_mass_tolerance")
              << " " << param_.getValue("precursor_mass_tolerance_unit")
      ///<< "\t"<<param_.getValue("rt_tolerance")
              << "\t" << param_.getValue("missed_cleavages")
              << "\t" << param_.getValue("taxonomy")
              << "\t" << param_.getValue("tmp_dir") << "---"
              << std::endl;
#endif

    FASTAFile fasta_file;
    std::vector<FASTAFile::FASTAEntry> entries;
    fasta_file.load(db_path, entries);
    ProteaseDigestion digest;
    digest.setMissedCleavages(param_.getValue("missed_cleavages"));
    std::map<String, std::vector<std::pair<String, Size> > > tmp_peptide_map;

    bool store_pep_seqs = (param_.getValue("store_peptide_sequences") == "true") ? true : false;

    // first get all protein sequences and calculate digest
    for (UInt e = 0; e < entries.size(); ++e)
    {

      // filter for taxonomy
      if (entries[e].description.toUpper().hasSubstring(((String)param_.getValue("taxonomy")).toUpper()))
      {
        // preprocess entry identifier
        filterTaxonomyIdentifier_(entries[e]);

        String& seq = entries[e].sequence;
        // check for unallowed characters
        if (seq.hasSubstring("X") || seq.hasSubstring("B") ||  seq.hasSubstring("Z"))
        {
          continue;
        }
        std::vector<double> prot_masses;
        // digest sequence
        AASequence aa_seq = AASequence::fromString(seq);
        std::vector<AASequence> vec;
        digest.digest(aa_seq, vec);

        // enter peptide sequences in map
        std::vector<AASequence>::iterator vec_iter = vec.begin();

        std::vector<String> peptide_seqs;
        for (; vec_iter != vec.end(); ++vec_iter)
        {

          // enter mod
          if (fixed_mods_)
          {
            // go through peptide sequence and check if AA is modified
            for (Size aa = 0; aa < vec_iter->size(); ++aa)
            {
              if (fixed_modifications_.find((vec_iter->toUnmodifiedString())[aa]) != fixed_modifications_.end())
              {
#ifdef DEBUG_PISP
                std::cout << "w/o Mod " << *vec_iter << " "
                          << vec_iter->getMonoWeight(Residue::Full, 1) << std::endl;
#endif
                std::vector<String>& mods = fixed_modifications_[(vec_iter->toUnmodifiedString())[aa]];
                for (Size m = 0; m < mods.size(); ++m)
                {
                  vec_iter->setModification(aa, mods[m]);
                }
#ifdef DEBUG_PISP
                std::cout << "set Mods " << *vec_iter << " "
                          << vec_iter->getMonoWeight(Residue::Full, 1) << std::endl;
#endif
              }
            }
          }

          // write peptide seq in temporary file, for rt prediction
          //seq_file << *vec_iter << "\n";
          double mass = vec_iter->getMonoWeight(Residue::Full, 1);
          prot_masses.push_back(mass);
          if (tmp_peptide_map.find(vec_iter->toUnmodifiedString()) != tmp_peptide_map.end())
          {
            tmp_peptide_map[vec_iter->toUnmodifiedString()].push_back(make_pair(entries[e].identifier,
                                                                                prot_masses.size() - 1));
          }
          else
          {
            std::vector<std::pair<String, Size> > tmp_vec;
            tmp_vec.push_back(make_pair(entries[e].identifier, prot_masses.size() - 1));
            tmp_peptide_map.insert(make_pair(vec_iter->toUnmodifiedString(), tmp_vec));
          }
          if (sequences_.count(*vec_iter) == 0) // peptide sequences are considered only once
          {
            sequences_.insert(*vec_iter);
            masses_.push_back(mass);
          }
          if (store_pep_seqs)
            peptide_seqs.push_back(vec_iter->toUnmodifiedString());
        }
        prot_masses_.insert(make_pair(entries[e].identifier, prot_masses));
        if (store_pep_seqs)
          prot_peptide_seq_map_.insert(make_pair(entries[e].identifier, peptide_seqs));
      }

    }
    entries.clear();
#ifdef PIPS_DEBUG
    std::cout << "now make the rt and pt predictions" << std::endl;
#endif
    std::set<AASequence>::iterator seq_it = sequences_.begin();
    std::vector<String> peptide_sequences;
    UInt index = 0;
    Param rt_param;
    rt_param.setValue("HPLC:model_file", rt_model_path);
    Param dt_param;
    dt_param.setValue("dt_simulation_on", "true");
    dt_param.setValue("dt_model_file", dt_model_path);
    // this is needed, as too many sequences require too much memory for the rt and dt prediction
    Size max_peptides_per_run = (Int)param_.getValue("max_peptides_per_run");
    peptide_sequences.resize(std::min(sequences_.size(), max_peptides_per_run));

    //#ifdef PIPS_DEBUG
    std::cout << sequences_.size() << " peptides for predictions." << std::endl;
    //#endif
    std::set<AASequence>::iterator seq_it_end = sequences_.end();
    std::vector<AASequence> peptide_aa_sequences;
    peptide_aa_sequences.resize(std::min(sequences_.size(), max_peptides_per_run));
    if (seq_it != seq_it_end)
      --seq_it_end;
    for (; seq_it != sequences_.end(); ++seq_it)
    {
      peptide_sequences[index] = seq_it->toUnmodifiedString();
      peptide_aa_sequences[index] = *seq_it;
      ++index;
      if (index == max_peptides_per_run || seq_it == seq_it_end)
      {
        std::cout << distance(sequences_.begin(), seq_it) << " peptides done." << std::endl;
        // now make RTPrediction using the RTSimulation class of the simulator
        RTSimulation rt_sim;
        rt_sim.setParameters(rt_param);

        std::vector<double> rts2;
        rt_sim.wrapSVM(peptide_aa_sequences, rts2);
        for (Size index2 = 0; index2 < rts2.size(); ++index2)
        {
          for (Size seq_idx = 0; seq_idx < tmp_peptide_map[peptide_sequences[index2]].size(); ++seq_idx)
          {
            String acc = tmp_peptide_map[peptide_sequences[index2]][seq_idx].first;
            Size tmp_idx = tmp_peptide_map[peptide_sequences[index2]][seq_idx].second;
            if (rt_prot_map_.find(acc) != rt_prot_map_.end())
            {
              rt_prot_map_[acc][tmp_idx] = rts2[index2];
            }
            else
            {
              std::vector<double> rt_vec(prot_masses_[acc].size());
              rt_prot_map_.insert(make_pair(acc, rt_vec));
              rt_prot_map_[acc][tmp_idx] = rts2[index2];
            }
          }
        }
        rts2.clear();
        peptide_aa_sequences.clear();
        std::cout << "Now Detectablity" << std::endl;
        // now make DTPrediction using the DetectabilitySimulation class of the simulator
        DetectabilitySimulation dt_sim;
        dt_sim.setParameters(dt_param);
        std::vector<double> labels;
        std::vector<double> detectabilities;
        dt_sim.predictDetectabilities(peptide_sequences, labels, detectabilities);

        for (Size index2 = 0; index2 < detectabilities.size(); ++index2)
        {
          for (Size seq_idx = 0; seq_idx < tmp_peptide_map[peptide_sequences[index2]].size(); ++seq_idx)
          {
            String acc = tmp_peptide_map[peptide_sequences[index2]][seq_idx].first;
            Size tmp_idx = tmp_peptide_map[peptide_sequences[index2]][seq_idx].second;
            if (pt_prot_map_.find(acc) != pt_prot_map_.end())
            {
              pt_prot_map_[acc][tmp_idx] = detectabilities[index2];
            }
            else
            {
              std::vector<double> pt_vec(prot_masses_[acc].size());
              pt_prot_map_.insert(make_pair(acc, pt_vec));
              pt_prot_map_[acc][tmp_idx] =   detectabilities[index2];
            }
          }
        }
        peptide_sequences.clear();
        Int size = std::min((int)distance(seq_it, sequences_.end()) - 1, (int)max_peptides_per_run);
        if (size > 0)
          peptide_sequences.resize(size);
        peptide_aa_sequences.resize(size);
        index = 0;
      }
    }
    peptide_sequences.clear();
    peptide_aa_sequences.clear();
    tmp_peptide_map.clear();
    sequences_.clear();
    tmp_peptide_map.clear();
#ifdef PISP_DEBUG
    std::cout << "Finished predictions!" << std::endl;
#endif
    if (masses_.empty())
    {
      std::cout << "no masses entered" << std::endl;
      return;
    }
    std::sort(masses_.begin(), masses_.end());
    // now get minimal and maximal mass and create counter_-vectors
    // count mass occurences using bins
#ifdef PISP_DEBUG
    std::cout << "min\tmax " << masses_[0] << "\t" << *(masses_.end() - 1) << std::endl;
    std::cout << "prot_masses.size() " << prot_masses_.size() << std::endl;
#endif
    // if the precursor mass tolerance is given in Da
    // we have equidistant bins
    if (param_.getValue("precursor_mass_tolerance_unit") == "Da")
    {
      counter_.resize((UInt)(ceil((*(masses_.end() - 1)) - masses_[0]) / (double)param_.getValue("precursor_mass_tolerance")) + 2, 0);
      for (UInt i = 0; i < masses_.size(); ++i)
      {
        // get bin index
        double tmp = (masses_[i] - masses_[0]) / (double)param_.getValue("precursor_mass_tolerance");
        ++counter_[(Size) ceil(tmp)];
      }
      UInt max = 0;
      for (UInt i = 0; i < counter_.size(); ++i)
      {
        if (counter_[i] >  max)
          max = counter_[i];
      }
      // store maximal frequency
      f_max_ = max;
    }
    else // with ppm, we have bins of increasing size and need to save the bin limits
    {
      bin_masses_.clear();
      counter_.clear();
      // so first we calculate the boundings of the bins
      double curr_mass = masses_[0];
#ifdef PISP_DEBUG
      std::cout << "min_max_curr " << masses_[0] << " " << *(masses_.end() - 1) << " " << curr_mass << std::endl;
#endif
      UInt size = 0;
      while (curr_mass < *(masses_.end() - 1))
      {
        /// store the lower bound of the current bin
        bin_masses_.push_back(curr_mass);
        ++size;
        /// calculate lower bound for next bin
        curr_mass = curr_mass + curr_mass * (double)param_.getValue("precursor_mass_tolerance") / 1e06;
      }
#ifdef PISP_DEBUG
      std::cout << "bin_masses_.size() " <<  bin_masses_.size() << " " << size << std::endl;
#endif
      counter_.resize(size, 0);
#ifdef PISP_DEBUG
      std::cout << "masses_.size() " << masses_.size()
                << " counter_.size() " << counter_.size()
                << " bin_masses_.size() " << bin_masses_.size() << std::endl;

#endif
      std::vector<double>::iterator old_begin = bin_masses_.begin();
      // then we put the peptide masses into the right bins
      for (UInt i = 0; i < masses_.size(); ++i)
      {
        // #ifdef PISP_DEBUG
        //   std::cout << "i "<<i << std::endl;
        //   StopWatch timer;
        //   timer.start();
        // #endif
        std::vector<double>::iterator tmp_iter = old_begin;

        // #ifdef PISP_DEBUG
        //   timer.start();
        //   std::cout << "old_begin "<<*old_begin << " masses_[i] "<<masses_[i]<<std::endl;
        //   if(old_begin != bin_masses_.begin()) std::cout << "old_begin-1 "<<*(old_begin-1)<<std::endl;
        // #endif
        while (tmp_iter != bin_masses_.end() &&  *tmp_iter < masses_[i])
        {
          ++tmp_iter;
        }
        // #ifdef PISP_DEBUG
        //  timer.stop();
        //  std::cout << "while schleife\t";
        //  std::cout << timer.getCPUTime()<<"\t";
        //  timer.reset();
        // #endif


        old_begin = tmp_iter;
        if (tmp_iter == bin_masses_.end())
        {
          ++counter_[distance(bin_masses_.begin(), tmp_iter - 1)];
        }
        else if ((tmp_iter + 1) == bin_masses_.end() // last entry
                || fabs(*tmp_iter - masses_[i]) < fabs(*(tmp_iter + 1) - masses_[i])) // or closer to left entry
        {
          // increase counter for bin to the left
          ++counter_[distance(bin_masses_.begin(), tmp_iter)];
        }
        else
          ++counter_[distance(bin_masses_.begin(), tmp_iter + 1)]; // increase right counter
        // #ifdef PISP_DEBUG
        //                      timer.stop();
        //                      std::cout << timer.getCPUTime ()<<std::endl;
        // #endif
      }
      // determine maximal frequency for normalization
      UInt max = 0;
      for (UInt i = 0; i < counter_.size(); ++i)
      {
        if (counter_[i] >  max)
          max = counter_[i];
      }
      f_max_ = max;
#ifdef PISP_DEBUG
      std::cout << "f_max_ " << f_max_ << std::endl;
#endif

    }
    if (save)
    {
      savePreprocessedDBWithRT_(db_path, (String)param_.getValue("preprocessed_db_path"));
    }
  }

  void PrecursorIonSelectionPreprocessing::dbPreprocessing(String db_path, bool save)
  {
#ifdef PISP_DEBUG
    std::cout << "Parameters: " << param_.getValue("preprocessed_db_path")
              << "\t" << param_.getValue("precursor_mass_tolerance")
              << " " << param_.getValue("precursor_mass_tolerance_unit")
              << "\t" << param_.getValue("missed_cleavages")
              << "\t" << param_.getValue("taxonomy") << std::endl;
#endif

    FASTAFile fasta_file;
    std::vector<FASTAFile::FASTAEntry> entries;
    fasta_file.load(db_path, entries);
    ProteaseDigestion digest;
    digest.setMissedCleavages(param_.getValue("missed_cleavages"));

    // first get all protein sequences and calculate digest
    for (UInt e = 0; e < entries.size(); ++e)
    {
      // filter for taxonomy
      if (entries[e].description.toUpper().hasSubstring(((String)param_.getValue("taxonomy")).toUpper()))
      {
        // preprocess entry identifier
        filterTaxonomyIdentifier_(entries[e]);

        String& seq = entries[e].sequence;
        // check for unallowed characters
        if (seq.hasSubstring("X") || seq.hasSubstring("B") ||  seq.hasSubstring("Z"))
        {
          continue;
        }
        std::vector<double> prot_masses;
        // digest sequence
        AASequence aa_seq = AASequence::fromString(seq);
        std::vector<AASequence> vec;
        digest.digest(aa_seq, vec);

        // enter peptide sequences in map
        std::vector<AASequence>::iterator vec_iter = vec.begin();
        for (; vec_iter != vec.end(); ++vec_iter)
        {
          // enter mod
          if (fixed_mods_)
          {
            // go through peptide sequence and check if AA is modified
            for (Size aa = 0; aa < vec_iter->size(); ++aa)
            {
              if (fixed_modifications_.find((vec_iter->toUnmodifiedString())[aa]) != fixed_modifications_.end())
              {
#ifdef DEBUG_PISP
                std::cout << "w/o Mod " << *vec_iter << " "
                          << vec_iter->getMonoWeight(Residue::Full, 1) << std::endl;
#endif
                std::vector<String>& mods = fixed_modifications_[(vec_iter->toUnmodifiedString())[aa]];
                for (Size m = 0; m < mods.size(); ++m)
                {
                  vec_iter->setModification(aa, mods[m]);
                }
#ifdef DEBUG_PISP
                std::cout << "set Mods " << *vec_iter << " "
                          << vec_iter->getMonoWeight(Residue::Full, 1) << std::endl;
#endif
              }
            }
          }

          double mass = vec_iter->getMonoWeight(Residue::Full, 1);
          prot_masses.push_back(mass);
          if (sequences_.count(*vec_iter) == 0) // peptide sequences are considered only once
          {
            sequences_.insert(*vec_iter);
            masses_.push_back(mass);
          }
        }
        prot_masses_.insert(make_pair(entries[e].identifier, prot_masses));
      }

    }

    if (masses_.empty())
    {
      std::cout << "no masses entered" << std::endl;
      return;
    }
    std::sort(masses_.begin(), masses_.end());
    // now get minimal and maximal mass and create counter_-vectors
    // count mass occurrences using bins
#ifdef PISP_DEBUG
    std::cout << "min\tmax " << masses_[0] << "\t" << *(masses_.end() - 1) << std::endl;
    std::cout << "prot_masses.size() " << prot_masses_.size() << std::endl;
#endif
    // if the precursor mass tolerance is given in Da
    // we have equidistant bins
    if (param_.getValue("precursor_mass_tolerance_unit") == "Da")
    {
      counter_.resize((UInt)(ceil((*(masses_.end() - 1)) - masses_[0]) / (double)param_.getValue("precursor_mass_tolerance")) + 2, 0);
      for (UInt i = 0; i < masses_.size(); ++i)
      {
        // get bin index
        double tmp = (masses_[i] - masses_[0]) / (double)param_.getValue("precursor_mass_tolerance");
        ++counter_[(Size) ceil(tmp)];
      }
      UInt max = 0;
      for (UInt i = 0; i < counter_.size(); ++i)
      {
        if (counter_[i] >  max)
          max = counter_[i];
      }
      // store maximal frequency
      f_max_ = max;
    }
    else // with ppm, we have bins of increasing size and need to save the bin limits
    {
      bin_masses_.clear();
      counter_.clear();
      // so first we calculate the boundings of the bins
      double curr_mass = masses_[0];
#ifdef PISP_DEBUG
      std::cout << "min_max_curr " << masses_[0] << " " << *(masses_.end() - 1) << " " << curr_mass << std::endl;
#endif
      UInt size = 0;
      while (curr_mass < *(masses_.end() - 1))
      {
        /// store the lower bound of the current bin
        bin_masses_.push_back(curr_mass);
        ++size;
        /// calculate lower bound for next bin
        curr_mass = curr_mass + curr_mass * (double)param_.getValue("precursor_mass_tolerance") / 1e06;
      }
#ifdef PISP_DEBUG
      std::cout << "bin_masses_.size() " <<  bin_masses_.size() << " " << size << std::endl;
#endif
      if (bin_masses_.empty())
      {
        throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0);
      }
      counter_.resize(size, 0);
#ifdef PISP_DEBUG
      std::cout << "masses_.size() " << masses_.size()
                << " counter_.size() " << counter_.size()
                << " bin_masses_.size() " << bin_masses_.size() << std::endl;

#endif
      std::vector<double>::iterator old_begin = bin_masses_.begin();
      // then we put the peptide masses into the right bins
      for (UInt i = 0; i < masses_.size(); ++i)
      {
#ifdef PISP_DEBUG
        std::cout << "i " << i << std::endl;
        StopWatch timer;
        timer.start();
#endif
        std::vector<double>::iterator tmp_iter = old_begin;

#ifdef PISP_DEBUG
        timer.start();
        std::cout << "old_begin " << *old_begin << " masses_[i] " << masses_[i] << std::endl;
        if (old_begin != bin_masses_.begin())
          std::cout << "old_begin-1 " << *(old_begin - 1) << std::endl;
#endif
        while (tmp_iter != bin_masses_.end() &&  *tmp_iter < masses_[i])
        {
          ++tmp_iter;
        }
#ifdef PISP_DEBUG
        timer.stop();
        std::cout << "while schleife\t";
        std::cout << timer.getCPUTime() << "\t";
        timer.reset();
#endif


        old_begin = tmp_iter;
        if (tmp_iter == bin_masses_.end())
        {
          ++counter_[distance(bin_masses_.begin(), tmp_iter - 1)];
        }
        else if ((tmp_iter + 1) == bin_masses_.end() // last entry
                || fabs(*tmp_iter - masses_[i]) < fabs(*(tmp_iter + 1) - masses_[i])) // or closer to left entry
        {
          // increase counter for bin to the left
          ++counter_[distance(bin_masses_.begin(), tmp_iter)];
        }
        else
          ++counter_[distance(bin_masses_.begin(), tmp_iter + 1)]; // increase right counter
#ifdef PISP_DEBUG
        timer.stop();
        std::cout << timer.getCPUTime() << std::endl;
#endif
      }
      // determine maximal frequency for normalization
      UInt max = 0;
      for (UInt i = 0; i < counter_.size(); ++i)
      {
        if (counter_[i] >  max)
          max = counter_[i];
      }
      f_max_ = max;
#ifdef PISP_DEBUG
      std::cout << "f_max_ " << f_max_ << std::endl;
#endif

    }
    if (save)
    {
      savePreprocessedDB_(db_path, (String)param_.getValue("preprocessed_db_path"));
    }

  }

  void PrecursorIonSelectionPreprocessing::savePreprocessedDB_(String db_path, String path)
  {
    std::ofstream out(path.c_str());
    out.precision(10);
    if (!out)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path);
    }
    // header: db_name  precursor_mass_tolerance precursor_mass_tolerance_unit  taxonomy);
    Size pos1 = db_path.rfind("/") + 1;
    // get db-name
    Size pos2 = db_path.rfind(".");
    String db_name = db_path.substr(pos1, pos2 - pos1);
    out << db_name << "\t" << param_.getValue("precursor_mass_tolerance")  << "\t"
        << param_.getValue("precursor_mass_tolerance_unit")
        << "\t" << (String)param_.getValue("taxonomy");
    // first save protein_masses_map
    out << prot_masses_.size() << std::endl;
#ifdef PISP_DEBUG
    std::cout << prot_masses_.size() << " " << counter_.size() << " " << bin_masses_.size() << std::endl;
#endif
    std::map<String, std::vector<double> >::iterator pm_iter = prot_masses_.begin();
    for (; pm_iter != prot_masses_.end(); ++pm_iter)
    {
      out << pm_iter->second.size() << "\t" << pm_iter->first;
      for (UInt i = 0; i < pm_iter->second.size(); ++i)
      {
        out << "\t" << pm_iter->second[i];
      }
      out << "\n";
    }
    //  out.close();
#ifdef PISP_DEBUG
    std::cout << (path + "_counter") << std::endl;
    std::cout << "counter: "  << counter_.size() << "\t" << masses_[0] << "\t" << *(masses_.end() - 1) << "\n";
#endif
    //      // now save counter
    //      std::ofstream out2((path + "_counter").c_str());
    //      if (!out2)
    //      {
    //        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path+ "_counter");
    //      }
    //      out2.precision(10);

    // header: number of entries, min, max
    out << "###\n";
    out << counter_.size() << "\t" << masses_[0] << "\t" << *(masses_.end() - 1) << "\n";
    for (UInt i = 0; i < counter_.size(); ++i)
    {
      out << counter_[i] << "\t";
    }
    out << "\n";
    if (param_.getValue("precursor_mass_tolerance_unit") == "ppm")
    {
#ifdef PISP_DEBUG
      std::cout << (path + "_bin_masses") << std::endl;
#endif
      //                // now save bin_masses
      //                std::ofstream out3((path + "_bin_masses").c_str());
      //                if (!out3)
      //                {
      //                  throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path+ "_bin_masses");
      //                }

      //                out3.precision(10);
      // header: number of entries
      out << "###\n";
      out << bin_masses_.size() << "\n";
      for (UInt i = 0; i < bin_masses_.size(); ++i)
      {
        out << bin_masses_[i] << "\n";
      }
    }

  }

  void PrecursorIonSelectionPreprocessing::savePreprocessedDBWithRT_(String db_path, String path)
  {
    std::ofstream out(path.c_str());
    out.precision(10);
    if (!out)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path);
    }
    // header: db_name  precursor_mass_tolerance precursor_mass_tolerance_unit  taxonomy);
    Size pos1 = db_path.rfind("/") + 1;
    // get db-name
    Size pos2 = db_path.rfind(".");
    String db_name = db_path.substr(pos1, pos2 - pos1);
    out << db_name << "\t" << param_.getValue("precursor_mass_tolerance")  << "\t"
        << param_.getValue("precursor_mass_tolerance_unit")
        << "\t" << (String)param_.getValue("taxonomy");
    // first save protein_masses_map
    out << prot_masses_.size() << std::endl;
#ifdef PISP_DEBUG
    std::cout << prot_masses_.size() << " " << counter_.size() << " " << bin_masses_.size() << std::endl;
#endif
    FASTAFile fasta_file;
    std::vector<FASTAFile::FASTAEntry> entries;
    fasta_file.load(db_path, entries);
    ProteaseDigestion digest;
    digest.setMissedCleavages(param_.getValue("missed_cleavages"));

    // first get all protein sequences and calculate digest
    for (UInt e = 0; e < entries.size(); ++e)
    {
      // filter for taxonomy
      if (entries[e].description.toUpper().hasSubstring(((String)param_.getValue("taxonomy")).toUpper()))
      {
        // preprocess entry identifier
        filterTaxonomyIdentifier_(entries[e]);
#ifdef PISP_DEBUG
        std::cout << entries[e].identifier << std::endl;
#endif
        String& seq = entries[e].sequence;
        // check for unallowed characters
        if (seq.hasSubstring("X") || seq.hasSubstring("B") ||  seq.hasSubstring("Z"))
        {
          continue;
        }

        // digest sequence
        AASequence aa_seq = AASequence::fromString(seq);
        std::vector<AASequence> vec;
        digest.digest(aa_seq, vec);
        out << vec.size() << "\t" << entries[e].identifier;
        // enter peptide sequences in map
        std::vector<AASequence>::iterator vec_iter = vec.begin();
        for (; vec_iter != vec.end(); ++vec_iter)
        {
          // write peptide seq in temporary file, for rt prediction
          //seq_file << *vec_iter << "\n";
          double mass = vec_iter->getMonoWeight(Residue::Full, 1);
          // out : masse, rt, pt
          out << "\t" << mass << "," << getRT(entries[e].identifier, distance(vec.begin(), vec_iter))
              << "," << getPT(entries[e].identifier, distance(vec.begin(), vec_iter));
        }
        out << "\n";
      }

    }

#ifdef PISP_DEBUG
    std::cout << (path + "_counter") << std::endl;
    std::cout << "counter: "  << counter_.size() << "\t" << masses_[0] << "\t" << *(masses_.end() - 1) << "\n";
#endif
    //      // now save counter
    //      std::ofstream out2((path + "_counter").c_str());
    //      if (!out2)
    //      {
    //        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path+ "_counter");
    //      }
    //      out2.precision(10);

    // header: number of entries, min, max
    out << "###\n";
    out << counter_.size() << "\t" << masses_[0] << "\t" << *(masses_.end() - 1) << "\n";
    for (UInt i = 0; i < counter_.size(); ++i)
    {
      out << counter_[i] << "\t";
    }
    out << "\n";
    if (param_.getValue("precursor_mass_tolerance_unit") == "ppm")
    {
#ifdef PISP_DEBUG
      std::cout << (path + "_bin_masses") << std::endl;
#endif
      //                // now save bin_masses
      //                std::ofstream out3((path + "_bin_masses").c_str());
      //                if (!out3)
      //                {
      //                  throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path+ "_bin_masses");
      //                }

      //                out3.precision(10);
      // header: number of entries
      out << "###\n";
      out << bin_masses_.size() << "\n";
      for (UInt i = 0; i < bin_masses_.size(); ++i)
      {
        out << bin_masses_[i] << "\n";
      }
    }

  }

  void PrecursorIonSelectionPreprocessing::loadPreprocessedDB_(String path)
  {
    // first get protein_masses_map
    TextFile file;
    file.load(path, true);
    //#ifdef PISP_DEBUG
    std::cout << "load " << path << std::endl;
    //#endif
    TextFile::ConstIterator iter = file.begin();
    ++iter;
    for (; iter != file.end() && !iter->hasPrefix("###"); ++iter)
    {
      std::vector<String> parts;
      iter->split('\t', parts);
      std::vector<double> masses;
      masses.reserve(parts[0].toInt());
      std::vector<String> line_parts;
      std::vector<double> rts;
      std::vector<double> pts;
      for (UInt i = 2; i < parts.size(); ++i)
      {
        // check if rts are stored also
        if (parts[i].hasSubstring(","))
        {
          parts[i].split(',', line_parts);
          masses.push_back(line_parts[0].toDouble());
          if (line_parts.size() > 1)
          {
            rts.push_back(line_parts[1].toDouble());
            if (line_parts.size() == 3)
              pts.push_back(line_parts[2].toDouble());
          }
        }
        else
          masses.push_back(parts[i].toDouble());
      }
      if (parts[1].hasSubstring("."))
        parts[1] = parts[1].prefix(11);
      prot_masses_.insert(make_pair(parts[1], masses));

      if (!rts.empty())
      {
        rt_prot_map_.insert(make_pair(parts[1], rts));
      }

      if (!pts.empty())
      {
        pt_prot_map_.insert(make_pair(parts[1], pts));
      }

#ifdef PISP_DEBUG
      std::cout << parts[1] << " " << masses.size() << std::endl;
#endif
    }
#ifdef PISP_DEBUG
    std::cout << "loaded" << std::endl;
#endif
    //    if(pt_prot_map)
    //      file.clear();
    //      // now get the counter_ vector
    //      TextFile file2;
    //      file2.load(path+"_counter",true);
    //      iter = file2.begin();
    ++iter;
    std::vector<String> header;
    iter->split('\t', header);
    // store minimal mass
    masses_.push_back(header[1].toFloat());

    f_max_ = 0;

    ++iter;
    std::vector<String> freqs;
    iter->split('\t', freqs);

    for (std::vector<String>::const_iterator f_iter = freqs.begin(); f_iter != freqs.end(); ++f_iter)
    {
      counter_.push_back(f_iter->toInt());
      if ((UInt)f_iter->toInt() > f_max_)
        f_max_ = f_iter->toInt();
    }
    ++iter;

    if (param_.getValue("precursor_mass_tolerance_unit") == "ppm")
    {
      if (iter == file.end())
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "ppm is used as precursor_mass_tolerance_unit, which requires the file "
                                          + path + "_bin_masses" + ", that could not be found.");
      else if (iter->hasPrefix("###"))
      {
        ++iter;
#ifdef PISP_DEBUG
        std::cout << "load " << path << "_bin_masses" << std::endl;
#endif
        bin_masses_.reserve(iter->toInt());
        ++iter;
        for (; iter != file.end(); ++iter)
        {
          bin_masses_.push_back(iter->toDouble());
        }

      }
      else
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "ppm is used as precursor_mass_tolerance_unit, which requires the file "
                                          + path + "_bin_masses" + ", that could not be found.");
    }


    //      // load rt predictions
    //      if(param_.getValue("preprocessing:preprocessed_db_pred_rt_path")!= ""
    //          &&File().exists(param_.getValue("preprocessing:preprocessed_db_pred_rt_path")))
    //      {
    //        rt_map_.clear();
    //        TextFile rt_pred_file;
    //        rt_pred_file.load(param_.getValue("preprocessing:preprocessed_db_pred_rt_path"),true);
    //        TextFile::Iterator text_it = rt_pred_file.begin();
    //        for(;text_it != rt_pred_file.end();++text_it)
    //        {
    //          std::vector<String> parts;
    //          text_it->split(' ',parts);
    //          rt_map_.insert(make_pair(parts[0],parts[1].toDouble()));
    //        }
    //      }
    //      std::cout <<"rt_map.size: "<<rt_map_.size()<<std::endl;

  }

  double PrecursorIonSelectionPreprocessing::getRTProbability_(double min_obs_rt, double max_obs_rt, double theo_rt)
  {
    // first adapt gaussian RT error distribution to a normal distribution with \mu = 0
    Int theo_scan = getScanNumber_(theo_rt);
    if (theo_scan == -1)
      return 0.;

    double obs_scan_begin = getScanNumber_(min_obs_rt);
    if (obs_scan_begin != 0)
      obs_scan_begin -= 1;
    double obs_scan_end = getScanNumber_(max_obs_rt) + 1;

    if (obs_scan_begin == -1 || obs_scan_end == -1)
    {
      std::cerr << "Probably an error occured during RTProb-calc: scan = -1: "
                << obs_scan_begin << " " << obs_scan_end << std::endl;
      return 0.;
    }

    obs_scan_begin -= mu_;
    obs_scan_end -= mu_;
    double x1 = double(theo_scan) - obs_scan_end;
    double x2 = double(theo_scan) - obs_scan_begin;

    // boost::math::cdf computes the cumulative probs up to x (i.e. the area under the curve)
    boost::math::normal dist(0.0, sigma_);
    double prob = boost::math::cdf(dist, x1) - boost::math::cdf(dist, x2);
    if (x2 > x1)
      prob = boost::math::cdf(dist, x2) - boost::math::cdf(dist, x1);
    if ((prob < 0.) || (obs_scan_begin == obs_scan_end))
    {
      std::cerr << min_obs_rt << " " << obs_scan_begin << " " << max_obs_rt << " " << obs_scan_end << " "
                << theo_rt << " " << theo_scan << " " << mu_ << " " << x1 << " " << x2 << " " << prob << std::endl;
      if (x2 > x1)
        std::cerr <<  boost::math::cdf(dist, x2) << " - " << boost::math::cdf(dist, x1) << std::endl;
      else
        std::cerr << boost::math::cdf(dist, x1) << " - " << boost::math::cdf(dist, x2) << std::endl;
    }
    return prob;
  }

  Int PrecursorIonSelectionPreprocessing::getScanNumber_(double rt)
  {
    double min_rt = param_.getValue("rt_settings:min_rt");
    double max_rt = param_.getValue("rt_settings:max_rt");
    double rt_step_size = param_.getValue("rt_settings:rt_step_size");

    if (rt > max_rt || rt < min_rt)
      return -1;

    Int scan = (Int)floor((rt - min_rt) / rt_step_size);
    return scan;
  }

  void PrecursorIonSelectionPreprocessing::setGaussianParameters(double mu, double sigma)
  {
    sigma_ = sigma;
    mu_ = mu;
  }

  void PrecursorIonSelectionPreprocessing::updateMembers_()
  {
    sigma_ = param_.getValue("rt_settings:gauss_sigma");
    mu_ = param_.getValue("rt_settings:gauss_mean");
  }

  double PrecursorIonSelectionPreprocessing::getRTProbability(String prot_id, Size peptide_index, Feature& feature)
  {
    double theo_rt = 0.;
    if (rt_prot_map_.size() > 0)
    {
      if (rt_prot_map_.find(prot_id) != rt_prot_map_.end())
      {
        if (rt_prot_map_[prot_id].size() > peptide_index)
          theo_rt = rt_prot_map_[prot_id][peptide_index];
      }
    }
    if (theo_rt == 0.)
    {
      if (rt_prot_map_.find(prot_id) == rt_prot_map_.end())
      {
        std::cerr << " prot_id not in map " << prot_id << std::endl;
      }
      else
      {
        std::cerr << "protein in map, but " << peptide_index << " " << rt_prot_map_[prot_id].size() << std::endl;
      }
      std::cerr << "rt_map is empty, no rts predicted!" << std::endl;
    }
    double rt_begin = feature.getConvexHull().getBoundingBox().minPosition()[0];
    double rt_end =   feature.getConvexHull().getBoundingBox().maxPosition()[0];
    return getRTProbability_(rt_begin, rt_end, theo_rt);
  }

  double PrecursorIonSelectionPreprocessing::getRTProbability(double pred_rt, Feature& feature)
  {
    double rt_begin = feature.getConvexHull().getBoundingBox().minPosition()[0];
    double rt_end =   feature.getConvexHull().getBoundingBox().maxPosition()[0];

    return getRTProbability_(rt_begin, rt_end, pred_rt);
  }

} //namespace
