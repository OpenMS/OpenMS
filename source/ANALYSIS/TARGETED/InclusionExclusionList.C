// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/TARGETED/InclusionExclusionList.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/SIMULATION/RTSimulation.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>

#include <fstream>

namespace OpenMS
{
  InclusionExclusionList::InclusionExclusionList() :
    DefaultParamHandler("InclusionExclusionList")
  {

    defaults_.setValue("missed_cleavages", 0, "Number of missed cleavages used for protein digestion.\n");
    defaults_.setValue("RT:unit", "minutes", "Create lists with units as seconds instead of minutes");
    defaults_.setValidStrings("RT:unit", StringList::create("minutes,seconds"));
    defaults_.setValue("RT:use_relative", "true", "Use relative RT window, which depends on RT of precursor.");
    defaults_.setValidStrings("RT:use_relative", StringList::create("true,false"));
    defaults_.setValue("RT:window_relative", 0.05, "[for RT:use_relative == true] The relative factor X for the RT exclusion window, e.g. the window is calculated as [rt - rt*X, rt + rt*X].");
    defaults_.setMinFloat("RT:window_relative", 0.0);
    defaults_.setMaxFloat("RT:window_relative", 10.0);
    defaults_.setValue("RT:window_absolute", 90.0, "[for RT:use_relative == false] The absolute value X for the RT exclusion window in [sec], e.g. the window is calculated as [rt - X, rt + X].");
    defaults_.setMinFloat("RT:window_absolute", 0.0);
    defaults_.setValue("merge:mz_tol", 10.0, "Two inclusion/exclusion windows are merged when they (almost) overlap in RT (see 'rt_tol') and are close in m/z by this tolerance. Unit of this is defined in 'mz_tol_unit'.");
    defaults_.setMinFloat("merge:mz_tol", 0.0);
    defaults_.setValue("merge:mz_tol_unit", "ppm", "Unit of 'mz_tol'");
    defaults_.setValidStrings("merge:mz_tol_unit", StringList::create("ppm,Da"));
    defaults_.setValue("merge:rt_tol", 1.1, "Maximal RT delta (in seconds) which would allow two windows in RT to overlap (which causes merging the windows). Two inclusion/exclusion windows are merged when they (almost) overlap in RT and are close in m/z by this tolerance (see 'mz_tol'). Unit of this param is [seconds].");
    defaults_.setMinFloat("merge:rt_tol", 0.0);


    defaultsToParam_();
  }

  void InclusionExclusionList::mergeOverlappingWindows_(WindowList & list) const
  {
    bool rt_in_seconds = (param_.getValue("RT:unit") == "seconds");

    std::vector<BinaryTreeNode> tree;
    // local scope to save memory - we do not need the clustering stuff later
    {
      DoubleReal min_to_s_factor = rt_in_seconds ? 1.0 : (1.0 / 60.0);

      bool mz_as_ppm = (param_.getValue("merge:mz_tol_unit") == "ppm");

      WindowDistance_ llc(DoubleReal(param_.getValue("merge:rt_tol")) * min_to_s_factor, DoubleReal(param_.getValue("merge:mz_tol")), mz_as_ppm);
      SingleLinkage sl;
      DistanceMatrix<Real> dist;       // will be filled
      ClusterHierarchical ch;

      //ch.setThreshold(0.99);
      // clustering ; threshold is implicitly at 1.0, i.e. distances of 1.0 (== similiarity 0) will not be clustered
      ch.cluster<IEWindow, WindowDistance_>(list, llc, sl, tree, dist);
    }

    // extract the clusters
    ClusterAnalyzer ca;
    std::vector<std::vector<Size> > clusters;
    // count number of real tree nodes (not the -1 ones):
    Size node_count = 0;
    for (Size ii = 0; ii < tree.size(); ++ii)
    {
      if (tree[ii].distance >= 1)
        tree[ii].distance = -1;                               // manually set to disconnect, as SingleLinkage does not support it
      if (tree[ii].distance != -1)
        ++node_count;
    }
    ca.cut(list.size() - node_count, tree, clusters);


    WindowList list_new;

    Map<Size, Size> cluster_sizes;

    for (Size i_outer = 0; i_outer < clusters.size(); ++i_outer)
    {
      // for each cluster: create one new entry
      IEWindow w_new = list[clusters[i_outer][0]];   // init with 0th element
      // add all other elements
      for (Size i_inner = 1; i_inner < clusters[i_outer].size(); ++i_inner)
      {
        Size cl_index = clusters[i_outer][i_inner];
        w_new.MZ_ += list[cl_index].MZ_;
        w_new.RTmax_ = std::max(w_new.RTmax_, list[cl_index].RTmax_); // expand RT range
        w_new.RTmin_ = std::min(w_new.RTmin_, list[cl_index].RTmin_);
      }
      w_new.MZ_ /= clusters[i_outer].size(); // average m/z value
      list_new.push_back(w_new);

      ++cluster_sizes[clusters[i_outer].size()];   // for stats
    }

    LOG_INFO << "Clustered overlapping windows\nCluster sizes:\n";
    for (Map<Size, Size>::const_iterator it = cluster_sizes.begin(); it != cluster_sizes.end(); ++it)
    {
      LOG_INFO << "  size " << it->first << ": " << it->second << "x\n";
    }
    LOG_INFO << " --> Window count before: " << list.size() << "\n"
             << "     Window count after : " << list_new.size() << "\n";

    // replace with clustered version
    list = list_new;
  }

//   void InclusionExclusionList::loadTargets(FeatureMap<>& map, std::vector<IncludeExcludeTarget>& targets,TargetedExperiment& exp)
//   {

//   }

//   void InclusionExclusionList::loadTargets(std::vector<FASTAFile::FASTAEntry>& fasta_entries, std::vector<IncludeExcludeTarget>& targets,TargetedExperiment& exp, Size missed_cleavages)
//   {

//   }

  void InclusionExclusionList::writeTargets(const std::vector<FASTAFile::FASTAEntry> & fasta_entries,
                                            const String & out_path,
                                            const IntList & charges,
                                            const String rt_model_path)
  {
    WindowList result;

    EnzymaticDigestion digest;

    digest.setMissedCleavages(param_.getValue("missed_cleavages"));

    SimRandomNumberGenerator rnd_gen;
    RTSimulation rt_sim(rnd_gen);
    Param rt_param;
    rt_param.setValue("HPLC:model_file", rt_model_path);
    rt_sim.setParameters(rt_param);
    std::vector<AASequence> pep_seqs;
    std::vector<FASTAFile::FASTAEntry>::const_iterator entry_iter = fasta_entries.begin();
    for (; entry_iter != fasta_entries.end(); ++entry_iter)
    {
      // digest sequence
      AASequence aa_seq(entry_iter->sequence);
      std::vector<AASequence> vec;
      digest.digest(aa_seq, vec);

      // copy
      pep_seqs.insert(pep_seqs.begin(), vec.begin(), vec.end());

      // TODO: enter modifications

      // // enter mod
      // if(fixed_mods_)
      //    {
//                              // go through peptide sequence and check if AA is modified
//                              for(Size aa = 0; aa < vec_iter->size();++aa)
//                                  {
//                                      if(fixed_modifications_.find((vec_iter->toUnmodifiedString())[aa])!= fixed_modifications_.end())
//                                          {
// #ifdef DEBUG_PISP
//                                              std::cout << "w/o Mod "<<*vec_iter<<" "
//                                                                  <<vec_iter->getMonoWeight(Residue::Full,1)<<std::endl;
// #endif
//                                              std::vector<String> & mods = fixed_modifications_[(vec_iter->toUnmodifiedString())[aa]];
//                                              for(Size m = 0; m < mods.size();++m)
//                                                  {
//                                                      vec_iter->setModification(aa,mods[m]);
//                                                  }
// #ifdef DEBUG_PISP
//                                              std::cout << "set Mods "<<*vec_iter<<" "
//                                                                  <<vec_iter->getMonoWeight(Residue::Full,1)<<std::endl;
// #endif
//                                          }
//                                  }
//            }

    }
    std::vector<DoubleReal> rts;
    rt_sim.wrapSVM(pep_seqs, rts);
    DoubleReal min_to_s_factor = (param_.getValue("RT:unit") == "seconds") ? 1.0 : (1.0 / 60.0);

    bool relative_rt = (param_.getValue("RT:use_relative") == "true");

    DoubleReal rel_rt_window_size = param_.getValue("RT:window_relative");
    DoubleReal abs_rt_window_size = param_.getValue("RT:window_absolute");

    for (Size i = 0; i < pep_seqs.size(); ++i)
    {
      for (Size c = 0; c < charges.size(); ++c)
      {
        // calculate exclusion window
        DoubleReal mz = pep_seqs[i].getMonoWeight(Residue::Full, charges[c]) / (DoubleReal)charges[c];
        DoubleReal rt_start = std::max(0.0, relative_rt ? (rts[i] - rel_rt_window_size * rts[i]) : rts[i] - abs_rt_window_size);
        DoubleReal rt_stop =                relative_rt ? (rts[i] + rel_rt_window_size * rts[i]) : rts[i] + abs_rt_window_size;

        rt_start *= min_to_s_factor;
        rt_stop *= min_to_s_factor;

        result.push_back(IEWindow(rt_start, rt_stop, mz));
      }
    }

    mergeOverlappingWindows_(result);
    writeToFile_(out_path, result);
  }

  void InclusionExclusionList::writeTargets(const FeatureMap<> & map,
                                            const String & out_path)
  {
    WindowList result;

    bool relative_rt = (param_.getValue("RT:use_relative") == "true");

    DoubleReal rel_rt_window_size = param_.getValue("RT:window_relative");
    DoubleReal abs_rt_window_size = param_.getValue("RT:window_absolute");

    DoubleReal min_to_s_factor = (param_.getValue("RT:unit") == "seconds") ? 1.0 : (1.0 / 60.0);
    for (Size f = 0; f < map.size(); ++f)
    {
      DoubleReal rt_start = std::max(0.0, relative_rt ? (map[f].getRT() - map[f].getRT() * rel_rt_window_size) : map[f].getRT() - abs_rt_window_size);
      DoubleReal rt_stop =                relative_rt ? (map[f].getRT() + map[f].getRT() * rel_rt_window_size) : map[f].getRT() + abs_rt_window_size;

      rt_start *= min_to_s_factor;
      rt_stop *= min_to_s_factor;

      result.push_back(IEWindow(rt_start, rt_stop, map[f].getMZ()));
    }

    mergeOverlappingWindows_(result);
    writeToFile_(out_path, result);
  }

  void InclusionExclusionList::writeTargets(const std::vector<PeptideIdentification> & pep_ids,
                                            const String & out_path,
                                            const IntList & charges)
  {
    WindowList result;

    Size charge_invalid_count(0);

    DoubleReal min_to_s_factor = (param_.getValue("RT:unit") == "seconds") ? 1.0 : (1.0 / 60.0);
    bool relative_rt = (param_.getValue("RT:use_relative") == "true");
    DoubleReal rel_rt_window_size = param_.getValue("RT:window_relative");
    DoubleReal abs_rt_window_size = param_.getValue("RT:window_absolute");

    std::vector<PeptideIdentification>::const_iterator pep_id_iter = pep_ids.begin();
    for (; pep_id_iter != pep_ids.end(); ++pep_id_iter)
    {
      if (pep_id_iter->getHits().size() > 1)
      {
        Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, pep_id_iter->getHits().size());
      }
      if (!pep_id_iter->metaValueExists("RT"))
      {
        Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Peptide identification contains no RT information.");
      }
      DoubleReal rt = pep_id_iter->getMetaValue("RT");

      DoubleReal rt_start = std::max(0.0, relative_rt ? (rt - rt * rel_rt_window_size) : rt - abs_rt_window_size);
      DoubleReal rt_stop =                relative_rt ? (rt + rt * rel_rt_window_size) : rt + abs_rt_window_size;

      rt_start *= min_to_s_factor;
      rt_stop *= min_to_s_factor;

      std::vector<PeptideHit>::const_iterator pep_hit_iter = pep_id_iter->getHits().begin();
      for (; pep_hit_iter != pep_id_iter->getHits().end(); ++pep_hit_iter)
      {
        Int charge = pep_hit_iter->getCharge();
        if (charge == 0)
        {
          ++charge_invalid_count;
          //fix charge
          charge = 2;
        }

        bool charge_found = false;
        for (Size c = 0; c < charges.size(); ++c)
        {
          DoubleReal mz = pep_hit_iter->getSequence().getMonoWeight(Residue::Full, charges[c]) / (DoubleReal)charges[c];
          result.push_back(IEWindow(rt_start, rt_stop, mz));
          if (charges[c] == charge)
          {
            charge_found = true;
          }
        }
        if (!charge_found) // if not already done, consider annotated charge of peptide (unless its 0)
        {
          DoubleReal mz = pep_hit_iter->getSequence().getMonoWeight(Residue::Full, charge) / (DoubleReal)charge;
          result.push_back(IEWindow(rt_start, rt_stop, mz));
        }
      }
    }

    if (charge_invalid_count > 0)
      LOG_WARN << "Warning: " << charge_invalid_count << " peptides with charge=0 were found, and assumed to have charge=2.\n";

    mergeOverlappingWindows_(result);
    writeToFile_(out_path, result);
  }

  void InclusionExclusionList::writeToFile_(const String & out_path,
                                            const WindowList & windows) const
  {

    std::ofstream outs(out_path.c_str());
    outs.precision(8);
    if (!outs)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Cannot open output file '" + out_path + "'.");
    }

    for (Size i = 0; i < windows.size(); ++i)
    {
      outs << windows[i].MZ_ << "\t" << windows[i].RTmin_ << "\t" << windows[i].RTmax_ << "\n";
    }


    outs.close();

  }

} // namespace
