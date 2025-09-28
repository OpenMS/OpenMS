// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <fstream>
#include <map>

namespace OpenMS
{

  ILPDCWrapper::ILPDCWrapper() = default;

  ILPDCWrapper::~ILPDCWrapper()  = default;

  double ILPDCWrapper::compute(const FeatureMap& fm, PairsType& pairs, Size verbose_level) const
  {
    if (fm.empty())
    {
      OPENMS_LOG_INFO << "ILPDC wrapper received empty feature list. Nothing to compute! Exiting..." << std::endl;
      return -1;
    }

    PairsType pairs_clique_ordered;
    pairs_clique_ordered.reserve(pairs.size());
    typedef std::vector<std::pair<Size, Size> > BinType;
    BinType bins;
    // check number of components for complete putative edge graph (usually not all will be set to 'active' during ILP):
    {
      //
      // find groups of edges
      //
      Size group_count(0);
      std::map<Size, Size> f2g; // feature id to connected group
      std::map<Size, std::set<Size> > g2pairs; // group id to all pairs involved
      std::map<Size, std::set<Size> > g2f; // group id to all features involved

      for (Size i = 0; i < pairs.size(); ++i)
      {
        Size f1 = pairs[i].getElementIndex(0);
        Size f2 = pairs[i].getElementIndex(1);
        if ((f2g.find(f1) != f2g.end()) && (f2g.find(f2) != f2g.end())) // edge connects two distinct groups
        {
          Size group1 = f2g[f1];
          Size group2 = f2g[f2];
          if (group2 != group1)
          {
            // point group2 to group1
            g2pairs[group1].insert(g2pairs[group2].begin(), g2pairs[group2].end());
            g2pairs.erase(group2);
            for (std::set<Size>::const_iterator its = g2f[group2].begin(); its != g2f[group2].end(); ++its)
            {
              g2f[group1].insert(*its);
              f2g[*its] = group1; // reassign features of group2 to group1
            }
            g2f.erase(group2);
          }
        }
        else if ((f2g.find(f1) != f2g.end()) && !(f2g.find(f2) != f2g.end())) // only f1 is part of a group
        {
          Size group1 = f2g[f1];
          f2g[f2] = group1;
          g2f[group1].insert(f2);
        }
        else if (!(f2g.find(f1) != f2g.end()) && f2g.find(f2) != f2g.end()) // only f2 is part of a group
        {
          Size group2 = f2g[f2];
          f2g[f1] = group2;
          g2f[group2].insert(f1);
        }
        else // neither feature has a group: make a new one
        {
          Size group = ++group_count;
          f2g[f1] = group;
          f2g[f2] = group;
          g2f[group].insert(f1);
          g2f[group].insert(f2);
        }
        // append current edge to common group
        g2pairs[f2g[f1]].insert(i);
      }

      if (g2pairs.size() != g2f.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Clique construction failed! Unequal number of groups produced!", String(g2pairs.size()) + "!=" + String(g2f.size()));
      }

      std::map<Size, Size> hist_component_sum;
      // now walk though groups and see the size:
      for (std::map<Size, std::set<Size> >::const_iterator it = g2f.begin(); it != g2f.end(); ++it)
      {
        ++hist_component_sum[it->second.size()]; // e.g. component 2 has size 4; thus increase count for size 4
      }
      if (verbose_level > 1)
      {
        OPENMS_LOG_INFO << "Components:\n";
        OPENMS_LOG_INFO << "  Size 1 occurs ?x\n";
        for (std::map<Size, Size>::const_iterator it = hist_component_sum.begin(); it != hist_component_sum.end(); ++it)
        {
          OPENMS_LOG_INFO << "  Size " << it->first << " occurs " << it->second << "x\n";
        }
      }

      /* partition the cliques into bins, one given to the ILP at a time */
      UInt pairs_per_bin = 1000;
      UInt big_clique_bin_threshold = 200;

      Size start(0);
      Size count(0);
      for (std::map<Size, std::set<Size> >::const_iterator it = g2pairs.begin(); it != g2pairs.end(); ++it)
      {
        Size clique_size = it->second.size();
        if (count > pairs_per_bin || clique_size > big_clique_bin_threshold)
        {
          if (count > 0) // either bin is full or we have to close it due to big clique
          {
            if (verbose_level > 2)
              OPENMS_LOG_INFO << "Overstepping border of " << pairs_per_bin << " by " << SignedSize(count - pairs_per_bin) << " elements!\n";
            bins.push_back(std::make_pair(start, pairs_clique_ordered.size()));
            start = pairs_clique_ordered.size();
            count = 0;
          }
          if (clique_size > big_clique_bin_threshold) // extra bin for this big clique
          {
            for (std::set<Size>::const_iterator i_p = it->second.begin(); i_p != it->second.end(); ++i_p)
            {
              pairs_clique_ordered.push_back(pairs[*i_p]);
            }
            if (verbose_level > 2)
              OPENMS_LOG_INFO << "Extra bin for big clique (" << clique_size << ") prepended to schedule\n";
            bins.insert(bins.begin(), std::make_pair(start, pairs_clique_ordered.size()));
            start = pairs_clique_ordered.size();
            continue; // next clique (this one is already processed)
          }
        }
        count += clique_size;
        for (std::set<Size>::const_iterator i_p = it->second.begin(); i_p != it->second.end(); ++i_p)
        {
          pairs_clique_ordered.push_back(pairs[*i_p]);
        }
      }
      if (count > 0)
        bins.push_back(std::make_pair(start, pairs_clique_ordered.size()));
    }

    if (pairs_clique_ordered.size() != pairs.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pairs_clique_ordered.size() - pairs.size());
    }
    /* swap pairs, such that edges are order by cliques (so we can make clean cuts) */
    pairs.swap(pairs_clique_ordered);

    //PairsType pt2 = pairs;
    StopWatch time1;
    time1.start();

    // split problem into slices and have each one solved by the ILPS
    double score = 0;
// OMP currently causes spurious segfaults in Release mode; OMP fix applied, however: disable if problem persists
//#ifdef _OPENMP
//#pragma omp parallel for schedule(dynamic, 1), reduction(+: score)
//#endif
    for (SignedSize i = 0; i < static_cast<SignedSize>(bins.size()); ++i)
    {
      score = computeSlice_(fm, pairs, bins[i].first, bins[i].second, verbose_level);
    }
    time1.stop();
    OPENMS_LOG_INFO << " Branch and cut took " << time1.getClockTime() << " seconds, "
             << " with objective value: " << score << "."
             << std::endl;

    return score;
  }

  void ILPDCWrapper::updateFeatureVariant_(FeatureType_& f_set, const String& rota_l, const Size& v) const
  {
    f_set[rota_l].insert(v);
  }

  double ILPDCWrapper::computeSlice_(const FeatureMap& fm,
                                     PairsType& pairs,
                                     const PairsIndex margin_left,
                                     const PairsIndex margin_right,
                                     const Size /* verbose_level */) const
  {
    // feature --> variants set  (with scores)
    typedef std::map<Size, FeatureType_> r_type;
    r_type features;


    LPWrapper build;
    build.setObjectiveSense(LPWrapper::MAX); // maximize

    // add ALL edges first. Their result is what is interesting to us later
    for (PairsIndex i = margin_left; i < margin_right; ++i)
    {
      // log scores are good for addition in ILP - but they are < 0, thus not suitable for maximizing
      // ... so we just add normal probabilities...
      double score = exp(getLogScore_(pairs[i], fm));
      pairs[i].setEdgeScore(score * pairs[i].getEdgeScore()); // multiply with preset score

      // create the column representing the edge
      Int index = build.addColumn();
      build.setColumnBounds(index, 0, 1, LPWrapper::DOUBLE_BOUNDED);
      build.setColumnType(index, LPWrapper::INTEGER); // integer variable
      build.setObjective(index, pairs[i].getEdgeScore());

      // create feature variants set
      String rota_l = String(pairs[i].getElementIndex(0)) + pairs[i].getCompomer().getAdductsAsString(0) + "_" + pairs[i].getCharge(0);
      updateFeatureVariant_(features[pairs[i].getElementIndex(0)], rota_l, index);
      String rota_r = String(pairs[i].getElementIndex(1)) + pairs[i].getCompomer().getAdductsAsString(1) + "_" + pairs[i].getCharge(1);
      updateFeatureVariant_(features[pairs[i].getElementIndex(1)], rota_r, index);
    }

    // ADD Features (multiple variants of one feature are constrained to size=1)
    Size count(0); // each entry is a feature idx --->    Map["AdductCgf"]->adjacentEdges
    for (r_type::iterator it = features.begin(); it != features.end(); ++it)
    {
      ++count;
      std::vector<Int> columns;
      std::vector<double> elements;
      for (FeatureType_::const_iterator iti = it->second.begin(); iti != it->second.end(); ++iti)
      {
        Int index = build.addColumn();
        build.setColumnBounds(index, 0, 1, LPWrapper::DOUBLE_BOUNDED);
        build.setColumnType(index, LPWrapper::INTEGER); // integer variable
        build.setObjective(index, 0); // obj value of feature must be a constant, as it must be neutral
        columns.push_back(index);
        elements.push_back(1.0);

        /* allow connected edges only if this variant of the feature is chosen */
        /* get adjacent edges */
        std::vector<Int> columns_e;
        std::vector<double> elements_e;
        for (std::set<Size>::const_iterator it_e = iti->second.begin(); it_e != iti->second.end(); ++it_e)
        {
          columns_e.push_back((Int) * it_e);
          elements_e.push_back(-1.0);
        }
        columns_e.push_back((Int) index);
        elements_e.push_back(iti->second.size()); // factor of variant is number of adjacent edges
        String se = String("cv") + index;
        build.addRow(columns_e, elements_e, se, 0, 10000, LPWrapper::LOWER_BOUND_ONLY);
      }
      String s = String("c") + count;
      // only allow exactly one charge variant
      build.addRow(columns, elements, s, 1, 1, LPWrapper::FIXED);
    }

    LPWrapper::SolverParam param;
    param.enable_mir_cuts = true;
    param.enable_cov_cuts = true;
    param.enable_feas_pump_heuristic = true;
    param.enable_binarization = false;
    param.enable_clq_cuts = true;
    param.enable_gmi_cuts = true;
    param.enable_presolve = true;

    build.solve(param);

    for (UInt iColumn = 0; iColumn < margin_right - margin_left; ++iColumn)
    {
      double value = build.getColumnValue(iColumn);
      if (fabs(value) > 0.5)
      {
        pairs[margin_left + iColumn].setActive(true);
      }
      else
      {
        // DEBUG
        //std::cerr << " edge " << iColumn << " with " << value << "\n";
      }
    }

    return build.getObjectiveValue();

  }

  // old version, slower, as ILP has different layout (i.e, the same as described in paper)

  double ILPDCWrapper::computeSliceOld_(const FeatureMap& fm,
                                        PairsType& pairs,
                                        const PairsIndex margin_left,
                                        const PairsIndex margin_right,
                                        const Size verbose_level) const
  {
    LPWrapper build;
    build.setObjectiveSense(LPWrapper::MAX); // maximize

    //------------------------------------objective function-----------------------------------------------
    // find maximal objective value
    double score_min = 10e10f, score_max = -10e10f;

    // fill in objective values
    std::ostringstream namebuf;

    for (PairsIndex i = margin_left; i < margin_right; ++i)
    {

      // log scores are good for addition in ILP - but they are < 0, thus not suitable for maximizing
      // ... so we just add normal probabilities...
      double score = exp(getLogScore_(pairs[i], fm));
      pairs[i].setEdgeScore(score * pairs[i].getEdgeScore()); // multiply with preset score
      namebuf.str("");
      namebuf << "x#" << i;
      // create the new variable object
      Int index = build.addColumn();
      build.setColumnBounds(index, 0, 1, LPWrapper::DOUBLE_BOUNDED);
      build.setColumnType(index, LPWrapper::INTEGER); // integer variable
      build.setObjective(index, pairs[i].getEdgeScore());
      if (score_min > score)
        score_min = score;
      if (score_max < score)
        score_max = score;

      // DEBUG:
      //std::cerr << "MIP: edge#"<< i << " score: " << pairs[i].getEdgeScore() << " adduct:" << pairs[i].getCompomer().getAdductsAsString() << "\n";
    }
    if (verbose_level > 2)
      OPENMS_LOG_INFO << "score_min: " << score_min << " score_max: " << score_max << "\n";

    //------------------------------------adding constraints--------------------------------------------------

    bool is_conflicting;
    std::vector<int> conflict_idx(4);

    for (PairsIndex i = margin_left; i < margin_right; ++i)
    {
      const Compomer& ci = pairs[i].getCompomer();

      // TODO: only go until next clique...
      for (PairsIndex j = i + 1; j < margin_right; ++j)
      {
        const Compomer& cj = pairs[j].getCompomer();

        is_conflicting = false;
        // add pairwise constraints (one for each two conflicting ChargePairs)
        // if features are identical they must have identical charges (because any single
        // feature can only have one unique charge)


        //outgoing edges (from one feature)
        if (pairs[i].getElementIndex(0) == pairs[j].getElementIndex(0))
        {
          if ((pairs[i].getCharge(0) != pairs[j].getCharge(0))  ||
              ci.isConflicting(cj, Compomer::LEFT, Compomer::LEFT))
          {
            is_conflicting = true;
            ++conflict_idx[0];
          }
        }

        //incoming edges (into one feature)
        if (pairs[i].getElementIndex(1) == pairs[j].getElementIndex(1))
        {
          if ((pairs[i].getCharge(1) != pairs[j].getCharge(1))  ||
              ci.isConflicting(cj, Compomer::RIGHT, Compomer::RIGHT))
          {
            is_conflicting = true;
            ++conflict_idx[1];
          }
        }

        //incoming/outgoing edge (from one feature)
        if (pairs[i].getElementIndex(1) == pairs[j].getElementIndex(0))
        {
          if ((pairs[i].getCharge(1) != pairs[j].getCharge(0))  ||
              ci.isConflicting(cj, Compomer::RIGHT, Compomer::LEFT))
          {
            is_conflicting = true;
            ++conflict_idx[2];
          }
        }

        //incoming/outgoing edge (from one feature) -- this should only happen to additionally inferred edges
        if (pairs[i].getElementIndex(0) == pairs[j].getElementIndex(1))
        {
          if ((pairs[i].getCharge(0) != pairs[j].getCharge(1))  ||
              ci.isConflicting(cj, Compomer::LEFT, Compomer::RIGHT))
          {
            is_conflicting = true;
            ++conflict_idx[3];
          }
        }


        if (is_conflicting)
        {
          String s = String("C") + i + "." + j;

          // Now build rows: two variables, with indices 'columns', factors '1', and 0-1 bounds.
          std::vector<double> element(2, 1.0);
          std::vector<int> columns;
          columns.push_back(int(i - margin_left));
          columns.push_back(int(j - margin_left));
          build.addRow(columns, element, s, 0., 1., LPWrapper::DOUBLE_BOUNDED);
        }
      }
    }

    if (verbose_level > 2)
      OPENMS_LOG_INFO << "node count: " << fm.size() << "\n";
    if (verbose_level > 2)
      OPENMS_LOG_INFO << "edge count: " << pairs.size() << "\n";
    if (verbose_level > 2)
      OPENMS_LOG_INFO << "constraint count: " << (conflict_idx[0] + conflict_idx[1] + conflict_idx[2] + conflict_idx[3]) << " = " << conflict_idx[0] << " + " << conflict_idx[1] << " + " << conflict_idx[2] << " + " << conflict_idx[3] << "(0 or inferred)" << std::endl;

    //---------------------------------------------------------------------------------------------------------
    //----------------------------------------Solving and querying result--------------------------------------
    //---------------------------------------------------------------------------------------------------------

    if (verbose_level > 0)
      OPENMS_LOG_INFO << "Starting to solve..." << std::endl;
    LPWrapper::SolverParam param;
    param.enable_mir_cuts = true;
    param.enable_cov_cuts = true;
    param.enable_feas_pump_heuristic = true;
    param.enable_binarization = false;
    param.enable_clq_cuts = true;
    param.enable_gmi_cuts = true;
    param.enable_presolve = true;
    StopWatch time1;
    time1.start();
    build.solve(param);
    time1.stop();
    if (verbose_level > 0)
      OPENMS_LOG_INFO << " Branch and cut took " << time1.getClockTime() << " seconds, "
               << " with objective value: " << build.getObjectiveValue() << "."
               << " Status: " << (!build.getStatus() ? " Finished" : " Not finished")
               << std::endl;

    // variable values
    UInt active_edges = 0;
    std::map<String, Size> count_cmp;

    for (Int iColumn = 0; iColumn < build.getNumberOfColumns(); ++iColumn)
    {
      double value = build.getColumnValue(iColumn);
      if (fabs(value) > 0.5)
      {
        ++active_edges;
        pairs[margin_left + iColumn].setActive(true);
        // for statistical purposes: collect compomer distribution
        String cmp = pairs[margin_left + iColumn].getCompomer().getAdductsAsString();
        ++count_cmp[cmp];
      }
      else
      {
        // DEBUG
        //std::cerr << " edge " << iColumn << " with " << value << "\n";
      }
    }
    if (verbose_level > 2)
      OPENMS_LOG_INFO << "Active edges: " << active_edges << " of overall " << pairs.size() << std::endl;

    for (std::map<String, Size>::const_iterator it = count_cmp.begin(); it != count_cmp.end(); ++it)
    {
      //std::cout << "Cmp " << it->first << " x " << it->second << "\n";
    }

    double opt_value = build.getObjectiveValue();

    //objective function value of optimal(?) solution
    return opt_value;
  } // !compute_slice

  double ILPDCWrapper::getLogScore_(const PairsType::value_type& pair, const FeatureMap& fm) const
  {
    double score;
    String e;
    if (getenv("M") != nullptr)
      e = String(getenv("M"));
    if (e.empty())
    {
      //std::cout << "1";
      score = pair.getCompomer().getLogP();
      /*double charge_enhance = 0;

      if (pairs[i].getCharge(0) == fm[pairs[i].getElementIndex(0)].getCharge())
          charge_enhance += log(0.9); else charge_enhance += log(0.1);

      if (pairs[i].getCharge(1) == fm[pairs[i].getElementIndex(1)].getCharge())
          charge_enhance += log(0.9); else charge_enhance += log(0.1);

      score += charge_enhance;
      */

    }
    else
    {
      //std::cout << "2";
      double rt_diff =  fabs(fm[pair.getElementIndex(0)].getRT() - fm[pair.getElementIndex(1)].getRT());
      // enhance correct charge
      double charge_enhance = ((pair.getCharge(0) == fm[pair.getElementIndex(0)].getCharge())
                              &&
                               (pair.getCharge(1) == fm[pair.getElementIndex(1)].getCharge()))
                              ? 100 : 1;
      score = charge_enhance * (1 / (pair.getMassDiff() + 1) + 1 / (rt_diff + 1));
    }

    //std::cout << "logscore: " << score << "\n";

    return score;
  }

}
