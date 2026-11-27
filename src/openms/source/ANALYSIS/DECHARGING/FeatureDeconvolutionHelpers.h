// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Fabian Aicheler $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Private helpers shared by FeatureDeconvolution and MetaboliteFeatureDeconvolution.

      Both classes infer additional edges, filter edges by intensity and check the solution in the
      same way. The algorithms live here once; the classes' member functions forward to them. Not installed.
    */
    namespace FeatureDeconvolutionHelpers
    {
      /**
        @brief 1-sided Compomer for a feature

        Holds information on an explicit (with H+) 1-sided Compomer of a feature.
      */
      struct CmpInfo
      {
        std::string s_comp; ///< formula as String
        Size idx_cp{}; ///< index into compomer vector
        UInt side_cp{}; ///< side of parent compomer (LEFT or RIGHT)

        // C'tor
        CmpInfo() :
          s_comp() {}

        // C'tor
        CmpInfo(std::string& s, Size idx, UInt side) :
          s_comp(s), idx_cp(idx), side_cp(side) {}

        // Comparator
        bool operator<(const CmpInfo& other) const
        {
          return s_comp < other.s_comp;
        }

        bool operator==(const CmpInfo& other) const
        {
          return s_comp == other.s_comp;
        }
      };

      /**
        @brief test if "simple" edges have alternative (more difficult explanation) supported by neighboring edges,
               e.g. (.) -> (H+) might be augmented to (Na+) -> (H+Na+)

        @param[in,out] edges The edges; inferred edges are appended
        @param[in] feature_adducts 1-sided compomers of each feature (CmpInfoType derives from CmpInfo)
        @param[in] is_neg Negative ionization mode?
        @param[in] class_name Name of the calling class (for exception messages)
        @param[in] function Name of the calling function, reported by the exceptions

        @exception Exception::InvalidValue if an inferred edge has wrong charges
      */
      template <typename CmpInfoType>
      void inferMoreEdges(ILPDCWrapper::PairsType& edges, std::map<Size, std::set<CmpInfoType> >& feature_adducts, const bool is_neg, const char* class_name, const char* function)
      {
        Adduct default_adduct;

        if (is_neg)
        {
          default_adduct = Adduct(-1, 1, -Constants::PROTON_MASS_U, "H-1", std::log(1.0),0);
        }
        else
        {
          default_adduct = Adduct(1, 1, Constants::PROTON_MASS_U, "H1", std::log(1.0), 0);
        }

        int left_charges, right_charges;

        Size edges_size = edges.size();
        for (Size i = 0; i < edges_size; ++i)
        {
          Size idx0 = edges[i].getElementIndex(0);
          Size idx1 = edges[i].getElementIndex(1);

          std::set<CmpInfoType> result;

          // look at the intersection of the two adjacent features
          // (we thus require the non-H adduct to be present on both sides of the edge,
          //  if one is deemed enough just change to union)
          std::set_intersection(feature_adducts[idx0].begin(), feature_adducts[idx0].end(),
                                feature_adducts[idx1].begin(), feature_adducts[idx1].end(),
                                std::inserter(result, result.begin()));
          typename std::set<CmpInfoType>::iterator result_it = result.begin();

          // add new edge with each adduct of the intersection
          while (result_it != result.end())
          {
            Compomer::CompomerSide to_add = edges[result_it->idx_cp].getCompomer().removeAdduct(default_adduct).getComponent()[result_it->side_cp];
            // we currently do not punish additional two-side adducts
            for (Compomer::CompomerSide::iterator it = to_add.begin(); it != to_add.end(); ++it)
            {
              it->second.setLogProb(0);
            }
            ChargePair cp(edges[i]); // make a copy
            Compomer new_cmp = cp.getCompomer().removeAdduct(default_adduct);

            new_cmp.add(to_add, Compomer::LEFT);
            new_cmp.add(to_add, Compomer::RIGHT);

            //We again need to consider inverted behavior (but cp.getCharge(x) gets negative charges as assigned before!
            if (is_neg)
            {
              left_charges =  -new_cmp.getPositiveCharges();
              right_charges =  -new_cmp.getNegativeCharges();
            }
            else
            {
              left_charges = new_cmp.getNegativeCharges();
              right_charges = new_cmp.getPositiveCharges();
            }


            // refill with default adducts (usually H+):
            if (((cp.getCharge(0) - left_charges) % default_adduct.getCharge() == 0) &&
                ((cp.getCharge(1) - right_charges) % default_adduct.getCharge() == 0)) // for singly charged default_adducts this should always be true
            {
              int hc_left, hc_right;
              hc_left = (cp.getCharge(0) - left_charges) / default_adduct.getCharge();
              hc_right = (cp.getCharge(1) - right_charges) / default_adduct.getCharge();

              // we have not stepped over the charge capacity of the features
              if (hc_left >= 0 && hc_right >= 0)
              {
                // fill up with defaults:
                if (hc_left > 0)
                  new_cmp.add(default_adduct * hc_left, Compomer::LEFT);
                if (hc_right > 0)
                  new_cmp.add(default_adduct * hc_right, Compomer::RIGHT);

                // charge constraints of feature still fulfilled? (taking ionization mode into account)
                int left_charge, right_charge;
                if (is_neg)
                {
                  left_charge =  -new_cmp.getPositiveCharges();
                  right_charge =  -new_cmp.getNegativeCharges();
                }
                else
                {
                  left_charge = new_cmp.getNegativeCharges();
                  right_charge = new_cmp.getPositiveCharges();
                }

                if ((left_charge == cp.getCharge(0)) &&
                    (right_charge == cp.getCharge(1)))
                {
                  cp.setCompomer(new_cmp);
                  cp.setEdgeScore(0.99); //TODO how to score this new edge?
                  edges.push_back(cp); // add edge
                }
                else
                {
                  throw Exception::InvalidValue(__FILE__, __LINE__, function, std::string(class_name) + "::inferMoreEdges_(): Inferred edges with wrong(switched?) charges! Left neg_charge, left feature charge, right pos_charge, right feature charge",StringUtils::toStr(new_cmp.getNegativeCharges())+","+StringUtils::toStr(cp.getCharge(0))+","+StringUtils::toStr(new_cmp.getPositiveCharges())+","+StringUtils::toStr(cp.getCharge(1)));
                }
              }

            }
            else // have nonzero modulo.SHOULD NOT HAPPEN FOR DEFAULT CHARGE 1/-1 !!
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, function, std::string(class_name) + "::inferMoreEdges_(): Modulo returns leftover charge!",StringUtils::toStr(new_cmp.getNegativeCharges()));
            }

            ++result_it;
          }
        } // edge for

        OPENMS_LOG_INFO << "Inferring edges raised edge count from " << edges_size << " to " << edges.size() << "\n";
      }

      /// A function mostly for debugging
      inline void printEdgesOfConnectedFeatures(Size idx_1, Size idx_2, const ILPDCWrapper::PairsType& feature_relation)
      {
        std::cout << " +++++ printEdgesOfConnectedFeatures_ +++++\n";
        for (Size i = 0; i < feature_relation.size(); ++i)
        {
          if (
            ((feature_relation[i].getElementIndex(0) == idx_1) && (feature_relation[i].getElementIndex(1) == idx_2))
             ||
            ((feature_relation[i].getElementIndex(0) == idx_2) && (feature_relation[i].getElementIndex(1) == idx_1))
            )
          {
            std::cout << feature_relation[i].getCompomer() << " Edge: "  << i << " score: " <<  feature_relation[i].getEdgeScore() << "\n";
          }
        }
        std::cout << " ----- printEdgesOfConnectedFeatures_ -----\n";
        return;
      }

      /**
        @brief returns true if the intensity filter was passed or switched off (@p enable_intensity_filter is false)

        Filter for adding an edge only when the two features connected by it, fulfill the
        intensity criterion.
      **/
      inline bool intensityFilterPassed(const bool enable_intensity_filter, const Int q1, const Int q2, const Compomer& cmp, const Feature& f1, const Feature& f2)
      {
        if (!enable_intensity_filter)
          return true;

        if (q1 == q2)
        {
          Compomer cl; cl.add(cmp.getComponent()[0], Compomer::LEFT);
          Compomer cr; cr.add(cmp.getComponent()[1], Compomer::LEFT);
          if (((cl.getLogP() <= cr.getLogP()) && (f1.getIntensity() <= f2.getIntensity()))
             ||
              ((cl.getLogP() >= cr.getLogP()) && (f1.getIntensity() >= f2.getIntensity()))
              )
          {
            return true;
          }
          else
          {
            // forbid this edge?!
            std::cout << "intensity constraint: edge with intensity " << f1.getIntensity() << "(" << cmp.getAdductsAsString(Compomer::LEFT) << ") and " << f2.getIntensity() << "(" << cmp.getAdductsAsString(Compomer::RIGHT) << ") deleted\n";
            return false;
          }
        }
        return true;
      }

      /// test for obviously wrong parameter settings (too high lower charge bound) and warn
      inline void checkSolution(const ConsensusMap& cons_map)
      {
        Size ladders_total(0);
        Size ladders_with_odd(0);

        // checking number of charge ladders which have all gapped shapes, hinting at wrong lower-bound bound (should be lower)
        for (const ConsensusFeature& cfeature : cons_map)
        {
          if (cfeature.size() == 1)
            continue;

          ++ladders_total;
          IntList charges = cfeature.getMetaValue("distinct_charges").toIntList();

          for (Size i = 0; i < charges.size(); ++i)
          {
            if (charges[i] % 2 == 1)
            {
              ++ladders_with_odd;
              break;
            }
          }
        }

        // if more than 5% of charge ladder have only gapped, report
        if (ladders_with_odd < ladders_total * 0.95)
        {
          OPENMS_LOG_WARN << ".\n..\nWarning: a significant portion of your decharged molecules have gapped, even-numbered charge ladders (" << ladders_total - ladders_with_odd << " of " << ladders_total << ")";
          OPENMS_LOG_WARN <<"This might indicate a too low charge interval being tested.\n..\n.\n";
        }
      }
    } // namespace FeatureDeconvolutionHelpers
  } // namespace Internal
} // namespace OpenMS
