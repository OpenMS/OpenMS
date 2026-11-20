/*******************************************************************************
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef FDRWALK_H_
#define FDRWALK_H_

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <stdexcept>

#include "Globals.h"

namespace OpenMS { namespace Internal { namespace Percolator {


/*
 * The walk over a list of targets and decoys, best first, that estimates the
 * FDR at the end of each group of tied scores: the one loop behind
 * Scores::calcQvals, PosteriorEstimator::getQValues and the count-only grid
 * evaluation of Scores::countTargetsBySelection. Their results have to agree
 * to the last bit, so the expressions below must not be reordered.
 *
 * This is a reimplementation of
 *   Crux/src/app/AssignConfidenceApplication.cpp::compute_decoy_qvalues_mixmax
 * which itself was a reimplementation of Uri Keich's code written in R. With
 * pi0 == 1.0 it is the "traditional" estimate, (decoys + 1) / targets.
 *
 * n               number of entries
 * isTarget(i)     whether entry i is a target
 * tiedWithNext(i) whether entries i and i + 1 have the same score (i + 1 < n)
 * pi0             prior probability of an incorrect target; below 1 the
 *                 mix-max correction is applied
 * decoyFactor     expected number of incorrect targets per decoy
 * skipDecoysPlusOne  leaves out the + 1 of the decoy count
 * group(begin, end, numGroupTargets, numTargets, fdr)  is called for each
 *                 group of tied scores [begin, end) with the number of its
 *                 targets, the number of targets up to its end and the
 *                 estimate at its end, which is not capped at 1.
 *
 * The mix-max correction of a group with decoys uses N_{w<=z} and N_{z<=z},
 * the numbers of targets and decoys that score at most as much as the
 * (N_{z>=w} - 1)-th best decoy. They follow from the totals and the counts
 * above the group. With skipDecoysPlusOne that decoy does not exist for the
 * best group with decoys if it holds a single decoy: std::out_of_range.
 */
template <class IsTarget, class TiedWithNext, class Group>
void walkFdrEstimates(std::size_t n,
                      IsTarget isTarget,
                      TiedWithNext tiedWithNext,
                      double pi0,
                      double decoyFactor,
                      bool skipDecoysPlusOne,
                      Group group) {
  int numTargets = 0, numDecoys = 0;
  if (pi0 < 1.0) {
    for (std::size_t ix = 0; ix < n; ++ix) {
      (isTarget(ix) ? numTargets : numDecoys) += 1;
    }
  }
  // the counts at the last tie group with decoys
  int lastCntW = 0, lastCntZ = 0;

  double estPx_lt_zj = 0.0;
  double E_f1_mod_run_tot = 0.0;
  double fdrValue = 0.0;

  int n_z_ge_w = 1, n_w_ge_w = 0;  // N_{z>=w} and N_{w>=w}
  if (skipDecoysPlusOne) n_z_ge_w = 0;
  int numDecoysAbove = 0;  // decoys ranked at or above the current entry

  int decoyQueue = 0, targetQueue = 0;  // handles ties
  std::size_t groupBegin = 0;
  for (std::size_t ix = 0; ix < n; ++ix) {
    if (isTarget(ix)) {
      ++n_w_ge_w;  // target PSM
      ++targetQueue;
    } else {
      ++n_z_ge_w;  // decoy PSM
      ++numDecoysAbove;
      ++decoyQueue;
    }

    // handles ties
    if (ix + 1 == n || !tiedWithNext(ix)) {
      if (pi0 < 1.0 && decoyQueue > 0) {
        const int cntW = numTargets - (n_w_ge_w - targetQueue);
        const int cntZ = numDecoys - (numDecoysAbove - decoyQueue);
        // rank of the decoy whose counts are used, 1 for the best one
        const int rank = n_z_ge_w - 1;
        int cnt_w = cntW, cnt_z = cntZ;
        if (rank <= numDecoysAbove - decoyQueue) {
          // with skipDecoysPlusOne the decoy ranked just above a group with
          // a single decoy, i.e. in the last group with decoys
          if (rank < 1) {
            throw std::out_of_range(
                "walkFdrEstimates: no mix-max counts above the best decoy");
          }
          cnt_w = lastCntW;
          cnt_z = lastCntZ;
        }
        lastCntW = cntW;
        lastCntZ = cntZ;
        estPx_lt_zj = (double)(cnt_w - pi0 * cnt_z) / ((1.0 - pi0) * cnt_z);
        estPx_lt_zj = estPx_lt_zj > 1 ? 1 : estPx_lt_zj;
        estPx_lt_zj = estPx_lt_zj < 0 ? 0 : estPx_lt_zj;
        E_f1_mod_run_tot += decoyQueue * estPx_lt_zj * (1.0 - pi0);
        if (VERB > 4) {
          std::cerr << "Mix-max num negatives correction: "
                    << (1.0 - pi0) * n_z_ge_w << " vs. " << E_f1_mod_run_tot
                    << std::endl;
        }
      }

      fdrValue = (n_z_ge_w * pi0 + E_f1_mod_run_tot) /
                 (double)((std::max)(1, n_w_ge_w)) * decoyFactor;
      group(groupBegin, ix + 1, targetQueue, n_w_ge_w, fdrValue);
      groupBegin = ix + 1;
      decoyQueue = 0;
      targetQueue = 0;
    }
  }
}

}}}  // namespace OpenMS::Internal::Percolator

#endif /* FDRWALK_H_ */
