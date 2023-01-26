// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  class OPENMS_DLLAPI FeatureOverlapFilter
  {
    public:   
    /*
        @brief Filter overlapping features using a spatial datastructure (quadtree). 
               Retains only the best feature in each cluster of overlapping features.

        @param FeatureComparator must implement the concept of a less comparator.
               If several features overlap, the feature that evaluates as "smallest" is considered the best (according to the passed comparator) and is kept.
               The other overlapping features are removed and FeatureOverlapCallback evaluated on them.
               Default: overall feature quality.

        @param FeatureOverlapCallback(best_in_cluster, f) is called if a feature f overlaps with a feature best_in_cluster.
               FeatureOverlapCallback provides a customization point to e.g.:
              - transfer information from the soon-to-be-removed feature f over to the best_in_cluster feature
              - gather overlap statistics
              - help in debugging
              - etc.
              in form of a callable.
              If the FeatureOverlapCallback returns false, the overlapping feature will be treated as not overlapping with best_in_cluster (and not removed).
              Default: function that just returns true.

        @ingroup Datareduction
    */
    static void filter(FeatureMap& fmap, 
      std::function<bool(const Feature&, const Feature&)> FeatureComparator = [](const Feature& left, const Feature& right){ return left.getOverallQuality() > right.getOverallQuality(); },
      std::function<bool(Feature&, Feature&)> FeatureOverlapCallback = [](Feature&, Feature&){ return true; },
      bool check_overlap_at_trace_level = true);
  };

}


