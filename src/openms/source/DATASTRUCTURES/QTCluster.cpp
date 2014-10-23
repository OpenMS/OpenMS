// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <numeric>

#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace std;

namespace OpenMS
{

  QTCluster::QTCluster()
  {
  }

  QTCluster::QTCluster(GridFeature * center_point, Size num_maps,
                       double max_distance, bool use_IDs) :
    center_point_(center_point), neighbors_(), max_distance_(max_distance),
    num_maps_(num_maps), quality_(0.0), changed_(false), use_IDs_(use_IDs),
    annotations_(),
    valid_(true)
  {
    if (use_IDs)
      annotations_ = center_point->getAnnotations();
  }

  QTCluster::~QTCluster()
  {
  }

  double QTCluster::getCenterRT() const
  {
    return center_point_->getRT();
  }

  double QTCluster::getCenterMZ() const
  {
    return center_point_->getMZ();
  }

  Size QTCluster::size() const
  {
    return neighbors_.size() + 1;     // + 1 for the center
  }

  bool QTCluster::operator<(QTCluster & cluster)
  {
    return this->getQuality() < cluster.getQuality();
  }

  void QTCluster::add(GridFeature * element, double distance)
  {
    // maybe TODO: check here if distance is smaller than max. distance?
    // maybe TODO: check here if peptide annotations are compatible?
    // (currently, both is done in QTClusterFinder)
    Size map_index = element->getMapIndex();
    if (map_index != center_point_->getMapIndex())
    {
      neighbors_[map_index].insert(make_pair(distance, element));
      changed_ = true;
    }
  }

  void QTCluster::getElements(OpenMSBoost::unordered_map<Size, GridFeature *> & elements)
  {
    elements.clear();
    elements[center_point_->getMapIndex()] = center_point_;

    if (neighbors_.empty())
    {
      return;
    }

    // if necessary, compute the optimal annotation for the cluster first:
    if (changed_ && use_IDs_ && center_point_->getAnnotations().empty())
    {
      optimizeAnnotations_();
    }

    if (annotations_.empty() || !center_point_->getAnnotations().empty())
    {
      // no need to take annotations into account:
      for (NeighborMap::const_iterator it = neighbors_.begin();
           it != neighbors_.end(); ++it)
      {
        elements[it->first] = it->second.begin()->second;
      }
    }
    else     // find elements that are compatible with the optimal annotation:
    {
      for (NeighborMap::const_iterator n_it = neighbors_.begin();
           n_it != neighbors_.end(); ++n_it)
      {
        for (std::multimap<double, GridFeature *>::const_iterator df_it =
               n_it->second.begin(); df_it != n_it->second.end(); ++df_it)
        {
          const set<AASequence> & current = df_it->second->getAnnotations();
          if (current.empty() || (current == annotations_))
          {
            elements[n_it->first] = df_it->second;
            break;             // found the best element for this input map
          }
        }
      }
    }
  }

  bool QTCluster::update(const OpenMSBoost::unordered_map<Size, GridFeature *> & removed)
  {
    // check if the cluster center was removed:
    for (OpenMSBoost::unordered_map<Size, GridFeature *>::const_iterator rm_it = removed.begin();
         rm_it != removed.end(); ++rm_it)
    {
      if (rm_it->second == center_point_)
        return false;
    }
    // update the cluster contents:
    for (OpenMSBoost::unordered_map<Size, GridFeature *>::const_iterator rm_it = removed.begin();
         rm_it != removed.end(); ++rm_it)
    {
      NeighborMap::iterator pos = neighbors_.find(rm_it->first);
      if (pos == neighbors_.end())
        continue;                                  // no points from this map
      for (std::multimap<double, GridFeature *>::iterator feat_it =
             pos->second.begin(); feat_it != pos->second.end(); ++feat_it)
      {
        if (feat_it->second == rm_it->second)         // remove this neighbor
        {
          if (!use_IDs_ || (annotations_ == rm_it->second->getAnnotations()))
          {
            changed_ = true;
          }
          // else: removed neighbor doesn't have optimal annotation, so it can't
          // be a "true" cluster element => no need to recompute the quality
          pos->second.erase(feat_it);
          break;
        }
      }
      if (pos->second.empty())       // only neighbor from this map was just removed
      {
        neighbors_.erase(pos);
      }
    }
    return true;
  }

  double QTCluster::getQuality()
  {
    if (changed_)
    {
      computeQuality_();
      changed_ = false;
    }
    return quality_;
  }

  void QTCluster::computeQuality_()
  {
    Size num_other = num_maps_ - 1;
    double internal_distance = 0.0;
    if (!use_IDs_ || !center_point_->getAnnotations().empty() ||
        neighbors_.empty())
    {
      // if the cluster center is annotated with peptide IDs, the neighbors can
      // consist only of features with compatible IDs, so we don't need to check
      // again here
      Size counter = 0;
      for (NeighborMap::iterator it = neighbors_.begin();
           it != neighbors_.end(); ++it)
      {
        internal_distance += it->second.begin()->first;
        counter++;
      }
      // add max. distance for missing cluster elements:
      internal_distance += (num_other - counter) * max_distance_;
    }
    else     // find the annotation that gives the best quality
    {
      internal_distance = optimizeAnnotations_();
    }

    // normalize:
    internal_distance /= num_other;
    quality_ = (max_distance_ - internal_distance) / max_distance_;
  }

  const set<AASequence> & QTCluster::getAnnotations()
  {
    if (changed_ && use_IDs_ && center_point_->getAnnotations().empty() &&
        !neighbors_.empty())
      optimizeAnnotations_();
    return annotations_;
  }

  double QTCluster::optimizeAnnotations_()
  {
    // mapping: peptides -> best distance per input map
    map<set<AASequence>, vector<double> > seq_table;

    for (NeighborMap::iterator n_it = neighbors_.begin();
         n_it != neighbors_.end(); ++n_it)
    {
      Size map_index = n_it->first;
      for (std::multimap<double, GridFeature *>::iterator df_it =
             n_it->second.begin(); df_it != n_it->second.end(); ++df_it)
      {
        double dist = df_it->first;
        const set<AASequence> & current = df_it->second->getAnnotations();
        map<set<AASequence>, vector<double> >::iterator pos =
          seq_table.find(current);
        if (pos == seq_table.end())         // new set of annotations
        {
          seq_table[current].resize(num_maps_, max_distance_);
          seq_table[current][map_index] = dist;
        }
        else         // new dist. value for this input map
        {
          pos->second[map_index] = min(dist, pos->second[map_index]);
        }
        if (current.empty())         // unannotated feature
        {
          // no need to check further (annotation-specific distances are worse
          // than this unspecific one):
          break;
        }
      }
    }

    // combine annotation-specific and unspecific distances:
    map<set<AASequence>, vector<double> >::iterator unspecific =
      seq_table.find(set<AASequence>());
    if (unspecific != seq_table.end())
    {
      for (map<set<AASequence>, vector<double> >::iterator it =
             seq_table.begin(); it != seq_table.end(); ++it)
      {
        if (it == unspecific)
          continue;
        for (Size i = 0; i < num_maps_; ++i)
        {
          it->second[i] = min(it->second[i], unspecific->second[i]);
        }
      }
    }

    // compute distance totals -> best annotation set has smallest value:
    map<set<AASequence>, vector<double> >::iterator best_pos =
      seq_table.begin();
    double best_total = num_maps_ * max_distance_;
    for (map<set<AASequence>, vector<double> >::iterator it =
           seq_table.begin(); it != seq_table.end(); ++it)
    {
      double total = accumulate(it->second.begin(), it->second.end(), 0.0);
      if (total < best_total)
      {
        best_pos = it;
        best_total = total;
      }
    }
    if (best_pos != seq_table.end())
      annotations_ = best_pos->first;
    // one "max_dist." too many (from the input map of the cluster center):
    return best_total - max_distance_;
  }

}
