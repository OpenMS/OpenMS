// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <vector>
#include <set>
#include <map> // for multimap<>
#include <algorithm> // for min
#include <numeric> // for make_pair

using std::map;
using std::vector;
using std::make_pair;
using std::min;
using std::set;
using std::iterator;

namespace OpenMS
{

  QTCluster::QTCluster()
  {
  }

  QTCluster::QTCluster(OpenMS::GridFeature* center_point, Size num_maps,
                       double max_distance, bool use_IDs, Int x_coord, 
                       Int y_coord) :
    center_point_(center_point),
    neighbors_(),
    tmp_neighbors_(nullptr),
    max_distance_(max_distance),
    num_maps_(num_maps),
    quality_(0.0),
    changed_(false),
    use_IDs_(use_IDs),
    valid_(true),
    collect_annotations_(false),
    finalized_(true),
    x_coord_(x_coord),
    y_coord_(y_coord),
    annotations_()
  {
    if (use_IDs)
    {
      annotations_ = center_point->getAnnotations();
    }
    if (use_IDs_ && center_point_->getAnnotations().empty())
    { 
      collect_annotations_ = true;
    }
  }

  GridFeature* QTCluster::getCenterPoint() 
  {
    return center_point_;
  }

  double QTCluster::getCenterRT() const
  {
    return center_point_->getRT();
  }

  double QTCluster::getCenterMZ() const
  {
    return center_point_->getMZ();
  }

  Int QTCluster::getXCoord() const
  {
    return x_coord_;
  }

  Int QTCluster::getYCoord() const
  {
    return y_coord_;
  }

  void QTCluster::setInvalid()
  {
    // this cluster is considered invalid, it will never be used again in the
    // algorithm. This means we can clean up a bit and save some memory.
    valid_ = false;
    neighbors_.clear();
    annotations_.clear();
  }

  Size QTCluster::size() const
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")
    return neighbors_.size() + 1; // + 1 for the center
  }

  bool QTCluster::operator<(QTCluster& cluster)
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")
    return this->getQuality() < cluster.getQuality();
  }

  void QTCluster::add(OpenMS::GridFeature* element, double distance)
  {
    OPENMS_PRECONDITION(!finalized_,
        "Cannot perform operation on cluster that is not initialized")
    // ensure we only add compatible peptide annotations
    OPENMS_PRECONDITION(distance <= max_distance_,
        "Distance cannot be larger than max_distance")
    // collect_annotations_ implies tmp_neighbors_ != NULL, 
    OPENMS_PRECONDITION(!collect_annotations_ || tmp_neighbors_ != NULL, 
        "Initialize the cluster first before adding elements")

    Size map_index = element->getMapIndex();

    // Ensure we only add compatible peptide annotations. If the cluster center
    // has an annotation, then each added neighbor should have the same
    // annotation. If the center element has no annotation we add all elements
    // and select the optimal annotation later, using optimizeAnnotations_ 
    if (use_IDs_)
    {
      bool one_empty = (center_point_->getAnnotations().empty() || element->getAnnotations().empty());
      if (!one_empty) // both are annotated
      {
        if (center_point_->getAnnotations() != element->getAnnotations()) 
        {
          // Both annotations are non-empty and are unequal, we don't add
          return;
        }
      }
    }

    // We have to store annotations in a temporary map first if we collect all
    // annotations
    if (collect_annotations_ && map_index != center_point_->getMapIndex())
    {
      (*tmp_neighbors_)[map_index].insert(make_pair(distance, element));
      changed_ = true;
    }

    // Store best (closest) element:
    // Only add the element if either no element is present for the map or if
    // the element is closer than the current element for that map
    if (map_index != center_point_->getMapIndex())
    {
      if (neighbors_.find(map_index) == neighbors_.end() ||
          distance < neighbors_[map_index].first)
      {
        neighbors_[map_index] = make_pair(distance, element);
        changed_ = true;
      }
    }
  }

  void QTCluster::getElements(OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>& elements)
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")

    elements.clear();
    elements[center_point_->getMapIndex()] = center_point_;

    if (neighbors_.empty())
    {
      return;
    }

    // since we are finalized, we do not need to care about the annotation
    for (NeighborMap::const_iterator it = neighbors_.begin(); it != neighbors_.end(); ++it)
    {
      elements[it->first] = it->second.second;
    }
  }

  bool QTCluster::update(const OpenMSBoost::unordered_map<Size,
      OpenMS::GridFeature*>& removed)
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")

    // check if the cluster center was removed:
    for (OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>::const_iterator
        rm_it = removed.begin(); rm_it != removed.end(); ++rm_it)
    {
      // If center point was removed, then we are done and no more work is
      // required
      if (rm_it->second == center_point_)
      {
        this->setInvalid();
        return false;
      }
    }

    // update cluster contents, remove those elements we find in our cluster
    for (OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>::const_iterator
        rm_it = removed.begin(); rm_it != removed.end(); ++rm_it)
    {
      NeighborMap::iterator pos = neighbors_.find(rm_it->first);
      if (pos == neighbors_.end())
      {
        continue; // no points from this map
      }

      const NeighborPairType current_feature = pos->second;
      if (current_feature.second == rm_it->second) // remove this neighbor
      {
        changed_ = true;
        neighbors_.erase(pos);
      }
    }
    return changed_;
  }

  double QTCluster::getQuality()
  {
    // this should work for finalized and non finalized states
    if (changed_)
    {
      computeQuality_();
      changed_ = false;
    }
    return quality_;
  }

  void QTCluster::computeQuality_()
  {
    // ensure cluster is not finalized as we cannot call optimizeAnnotations_
    // in that case
    OPENMS_PRECONDITION(!finalized_,
        "Cannot perform operation on cluster that is finalized")

    Size num_other = num_maps_ - 1;
    double internal_distance = 0.0;
    if (!use_IDs_ || !center_point_->getAnnotations().empty() ||
        neighbors_.empty())
    {
      // if the cluster center is annotated with peptide IDs, the neighbors can
      // consist only of features with compatible IDs, so we don't need to
      // check again here
      Size counter = 0;
      for (NeighborMap::iterator it = neighbors_.begin(); it != neighbors_.end(); ++it)
      {
        internal_distance += it->second.first;
        counter++;
      }
      // add max. distance for missing cluster elements:
      internal_distance += (num_other - counter) * max_distance_;
    }
    else // find the annotation that gives the best quality
    {
      internal_distance = optimizeAnnotations_();
    }

    // normalize:
    internal_distance /= num_other;
    quality_ = (max_distance_ - internal_distance) / max_distance_;
  }

  OpenMSBoost::unordered_map<Size, std::vector<GridFeature*> > QTCluster::getAllNeighbors() 
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")

    OpenMSBoost::unordered_map<Size, std::vector<GridFeature*> > tmp;
    for (NeighborMap::iterator it = neighbors_.begin(); it != neighbors_.end(); ++it)
    {
      tmp[ it->first ].push_back(it->second.second);
    }
    return tmp;
  }

  const set<AASequence>& QTCluster::getAnnotations()
  {
    if (changed_ && use_IDs_ && center_point_->getAnnotations().empty() && !neighbors_.empty())
    {
      optimizeAnnotations_();
    }
    return annotations_;
  }

  double QTCluster::optimizeAnnotations_()
  {
    OPENMS_PRECONDITION(collect_annotations_,
        "QTCluster::optimizeAnnotations_ should only be called if we use collect_annotations_")
    OPENMS_PRECONDITION(tmp_neighbors_ != NULL,
        "QTCluster::optimizeAnnotations_ needs to have working tmp_neighbors_")
    OPENMS_PRECONDITION(!finalized_,
        "QTCluster::optimizeAnnotations_ cannot work on finalized cluster")

    // mapping: peptides -> best distance per input map
    map<set<AASequence>, vector<double> > seq_table;

    for (NeighborMapMulti::iterator n_it = tmp_neighbors_->begin();
        n_it != tmp_neighbors_->end(); ++n_it)
    {
      Size map_index = n_it->first;
      for (NeighborListType::iterator df_it = n_it->second.begin(); 
          df_it != n_it->second.end(); ++df_it)
      {
        double dist = df_it->first;
        const set<AASequence>& current = df_it->second->getAnnotations();
        map<set<AASequence>, vector<double> >::iterator pos =
          seq_table.find(current);
        if (pos == seq_table.end())
        {
          // new set of annotations, fill vector with max distance for all maps
          seq_table[current].resize(num_maps_, max_distance_);
          seq_table[current][map_index] = dist;
        }
        else 
        {
          // new dist. value for this input map
          pos->second[map_index] = min(dist, pos->second[map_index]);
        }
        if (current.empty()) // unannotated feature
        {
          // no need to check further (annotation-specific distances are worse
          // than this unspecific one):
          break;
        }
      }
    }

    // combine annotation-specific and unspecific distances 
    // (all unspecific ones are grouped as empty set<AASequence>):
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
      double total = std::accumulate(it->second.begin(), it->second.end(), 0.0);
      if (total < best_total)
      {
        best_pos = it;
        best_total = total;
      }
    }

    if (best_pos != seq_table.end())
    {
      annotations_ = best_pos->first;
    }

    // report elements that are compatible with the optimal annotation:
    neighbors_.clear();
    for (NeighborMapMulti::const_iterator n_it = tmp_neighbors_->begin();
         n_it != tmp_neighbors_->end(); ++n_it)
    {
      for (std::multimap<double, GridFeature*>::const_iterator df_it =
             n_it->second.begin(); df_it != n_it->second.end(); ++df_it)
      {
        const set<AASequence>& current = df_it->second->getAnnotations();
        if (current.empty() || (current == annotations_))
        {
          neighbors_[n_it->first] = make_pair(df_it->first, df_it->second);
          break; // found the best element for this input map
        }
      }
    }

    // one "max_dist." too many (from the input map of the cluster center):
    return best_total - max_distance_;
  }

  void QTCluster::finalizeCluster()
  {
    OPENMS_PRECONDITION(tmp_neighbors_ != NULL,
        "Try to finalize QTCluster that was not initialized")
    OPENMS_PRECONDITION(!finalized_,
        "Try to finalize QTCluster that was not initialized")

    // calls computeQuality_ if something changed since initialization. In
    // turn, computeQuality_ calls optimizeAnnotations_ if necessary which
    // ensures that the neighbors_ hash is populated correctly.
    getQuality();

    finalized_ = true;

    // delete memory again
    delete tmp_neighbors_;
    tmp_neighbors_ = nullptr;
  }

  void QTCluster::initializeCluster()
  {
    OPENMS_PRECONDITION(tmp_neighbors_ == NULL,
        "Try to initialize QTCluster that was not finalized")
    OPENMS_PRECONDITION(finalized_,
        "Try to initialize QTCluster that was not finalized")

    finalized_ = false;

    if (tmp_neighbors_ != nullptr)
    {
      // delete memory again (should never actually happen but lets make sure
      // we release the memory under all circumstances)
      delete tmp_neighbors_;
      tmp_neighbors_ = nullptr;
    }

    // create empty map 
    tmp_neighbors_ = new NeighborMapMulti();
  }

  QTCluster::~QTCluster()
  {
    if (tmp_neighbors_ != nullptr)
    {
      // delete memory again (should never actually happen but lets make sure
      // we release the memory under all circumstances)
      delete tmp_neighbors_;
      tmp_neighbors_ = nullptr;
    }
  }

} // namespace OpenMS
