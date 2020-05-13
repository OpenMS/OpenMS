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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/QTCluster.h>

#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <numeric> // for make_pair

using std::map;
using std::vector;
using std::make_pair;
using std::min;
using std::set;
using std::iterator;

namespace OpenMS
{
  QTCluster::QTCluster(QTCluster::Data_* data, OpenMS::GridFeature* center_point, 
                       Size num_maps, double max_distance,
                       bool use_IDs, Int x_coord, Int y_coord, Size id) :
    quality_(0.0),
    valid_(true),
    changed_(false),
    use_IDs_(use_IDs),
    collect_annotations_(false),
    finalized_(true),
    data_(data)
  {
    data_->center_point_ = center_point;
    data_->id_ = id;
    data_->max_distance_ = max_distance;
    data_->num_maps_ = num_maps;
    data_->x_coord_ = x_coord;
    data_->y_coord_ = y_coord;

    if (use_IDs)
    {
      data_->annotations_ = center_point->getAnnotations();
    }
    if (use_IDs_ && data_->center_point_->getAnnotations().empty())
    { 
      collect_annotations_ = true;
    }
  }

  GridFeature* QTCluster::getCenterPoint() const
  {
    return data_->center_point_;
  }

  Size QTCluster::getId() const
  {
    return data_->id_;
  }

  double QTCluster::getCenterRT() const
  {
    return data_->center_point_->getRT();
  }

  double QTCluster::getCenterMZ() const
  {
    return data_->center_point_->getMZ();
  }

  Int QTCluster::getXCoord() const
  {
    return data_->x_coord_;
  }

  Int QTCluster::getYCoord() const
  {
    return data_->y_coord_;
  }

  void QTCluster::setInvalid()
  {
    // this cluster is considered invalid, it will never be used again in the
    // algorithm. This means we can clean up a bit and save some memory.
    valid_ = false;
    data_->annotations_.clear();
  }

  Size QTCluster::size() const
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")
    return data_->neighbors_.size() + 1; // + 1 for the center
  }

  bool QTCluster::operator<(QTCluster const& rhs)
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")

    return quality_ < rhs.quality_;
  }

  void QTCluster::add(OpenMS::GridFeature* element, double distance)
  {
    // get references on member that is used in this function
    NeighborMapMulti & tmp_neighbors_ = data_->tmp_neighbors_;

    OPENMS_PRECONDITION(!finalized_,
        "Cannot perform operation on cluster that is not initialized")
    // ensure we only add compatible peptide annotations
    OPENMS_PRECONDITION(distance <= data_->max_distance_,
        "Distance cannot be larger than max_distance")
    // collect_annotations_ implies tmp_neighbors_ != NULL, 
    OPENMS_PRECONDITION(!collect_annotations_ || !tmp_neighbors_.empty(), 
        "Initialize the cluster first before adding elements")

    Size map_index = element->getMapIndex();

    // get references on more members that are used in this function
    OpenMS::GridFeature & center_point = *(data_->center_point_);
    NeighborMap & neighbors_ = data_->neighbors_;
    
    // Ensure we only add compatible peptide annotations. If the cluster center
    // has an annotation, then each added neighbor should have the same
    // annotation. If the center element has no annotation we add all elements
    // and select the optimal annotation later, using optimizeAnnotations_ 
    if (use_IDs_)
    {
      bool one_empty = (center_point.getAnnotations().empty() || element->getAnnotations().empty());
      if (!one_empty) // both are annotated
      {
        if (center_point.getAnnotations() != element->getAnnotations()) 
        {
          // Both annotations are non-empty and are unequal, we don't add
          return;
        }
      }
    }

    // We have to store annotations in a temporary map first if we collect all
    // annotations
    if (collect_annotations_ && map_index != center_point.getMapIndex())
    {
      tmp_neighbors_[map_index].insert(make_pair(distance, element));
      changed_ = true;
    }

    // Store best (closest) element:
    // Only add the element if either no element is present for the map or if
    // the element is closer than the current element for that map
    if (map_index != center_point.getMapIndex())
    {
      if (neighbors_.find(map_index) == neighbors_.end() ||
          distance < neighbors_[map_index].first)
      {
        neighbors_[map_index] = make_pair(distance, element);
        changed_ = true;
      }
    }
  }

  QTCluster::ClusterElementsMap QTCluster::getElements() const
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")

    // copy neighbor map
    NeighborMap neighbors_copy = data_->neighbors_;

    // add center point to the copy
    neighbors_copy[data_->center_point_->getMapIndex()].second = data_->center_point_;

    return neighbors_copy;
  }

  bool QTCluster::update(ClusterElementsMap const& removed)
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")

    // check if the cluster center was removed:
    for (NeighborMap::const_iterator
        rm_it = removed.begin(); rm_it != removed.end(); ++rm_it)
    {
      // If center point was removed, then we are done and no more work is
      // required
      if (rm_it->second.second == data_->center_point_)
      {
        this->setInvalid();
        return false;
      }
    }

    // get references on member that is used in this function
    NeighborMap & neighbors_ = data_->neighbors_;

    // update cluster contents, remove those elements we find in our cluster
    for (NeighborMap::const_iterator
        rm_it = removed.begin(); rm_it != removed.end(); ++rm_it)
    {
      NeighborMap::iterator pos = neighbors_.find(rm_it->first);
      if (pos == neighbors_.end())
      {
        continue; // no points from this map
      }

      const NeighborPairType current_feature = pos->second;
      if (current_feature.second == rm_it->second.second) // remove this neighbor
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

  double QTCluster::getCurrentQuality() const
  {
    // ensure cluster is finalized
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is finalized")
  
    return quality_;
  }

  void QTCluster::computeQuality_()
  {
    // ensure cluster is not finalized as we cannot call optimizeAnnotations_
    // in that case
    OPENMS_PRECONDITION(!finalized_,
        "Cannot perform operation on cluster that is finalized")

    // get references on member that is used in this function
    NeighborMap & neighbors_ = data_->neighbors_;

    // get copy of member that is used in this function
    double max_distance_ = data_->max_distance_;

    Size num_other = data_->num_maps_ - 1;
    double internal_distance = 0.0;
    if (!use_IDs_ || !data_->center_point_->getAnnotations().empty() ||
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

  OpenMSBoost::unordered_map<Size, std::vector<GridFeature*> > QTCluster::getAllNeighbors() const
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")

    // get references on member that is used in this function
    NeighborMap & neighbors_ = data_->neighbors_;

    OpenMSBoost::unordered_map<Size, std::vector<GridFeature*> > tmp;
    for (NeighborMap::iterator it = neighbors_.begin(); it != neighbors_.end(); ++it)
    {
      tmp[ it->first ].push_back(it->second.second);
    }
    return tmp;
  }

  QTCluster::NeighborMap const& QTCluster::getAllNeighborsDirect() const
  {
    OPENMS_PRECONDITION(finalized_,
        "Cannot perform operation on cluster that is not finalized")
        
    return data_->neighbors_;
  }

  const set<AASequence>& QTCluster::getAnnotations()
  {
    if (changed_ && use_IDs_ && data_->center_point_->getAnnotations().empty() && !data_->neighbors_.empty())
    {
      optimizeAnnotations_();
    }
    return data_->annotations_;
  }

  double QTCluster::optimizeAnnotations_()
  {
    OPENMS_PRECONDITION(collect_annotations_,
        "QTCluster::optimizeAnnotations_ should only be called if we use collect_annotations_")
    OPENMS_PRECONDITION(!data_->tmp_neighbors_.empty(),
        "QTCluster::optimizeAnnotations_ needs to have working tmp_neighbors_")
    OPENMS_PRECONDITION(!finalized_,
        "QTCluster::optimizeAnnotations_ cannot work on finalized cluster")

    // mapping: peptides -> best distance per input map
    map<set<AASequence>, vector<double> > seq_table;

    makeSeq_table_(seq_table);

    // get copies of members that are used in this function
    Size num_maps_ = data_->num_maps_;
    double max_distance_ = data_->max_distance_;
    
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
      data_->annotations_ = best_pos->first;
    }

    recompute_neighbors_();

    // one "max_dist." too many (from the input map of the cluster center):
    return best_total - max_distance_;
  }

  void QTCluster::recompute_neighbors_()
  {
    // get references on members that are used in this function
    NeighborMap & neighbors_ = data_->neighbors_;
    NeighborMapMulti & tmp_neighbors_ = data_->tmp_neighbors_;
    std::set<AASequence> & annotations_ = data_->annotations_;

  
    neighbors_.clear();
    for (NeighborMapMulti::const_iterator n_it = tmp_neighbors_.begin();
         n_it != tmp_neighbors_.end(); ++n_it)
    {
      for (std::multimap<double, GridFeature*>::const_iterator df_it =
             n_it->second.begin(); df_it != n_it->second.end(); ++df_it)
      {
        const std::set<AASequence>& current = df_it->second->getAnnotations();
        if (current.empty() || (current == annotations_))
        {
          neighbors_[n_it->first] = make_pair(df_it->first, df_it->second);
          break; // found the best element for this input map
        }
      }
    }
  }

  void QTCluster::makeSeq_table_(map<set<AASequence>, vector<double>> &seq_table) const
  {

    // get copies of members that are used in this function
    Size num_maps_ = data_->num_maps_;
    double max_distance_ = data_->max_distance_;

    // get reference on member that is used in this function
    NeighborMapMulti & tmp_neighbors_ = data_->tmp_neighbors_;

    for (NeighborMapMulti::iterator n_it = tmp_neighbors_.begin();
         n_it != tmp_neighbors_.end(); ++n_it)
    {
      Size map_index = n_it->first;
      for (NeighborListType::iterator df_it = n_it->second.begin();
          df_it != n_it->second.end(); ++df_it)
      {
        double dist = df_it->first;
        const std::set<AASequence>& current = df_it->second->getAnnotations();
        map<std::set<AASequence>, vector<double> >::iterator pos =
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
  }

  void QTCluster::finalizeCluster()
  {
    OPENMS_PRECONDITION(!data_->tmp_neighbors_.empty(),
        "Try to finalize QTCluster that was not initialized")
    OPENMS_PRECONDITION(!finalized_,
        "Try to finalize QTCluster that was not initialized")

    // calls computeQuality_ if something changed since initialization. In
    // turn, computeQuality_ calls optimizeAnnotations_ if necessary which
    // ensures that the neighbors_ hash is populated correctly.
    getQuality();

    finalized_ = true;

    data_->tmp_neighbors_.clear();
  }

  void QTCluster::initializeCluster()
  {
    OPENMS_PRECONDITION(data_->tmp_neighbors_.empty(),
        "Try to initialize QTCluster that was not finalized")
    OPENMS_PRECONDITION(finalized_,
        "Try to initialize QTCluster that was not finalized")

    finalized_ = false;

    data_->tmp_neighbors_.clear();
  }

  bool operator<(QTCluster const& q1, QTCluster const& q2)
  {
    return q1.getCurrentQuality() < q2.getCurrentQuality(); 
  }
} // namespace OpenMS
