// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <algorithm> // for set_intersection
#include <iterator>  // for inserter
#include <numeric>   // for make_pair

using std::make_pair;
using std::map;
using std::min;
using std::set;
using std::vector;

namespace OpenMS
{
  QTCluster::BulkData::BulkData(const OpenMS::GridFeature* const center_point, Size num_maps, double max_distance, Int x_coord, Int y_coord, Size id) :
      center_point_(center_point), id_(id), neighbors_(), tmp_neighbors_(), max_distance_(max_distance), num_maps_(num_maps), x_coord_(x_coord), y_coord_(y_coord), annotations_()
  {
  }

  QTCluster::QTCluster(QTCluster::BulkData* const data, bool use_IDs) : quality_(0.0), data_(data), valid_(true), changed_(false), use_IDs_(use_IDs), collect_annotations_(false), finalized_(true)
  {
    if (use_IDs)
    {
      data_->annotations_ = data_->center_point_->getAnnotations();
    }
    if (use_IDs_ && data_->center_point_->getAnnotations().size() != 1)
    {
      collect_annotations_ = true;
    }
  }

  const GridFeature* QTCluster::getCenterPoint() const
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
    OPENMS_PRECONDITION(finalized_, "Cannot perform operation on cluster that is not finalized")
    return data_->neighbors_.size() + 1; // + 1 for the center
  }

  bool QTCluster::operator<(const QTCluster& rhs) const
  {
    OPENMS_PRECONDITION(finalized_, "Cannot perform operation on cluster that is not finalized")

    return quality_ < rhs.quality_;
  }

  void QTCluster::add(const OpenMS::GridFeature* const element, double distance)
  {
    // get reference on member that is used in this function
    NeighborMapMulti& tmp_neighbors_ = data_->tmp_neighbors_;

    OPENMS_PRECONDITION(!finalized_, "Cannot perform operation on cluster that is not initialized")
    // ensure we only add compatible peptide annotations
    OPENMS_PRECONDITION(distance <= data_->max_distance_, "Distance cannot be larger than max_distance")

    Size map_index = element->getMapIndex();

    // get reference on member that is used in this function
    const OpenMS::GridFeature& center_point = *(data_->center_point_);

    // Ensure we only add compatible peptide annotations. If the cluster center
    // has an annotation, then each added neighbor should have at least one matching annotation.
    // If the center element has no annotation we add all elements
    // and select the optimal annotation later (as in the case of multiple annotations), using optimizeAnnotations_
    if (use_IDs_)
    {
      bool one_empty = (center_point.getAnnotations().empty() || element->getAnnotations().empty());
      if (!one_empty) // both are annotated
      {
        set<AASequence> intersect;
        // overlap of at least one sequence is enough
        std::set_intersection(center_point.getAnnotations().begin(), center_point.getAnnotations().end(), element->getAnnotations().begin(), element->getAnnotations().end(),
                              std::inserter(intersect, intersect.begin()));
        if (intersect.empty())
          return;
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
    // TODO This might be wrong now with multiple seqs.
    //  It might need consider the seqTable for every seq of the intersection!
    //  On the other hand this just fills data_->neighbors_ which says it only stores the BEST feature per map.
    if (map_index != center_point.getMapIndex())
    {
      NeighborMap& neighbors_ = data_->neighbors_;

      if (neighbors_.find(map_index) == neighbors_.end() || distance < neighbors_[map_index].distance)
      {
        neighbors_[map_index] = Neighbor {distance, element};
        changed_ = true;
      }
    }
  }

  QTCluster::Elements QTCluster::getElements() const
  {
    OPENMS_PRECONDITION(finalized_, "Cannot perform operation on cluster that is not finalized")

    // get the neighbors and then insert the cluster center
    Elements elements = getAllNeighbors();

    // add center point to the copy
    elements.push_back({data_->center_point_->getMapIndex(), data_->center_point_});

    // Named return value optimization (no additional copy or move when returning by value)
    return elements;
  }

  bool QTCluster::update(const Elements& removed)
  {
    OPENMS_PRECONDITION(finalized_, "Cannot perform operation on cluster that is not finalized")

    // check if the cluster center was removed:
    for (const auto& removed_element : removed)
    {
      // If center point was removed, then we are done and no more work is
      // required
      if (removed_element.feature == data_->center_point_)
      {
        this->setInvalid();
        return false;
      }
    }

    // get references on member that is used in this function
    NeighborMap& neighbors_ = data_->neighbors_;

    // update cluster contents, remove those elements we find in our cluster
    for (const auto& removed_element : removed)
    {
      NeighborMap::iterator pos = neighbors_.find(removed_element.map_index);
      if (pos == neighbors_.end())
      {
        continue; // no points from this map
      }

      const GridFeature* const current_feature = pos->second.feature;
      if (current_feature == removed_element.feature) // remove this neighbor
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
    OPENMS_PRECONDITION(finalized_, "Cannot perform operation on cluster that is finalized")

    return quality_;
  }

  void QTCluster::computeQuality_()
  {
    // ensure cluster is not finalized as we cannot call optimizeAnnotations_
    // in that case
    OPENMS_PRECONDITION(!finalized_, "Cannot perform operation on cluster that is finalized")

    // get references on member that is used in this function
    NeighborMap& neighbors_ = data_->neighbors_;

    // get copy of member that is used in this function
    double max_distance_ = data_->max_distance_;

    Size num_other = data_->num_maps_ - 1;
    double internal_distance = 0.0;
    if (!use_IDs_ || data_->center_point_->getAnnotations().size() == 1 || neighbors_.empty())
    {
      // if the cluster center is annotated with one peptide ID, the neighbors can
      // consist only of features with compatible IDs, so we don't need to
      // check again here
      for (const auto& neighbor : neighbors_)
      {
        internal_distance += neighbor.second.distance;
      }
      // add max. distance for missing cluster elements:
      internal_distance += double(num_other - neighbors_.size()) * max_distance_;
    }
    else // find the annotation that gives the best quality
    {
      internal_distance = optimizeAnnotations_();
    }

    // normalize:
    internal_distance /= num_other;
    quality_ = (max_distance_ - internal_distance) / max_distance_;
  }

  QTCluster::Elements QTCluster::getAllNeighbors() const
  {
    OPENMS_PRECONDITION(finalized_, "Cannot perform operation on cluster that is not finalized")

    // copy the important info about the neighbors
    Elements elements;
    for (const auto& neighbor : data_->neighbors_)
    {
      elements.push_back({neighbor.first, neighbor.second.feature});
    }

    // Named return value optimization (no additional copy or move when returning by value)
    return elements;
  }

  const set<AASequence>& QTCluster::getAnnotations()
  {
    if (changed_ && use_IDs_ && data_->center_point_->getAnnotations().size() != 1 && !data_->neighbors_.empty())
    {
      optimizeAnnotations_();
    }
    return data_->annotations_;
  }

  double QTCluster::optimizeAnnotations_()
  {
    OPENMS_PRECONDITION(collect_annotations_, "QTCluster::optimizeAnnotations_ should only be called if we use collect_annotations_")
    OPENMS_PRECONDITION(!data_->tmp_neighbors_.empty(), "QTCluster::optimizeAnnotations_ needs to have working tmp_neighbors_")
    OPENMS_PRECONDITION(!finalized_, "QTCluster::optimizeAnnotations_ cannot work on finalized cluster")

    // mapping: peptides -> best distance per input map
    map<AASequence, map<Size, double>> seq_table;

    makeSeqTable_(seq_table);

    // get copies of members that are used in this function
    Size num_maps_ = data_->num_maps_;
    double max_distance_ = data_->max_distance_;

    // combine annotation-specific and unspecific distances
    // (all unspecific ones are grouped as empty AASequence):
    auto unspecific = seq_table.find(AASequence());
    if (unspecific != seq_table.end())
    {
      for (auto it = seq_table.begin(); it != seq_table.end(); ++it) // OMS_CODING_TEST_EXCLUDE
      {
        if (it == unspecific)
          continue;
        // for all the maps for the "real" sequences, overwrite the distance, if an unspecific one is better or
        // add a new entry, so we can just "sum up" the distances for each seq later.
        for (const auto& mapidx_unspecdist : unspecific->second)
        {
          // try to add new entry
          auto mapidx_inserted = it->second.emplace(mapidx_unspecdist.first, mapidx_unspecdist.second);
          if (!mapidx_inserted.second) // an entry for that map idx already existed for the sequence, check minimum of both
          {
            mapidx_inserted.first->second = min(mapidx_inserted.first->second, mapidx_unspecdist.second);
          }
        }
      }
    }

    // compute distance totals -> best annotation set has smallest value:
    auto best_pos = seq_table.begin();
    double best_total = num_maps_ * max_distance_;
    for (auto it = seq_table.begin(); it != seq_table.end(); ++it) // OMS_CODING_TEST_EXCLUDE
    {
      OPENMS_PRECONDITION(num_maps_ - 1 >= it->second.size(), "num_maps bigger than map size -1 (center)");
      // init value is #missing maps times max_distance
      // above, unspecific distances were incorporated into the rest already.
      double total = std::accumulate(it->second.begin(), it->second.end(), double(num_maps_ - 1 - it->second.size()) * max_distance_,
                                     [](double val, const std::map<Size, double>::value_type& p) { return val + p.second; });
      if (total < best_total)
      {
        best_pos = it;
        best_total = total;
      }
    }

    if (best_pos != seq_table.end())
    {
      // TODO can we accumulate the union of possible annotations and set the best as "representative"?
      //  Probably in another member and function though (e.g. after finalize),
      //  since annotations_ is used in recomputeNeighbors to filter the cluster for matching
      //  features of this "best" annotation.
      //  Actually I think during creation of the consensusFeature later, all IDs of the linked features
      //  (from the original full data) are copied anyway.
      //  Then it would make sense to save the "best" annotation "distance-wise" from this algorithm, to be used during
      //  IDConflictResolution (which is based on only ID scores).
      //  OR already consider the ID scores here and make a more elaborate scoring.
      data_->annotations_ = {best_pos->first};
    }

    // only keep neighbors that fit with the best annotation!
    recomputeNeighbors_();

    return best_total;
  }

  void QTCluster::recomputeNeighbors_()
  {
    // get references on members that are used in this function
    NeighborMap& neighbors_ = data_->neighbors_;
    NeighborMapMulti& tmp_neighbors_ = data_->tmp_neighbors_;
    std::set<AASequence>& annotations_ = data_->annotations_;

    neighbors_.clear();
    for (NeighborMapMulti::const_iterator n_it = tmp_neighbors_.begin(); n_it != tmp_neighbors_.end(); ++n_it)
    {
      for (std::multimap<double, const GridFeature*>::const_iterator df_it = n_it->second.begin(); df_it != n_it->second.end(); ++df_it)
      {
        std::set<AASequence> intersect;
        const std::set<AASequence>& current = df_it->second->getAnnotations();
        std::set_intersection(current.begin(), current.end(), annotations_.begin(), annotations_.end(), std::inserter(intersect, intersect.begin()));
        // if no overlap with the re-calculated IDs in the center, do not re-add neighbor to the updated neighbors anymore.
        if (!intersect.empty() || current.empty())
        {
          neighbors_[n_it->first] = Neighbor {df_it->first, df_it->second};
          break; // found the best element for this input map
        }
      }
    }
  }

  void QTCluster::makeSeqTable_(map<AASequence, map<Size, double>>& seq_table) const
  {
    // get reference on member that is used in this function
    NeighborMapMulti& tmp_neighbors_ = data_->tmp_neighbors_;

    // for all maps contributing to this cluster
    for (NeighborMapMulti::iterator n_it = tmp_neighbors_.begin(); n_it != tmp_neighbors_.end(); ++n_it)
    {
      // for all neighbors relevant for this cluster in this map
      Size map_index = n_it->first;
      for (NeighborList::iterator df_it = n_it->second.begin(); df_it != n_it->second.end(); ++df_it)
      {
        double dist = df_it->first;
        // for all IDs/annotations of the neighboring feature (skipped if empty)
        for (const auto& current : df_it->second->getAnnotations())
        {
          auto seqit_inserted = seq_table.emplace(current, map<Size, double> {{map_index, dist}});
          // check if a minimum distance was already set for this ID
          if (!seqit_inserted.second)
          {
            // if so, check if a dist was annotated for that map_index already.
            auto distit_inserted = seqit_inserted.first->second.emplace(map_index, dist);

            // if so:
            // new dist. value for this input map
            // compare with old and set minimum
            if (!distit_inserted.second)
            {
              distit_inserted.first->second = min(dist, distit_inserted.first->second);
            }
          }
        }

        if (df_it->second->getAnnotations().empty()) // unannotated feature
        {
          auto seqit_inserted = seq_table.emplace(AASequence(), map<Size, double> {{map_index, dist}});
          // check if a minimum distance was already set for empty ID = unannotated
          if (!seqit_inserted.second)
          {
            // if so, check if a dist was annotated for that map_index already.
            auto distit_inserted = seqit_inserted.first->second.emplace(map_index, dist);

            // if so:
            // new dist. value for this input map
            // compare with old and set minimum
            if (!distit_inserted.second)
            {
              distit_inserted.first->second = min(dist, distit_inserted.first->second);
            }
          }
          // As opposed to above IDed features (which could lead to new additional annotations),
          // no need to check further here: all following (also annotation-specific) distances are worse
          // than this unspecific one, since multimap is sorted & dists are already corrected
          // with noID_penalty. If you don't want this to happen, set the penalty to one and unIDed ones
          // will always be added at the end):
          break;
        }
      }
    }
  }

  void QTCluster::finalizeCluster()
  {
    OPENMS_PRECONDITION(!finalized_, "Try to finalize QTCluster that was not initialized or already finalized")

    // calls computeQuality_ if something changed since initialization. In
    // turn, computeQuality_ calls optimizeAnnotations_ if necessary which
    // ensures that the neighbors_ hash is populated correctly.
    getQuality();

    finalized_ = true;

    data_->tmp_neighbors_.clear();
  }

  void QTCluster::initializeCluster()
  {
    OPENMS_PRECONDITION(data_->tmp_neighbors_.empty(), "Try to initialize QTCluster that was not finalized")
    OPENMS_PRECONDITION(finalized_, "Try to initialize QTCluster that was not finalized")

    finalized_ = false;

    data_->tmp_neighbors_.clear();
  }

} // namespace OpenMS
