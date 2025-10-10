// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <boost/numeric/conversion/cast.hpp>

namespace OpenMS
{

  /**
    @brief The representation of a group of transitions in a targeted proteomics experiment.

    The transition group carries information about the transitions (assays), the
    individual chromatograms as well as features found on these chromatograms.

    On the one hand, the MRMTransitionGroup provides a convenient way to store
    the mapping between the individual transitions (containing the meta-data) and
    the actual chromatographic data points (measured data) relating to it. In
    addition, the structure allows storage of features found (regions of the
    chromatograms) where a potential elution peak was detected (see MRMFeature).
    Note that these features are usually found on the full collection of
    chromatograms and therefore relate to the whole collection of chromatograms.

    Note that for the data structure to be consistent, it needs to have the
    same identifiers for the chromatograms as well as for the transitions.

    Since not all the functions in OpenMS will work with MSChromatogram data
    structures, this needs to accept also MSSpectrum as a type for raw data
    storage.

  */
  template <typename ChromatogramType, typename TransitionType>
  class MRMTransitionGroup
  {

public:

    ///Type definitions
    //@{
    /// List of MRM Features type
    typedef std::vector<MRMFeature> MRMFeatureListType;
    /// List of Reaction Monitoring transitions (meta data) type
    typedef std::vector<TransitionType> TransitionsType;
    /// Peak type
    typedef typename ChromatogramType::PeakType PeakType;
    //@}

    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor
    MRMTransitionGroup()
    {
    }

    /// Copy Constructor
    MRMTransitionGroup(const MRMTransitionGroup & rhs) :
      tr_gr_id_(rhs.tr_gr_id_),
      transitions_(rhs.transitions_),
      chromatograms_(rhs.chromatograms_),
      precursor_chromatograms_(rhs.precursor_chromatograms_),
      mrm_features_(rhs.mrm_features_),
      chromatogram_map_(rhs.chromatogram_map_),
      precursor_chromatogram_map_(rhs.precursor_chromatogram_map_),
      transition_map_(rhs.transition_map_)
    {
    }

    /// Destructor
    virtual ~MRMTransitionGroup()
    {
    }
    //@}

    MRMTransitionGroup & operator=(const MRMTransitionGroup & rhs)
    {
      if (&rhs != this)
      {
        tr_gr_id_ = rhs.tr_gr_id_;
        transitions_ = rhs.transitions_;
        chromatograms_ = rhs.chromatograms_;
        precursor_chromatograms_ = rhs.precursor_chromatograms_;
        mrm_features_ = rhs.mrm_features_;
        transition_map_ = rhs.transition_map_;
        chromatogram_map_ = rhs.chromatogram_map_;
        precursor_chromatogram_map_ = rhs.precursor_chromatogram_map_;
      }
      return *this;
    }

    inline Size size() const
    {
      return chromatograms_.size();
    }

    inline const String & getTransitionGroupID() const
    {
      return tr_gr_id_;
    }

    inline void setTransitionGroupID(const String & tr_gr_id)
    {
      tr_gr_id_ = tr_gr_id;
    }

    /// @name Transition access
    //@{
    inline const std::vector<TransitionType> & getTransitions() const
    {
      return transitions_;
    }

    inline std::vector<TransitionType> & getTransitionsMuteable()
    {
      return transitions_;
    }

    /** Add a transition
     *
     * Transitions are internally mapped using their nativeID,
     * i.e. TransitionType::getNativeID.
     *
     * When querying for a transition, make sure to use this key.
     */
    inline void addTransition(const TransitionType& transition, const String& key)
    {
      // store the index where to find the transition, using the key for lookup
      auto result = transition_map_.emplace(key, int(transitions_.size()));
      if (!result.second) // ouch: key was already used
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal error: Transition with nativeID was already present!", key);
      }
      transitions_.push_back(transition);
    }

    inline bool hasTransition(const String& key) const
    {
      return transition_map_.find(key) != transition_map_.end();
    }

    inline const TransitionType& getTransition(const String& key)
    {
      OPENMS_PRECONDITION(hasTransition(key), "Cannot retrieve transitions that does not exist")
      OPENMS_PRECONDITION(transitions_.size() > (size_t)transition_map_[key], "Mapping needs to be accurate")
      return transitions_[transition_map_[key]];
    }
    //@}


    /// @name (Fragment ion) chromatogram access
    //@{
    inline std::vector<ChromatogramType>& getChromatograms()
    {
      return chromatograms_;
    }

    inline const std::vector<ChromatogramType>& getChromatograms() const
    {
      return chromatograms_;
    }

    /** Add a chromatogram 
     *
     * Chromatograms are internally mapped using the provided key.
     * The ChromatogramType::getNativeID is a good choice.
     *
     * When querying for a chromatogram, make sure to use this key.
     */
    inline void addChromatogram(const ChromatogramType& chromatogram, const String& key)
    {
      // store the index where to find the chromatogram, using the key for lookup
      auto result = chromatogram_map_.emplace(key, int(chromatograms_.size()));
      if (!result.second) // ouch: key was already used
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal error: Chromatogram with nativeID was already present!", key);
      }
      chromatograms_.push_back(chromatogram);
    }

    inline bool hasChromatogram(const String& key) const
    {
      return chromatogram_map_.find(key) != chromatogram_map_.end();
    }

    inline ChromatogramType& getChromatogram(const String& key)
    {
      OPENMS_PRECONDITION(hasChromatogram(key), "Cannot retrieve chromatogram that does not exist")
      OPENMS_PRECONDITION(chromatograms_.size() > (size_t)chromatogram_map_[key], "Mapping needs to be accurate")
      return chromatograms_[chromatogram_map_.at(key)];
    }

    inline const ChromatogramType& getChromatogram(const String& key) const
    {
      OPENMS_PRECONDITION(hasChromatogram(key), "Cannot retrieve chromatogram that does not exist")
      OPENMS_PRECONDITION(chromatograms_.size() > (size_t)chromatogram_map_.at(key), "Mapping needs to be accurate")
      return chromatograms_[chromatogram_map_.at(key)];
    }
    //@}


    /// @name (Precursor ion) chromatogram access
    //@{
    inline std::vector<ChromatogramType>& getPrecursorChromatograms()
    {
      return precursor_chromatograms_;
    }

    inline const std::vector<ChromatogramType>& getPrecursorChromatograms() const
    {
      return precursor_chromatograms_;
    }

    /** Add a precursor chromatogram (extracted from an MS1 map)
     *
     * Precursor chromatograms are internally mapped using the provided key,
     * i.e. ChromatogramType::getNativeID.
     *
     * When querying for a chromatogram, make sure to use this key.
     *
     * @param chromatogram Chromatographic traces from the MS1 map to be added
     * @param key Unique identifier of the chromatogram, e.g. its nativeID
     */
    inline void addPrecursorChromatogram(const ChromatogramType& chromatogram, const String& key)
    {
      // store the index where to find the chromatogram, using the key for lookup
      auto result = precursor_chromatogram_map_.emplace(key, int(precursor_chromatograms_.size()));
      if (!result.second) // ouch: key was already used
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal error: Chromatogram with nativeID was already present!", key);
      }
      precursor_chromatograms_.push_back(chromatogram);
    }

    inline bool hasPrecursorChromatogram(const String& key) const
    {
      return precursor_chromatogram_map_.find(key) != precursor_chromatogram_map_.end();
    }

    inline ChromatogramType & getPrecursorChromatogram(const String& key)
    {
      OPENMS_PRECONDITION(hasPrecursorChromatogram(key), "Cannot retrieve precursor chromatogram that does not exist")
      OPENMS_PRECONDITION(precursor_chromatograms_.size() > (size_t)precursor_chromatogram_map_.at(key), "Mapping needs to be accurate")
      return precursor_chromatograms_[precursor_chromatogram_map_.at(key)];
    }

    inline const ChromatogramType & getPrecursorChromatogram(const String& key) const
    {
      OPENMS_PRECONDITION(hasPrecursorChromatogram(key), "Cannot retrieve precursor chromatogram that does not exist")
      OPENMS_PRECONDITION(precursor_chromatograms_.size() > (size_t)precursor_chromatogram_map_.at(key), "Mapping needs to be accurate")
      return precursor_chromatograms_[precursor_chromatogram_map_.at(key)];
    }
    //@}


    /// @name MRM feature access (positions in RT where a peak was found across all chromatograms)
    //@{
    inline const std::vector<MRMFeature> & getFeatures() const
    {
      return mrm_features_;
    }

    inline std::vector<MRMFeature> & getFeaturesMuteable()
    {
      return mrm_features_;
    }

    inline void addFeature(const MRMFeature & feature)
    {
      mrm_features_.push_back(feature);
    }

    inline void addFeature(MRMFeature && feature)
    {
      mrm_features_.push_back(std::move(feature));
    }
    //@}


    /// @name Helper functions
    //@{

    /// Check whether internal state is consistent, e.g. same number of chromatograms and transitions are present (no runtime overhead in release mode)
    inline bool isInternallyConsistent() const
    {
      OPENMS_PRECONDITION(transitions_.size() == chromatograms_.size(), "Same number of transitions as chromatograms are required")
      OPENMS_PRECONDITION(transition_map_.size() == chromatogram_map_.size(), "Same number of transitions as chromatograms mappings are required")
      OPENMS_PRECONDITION(isMappingConsistent_(), "Mapping needs to be consistent")
      return true;
    }

    /// Ensure that chromatogram native ids match their keys in the map
    inline bool chromatogramIdsMatch() const
    {
      for (std::map<String, int>::const_iterator it = chromatogram_map_.begin(); it != chromatogram_map_.end(); it++)
      {
        if (getChromatogram(it->first).getNativeID() != it->first)
        {
          return false;
        }
      }
      for (std::map<String, int>::const_iterator it = precursor_chromatogram_map_.begin(); it != precursor_chromatogram_map_.end(); it++)
      {
        if (getPrecursorChromatogram(it->first).getNativeID() != it->first)
        {
          return false;
        }
      }
      return true;
    }

    void getLibraryIntensity(std::vector<double> & result) const
    {
      for (typename TransitionsType::const_iterator it = transitions_.begin(); it != transitions_.end(); ++it)
      {
        result.push_back(it->getLibraryIntensity());
      }
      for (Size i = 0; i < result.size(); i++)
      {
        // the library intensity should never be below zero
        if (result[i] < 0.0)
        {
          result[i] = 0.0;
        }
      }
    }

    MRMTransitionGroup subset(std::vector<std::string> tr_ids) const
    {
      MRMTransitionGroup transition_group_subset;
      transition_group_subset.setTransitionGroupID(tr_gr_id_);

      for (const auto& tr : transitions_)
      {
        if (std::find(tr_ids.begin(), tr_ids.end(), tr.getNativeID()) != tr_ids.end())
        {
          if (this->hasTransition(tr.getNativeID()))
          {
            transition_group_subset.addTransition(tr, tr.getNativeID());
          }
          if (this->hasChromatogram(tr.getNativeID()))
          {
            transition_group_subset.addChromatogram(chromatograms_[chromatogram_map_.at(tr.getNativeID())], tr.getNativeID());
          }
        }
      }

      for (const auto& pc : precursor_chromatograms_)
      {
        // add precursor chromatograms if present
        transition_group_subset.addPrecursorChromatogram(pc, pc.getNativeID());
      }

      for (const auto& tgf : mrm_features_)
      {
        MRMFeature mf;
        mf.setIntensity(tgf.getIntensity());
        mf.setRT(tgf.getRT());
        mf.MetaInfoInterface::operator=(tgf);

        for (const auto& tr : transitions_)
        {
          if (std::find(tr_ids.begin(), tr_ids.end(), tr.getNativeID()) != tr_ids.end())
          {
            mf.addFeature(tgf.getFeature(tr.getNativeID()),tr.getNativeID());
          }
        }
        std::vector<String> pf_ids;
        tgf.getPrecursorFeatureIDs(pf_ids);
        for (const auto& pf_id : pf_ids)
        {
          mf.addPrecursorFeature(tgf.getPrecursorFeature(pf_id), pf_id);
        }
        transition_group_subset.addFeature(mf);
      }

      return transition_group_subset;
    }

    MRMTransitionGroup subsetDependent(std::vector<std::string> tr_ids) const
    {
      MRMTransitionGroup transition_group_subset;
      transition_group_subset.setTransitionGroupID(tr_gr_id_);

      for (typename TransitionsType::const_iterator tr_it = transitions_.begin(); tr_it != transitions_.end(); ++tr_it)
      {
        if (std::find(tr_ids.begin(), tr_ids.end(), tr_it->getNativeID()) != tr_ids.end())
        {
          transition_group_subset.addTransition(*tr_it, tr_it->getNativeID());
          transition_group_subset.addChromatogram(chromatograms_[chromatogram_map_.at(tr_it->getNativeID())], tr_it->getNativeID());
        }
      }

      for (std::vector< MRMFeature >::const_iterator tgf_it = mrm_features_.begin(); tgf_it != mrm_features_.end(); ++tgf_it)
      {
        transition_group_subset.addFeature(*tgf_it);
      }

      return transition_group_subset;
    }
    //@}

    /**
      @brief Returns the best feature by overall quality.

      For the given transition group, return the best feature as determined by
      the overall quality score. Requires the feature list to not be empty.

    */
    const MRMFeature& getBestFeature() const
    {
      OPENMS_PRECONDITION(!getFeatures().empty(), "Cannot get best feature for empty transition group")

      // Find the feature with the highest score
      Size bestf = 0;
      double highest_score = getFeatures()[0].getOverallQuality();
      for (Size it = 0; it < getFeatures().size(); it++)
      {
        if (getFeatures()[it].getOverallQuality() > highest_score)
        {
          bestf = it;
          highest_score = getFeatures()[it].getOverallQuality();
        }
      }
      return getFeatures()[bestf];
    }

protected:

    /// Checks that the mapping between chromatograms and transitions is consistent
    bool isMappingConsistent_() const
    {
      if (transition_map_.size() != chromatogram_map_.size()) 
      {
        return false;
      }
      for (std::map<String, int>::const_iterator it = chromatogram_map_.begin(); it != chromatogram_map_.end(); it++)
      {
        if (!hasTransition(it->first)) 
        {
          return false;
        }
      }
      return true;
    }

    /// transition group id (peak group id)
    String tr_gr_id_;

    /// transition list
    TransitionsType transitions_;

    /// chromatogram list
    std::vector<ChromatogramType> chromatograms_;

    /// precursor chromatogram list
    std::vector<ChromatogramType> precursor_chromatograms_;

    /// feature list
    MRMFeatureListType mrm_features_;

    std::map<String, int> chromatogram_map_;
    std::map<String, int> precursor_chromatogram_map_;
    std::map<String, int> transition_map_;

  };
}

