// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MRMTRANSITIONGROUP_C
#define OPENMS_KERNEL_MRMTRANSITIONGROUP_C

#include "MRMFeature.h"

namespace OpenMS
{
  /**
  @brief The representation of a transition group that has information about
  the individual chromatograms as well as the transitions it refers to.

  This means that the MRM Transition Group establishes the mapping between the
  individual Transition (containind the meta-data) and the Chromatogram data
  points (measured data).

  Since not all the functions in OpenMS will work with MSChromatogram data
  structures, this needs to accept also MSSpectrum as a type for raw data
  storage.
  */
  template <template <typename> class SpectrumType, typename PeakType, typename TransitionType>
  class MRMTransitionGroup
  {

public:

    ///Type definitions
    //@{
    /// List of MRM Features type
    typedef std::vector<MRMFeature> MRMFeatureListType;
    /// List of Reaction Monitoring transitions (meta data) type
    //typedef OpenSWATH::LightTransition TransitionType;
    typedef std::vector<TransitionType> TransitionsType;
    /// List of Chromatograms (measured data) type
    typedef SpectrumType<PeakType> SpectrumPeakType;
    //@}

    /// Constructor
    MRMTransitionGroup()
    {
    }

    /// Copy Constructor
    MRMTransitionGroup(const MRMTransitionGroup & rhs) :
      tr_gr_id_(rhs.tr_gr_id_),
      transitions_(rhs.transitions_),
      chromatograms_(rhs.chromatograms_),
      cons_features_(rhs.cons_features_),
      chromatogram_map_(rhs.chromatogram_map_),
      transition_map_(rhs.transition_map_)
    {
    }

    /// Destructor
    virtual ~MRMTransitionGroup()
    {
    }

    MRMTransitionGroup & operator=(const MRMTransitionGroup & rhs)
    {
      if (&rhs != this)
      {
        tr_gr_id_ = rhs.tr_gr_id_;
        transitions_ = rhs.transitions_;
        chromatograms_ = rhs.chromatograms_;
        cons_features_ = rhs.cons_features_;

        transition_map_ = rhs.transition_map_;
        chromatogram_map_ = rhs.chromatogram_map_;
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

    inline const std::vector<TransitionType> & getTransitions() const
    {
      return transitions_;
    }

    inline std::vector<TransitionType> & getTransitionsMuteable()
    {
      return transitions_;
    }

    inline void addTransition(const TransitionType & transition, String key)
    {
      transitions_.push_back(transition);
      transition_map_[key] = transitions_.size() - 1;
    }

    inline const TransitionType & getTransition(String key) 
    {
      return transitions_[transition_map_[key]];
    }

    inline bool hasTransition(String key)
    {
      return transition_map_.find(key) != transition_map_.end();
    }

    inline const std::vector<SpectrumPeakType> & getChromatograms() const
    {
      return chromatograms_;
    }

    inline std::vector<SpectrumPeakType> & getChromatograms()
    {
      return chromatograms_;
    }

    inline void addChromatogram(SpectrumPeakType & chromatogram, String key)
    {
      chromatograms_.push_back(chromatogram);
      chromatogram_map_[key] = chromatograms_.size() - 1;
    }

    inline SpectrumPeakType & getChromatogram(String key)
    {
      return chromatograms_[chromatogram_map_[key]];
    }

    inline bool hasChromatogram(String key)
    {
      return chromatogram_map_.find(key) != chromatogram_map_.end();
    }

    inline const std::vector<MRMFeature> & getFeatures() const
    {
      return cons_features_;
    }

    inline std::vector<MRMFeature> & getFeaturesMuteable()
    {
      return cons_features_;
    }

    inline void addFeature(MRMFeature & feature)
    {
      cons_features_.push_back(feature);
    }

    void getLibraryIntensity(std::vector<double> & result) const
    {
      for (typename TransitionsType::const_iterator it = transitions_.begin(); it != transitions_.end(); it++)
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

protected:

    /// transition group id (peak group id)
    String tr_gr_id_;

    /// transition list
    TransitionsType transitions_;

    /// chromatogram list
    std::vector<SpectrumPeakType> chromatograms_;

    /// feature list
    MRMFeatureListType cons_features_;

    std::map<String, int> chromatogram_map_;
    std::map<String, int> transition_map_;

  };
}
#endif
