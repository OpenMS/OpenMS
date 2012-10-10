// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MRMTRANSITIONGROUP_H
#define OPENMS_KERNEL_MRMTRANSITIONGROUP_H

#include <OpenMS/KERNEL/MRMFeature.h>

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
  template <typename SpectrumType, typename TransitionType>
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
    typedef typename SpectrumType::PeakType PeakType;
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

    inline const std::vector<SpectrumType> & getChromatograms() const
    {
      return chromatograms_;
    }

    inline std::vector<SpectrumType> & getChromatograms()
    {
      return chromatograms_;
    }

    inline void addChromatogram(SpectrumType & chromatogram, String key)
    {
      chromatograms_.push_back(chromatogram);
      chromatogram_map_[key] = chromatograms_.size() - 1;
    }

    inline SpectrumType & getChromatogram(String key)
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
    std::vector<SpectrumType> chromatograms_;

    /// feature list
    MRMFeatureListType cons_features_;

    std::map<String, int> chromatogram_map_;
    std::map<String, int> transition_map_;

  };
}
#endif
