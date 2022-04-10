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
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/METADATA/ID/IdentifiedMolecule.h>
#include <OpenMS/METADATA/ID/MetaData.h>
#include <OpenMS/METADATA/ID/Observation.h>
#include <OpenMS/METADATA/PeptideHit.h> // for "PeakAnnotation"
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    // @TODO: move "PeakAnnotation" out of "PeptideHit"
    typedef std::vector<PeptideHit::PeakAnnotation> PeakAnns;

    struct PeakAnnotations : public PeakAnns {
      PeakAnnotations() : PeakAnns()
      {
      }
      PeakAnnotations(const PeakAnnotations& other) : PeakAnns(other)
      {
      }
      PeakAnnotations(const PeakAnns& other) : PeakAnns(other)
      {
      }
    };

    typedef std::map<std::optional<ProcessingStepRef>, PeakAnnotations> PeakAnnSteps;
    struct PeakAnnotationSteps : public PeakAnnSteps {
      PeakAnnotationSteps() : PeakAnnSteps()
      {
      }
      PeakAnnotationSteps(const PeakAnnotationSteps& other) : PeakAnnSteps(other)
      {
      }
      PeakAnnotationSteps(const PeakAnnSteps& other) : PeakAnnSteps(other)
      {
      }
    };


    /// Comparator for adducts
    // @TODO: this allows adducts with duplicate names, but requires different
    // sum formulas/charges - is this what we want?
    struct AdductCompare {
      bool operator()(const AdductInfo& left, const AdductInfo& right) const
      {
        return (std::make_pair(left.getCharge(), left.getEmpiricalFormula()) < std::make_pair(right.getCharge(), right.getEmpiricalFormula()));
      }
    };

    typedef std::set<AdductInfo, AdductCompare> Adds;
    struct Adducts : public Adds {
      Adducts() : Adds()
      {
      }
      Adducts(const Adducts& other) : Adds(other)
      {
      }
      Adducts(const Adds& other) : Adds(other)
      {
      }
    };


    typedef IteratorWrapper<Adducts::iterator, AdductInfo> AddRef;

    struct AdductRef : public AddRef {
      AdductRef() : AddRef()
      {
      }
      AdductRef(const AdductRef& other) : AddRef(other)
      {
      }
      AdductRef(const AddRef& other) : AddRef(other)
      {
      }
      AdductRef(std::set<OpenMS::AdductInfo, OpenMS::IdentificationDataInternal::AdductCompare>::iterator other) : AddRef(other)
      {
      }
      AdductRef operator=(const AdductRef& other)
      {
        return AddRef::operator=(other);
      }
    };


    typedef std::optional<AdductRef> AdductOpt;

    /// Representation of a search hit (e.g. peptide-spectrum match).
    struct ObservationMatch : public ScoredProcessingResult {
      IdentifiedMolecule identified_molecule_var;

      ObservationRef observation_ref;

      Int charge;

      AdductOpt adduct_opt; ///< optional reference to adduct

      // peak annotations (fragment ion matches), potentially from different
      // data processing steps:
      PeakAnnotationSteps peak_annotations;

      explicit ObservationMatch(IdentifiedMolecule identified_molecule_var, ObservationRef observation_ref, Int charge = 0, const std::optional<AdductRef>& adduct_opt = std::nullopt,
                                const AppliedProcessingSteps& steps_and_scores = AppliedProcessingSteps(), const PeakAnnotationSteps& peak_annotations = PeakAnnotationSteps()) :
          ScoredProcessingResult(steps_and_scores),
          identified_molecule_var(identified_molecule_var), observation_ref(observation_ref), charge(charge), adduct_opt(adduct_opt), peak_annotations(peak_annotations)
      {
      }

      ObservationMatch(const ObservationMatch&) = default;

      // ONLY FOR PYOPENMS
      ObservationMatch();

      ObservationMatch& merge(const ObservationMatch& other)
      {
        ScoredProcessingResult::merge(other);
        if (charge == 0)
        {
          charge = other.charge;
        }
        else if (charge != other.charge)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Trying to overwrite ObservationMatch charge with conflicting value.", String(charge));
        }

        if (!adduct_opt)
        {
          adduct_opt = other.adduct_opt;
        }
        else if (adduct_opt != other.adduct_opt)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Trying to overwrite ObservationMatch adduct_opt with conflicting value.", (*adduct_opt)->getName());
        }

        peak_annotations.insert(other.peak_annotations.begin(), other.peak_annotations.end());
        return *this;
      }
    };

    // all matches for the same observation should be consecutive, so make sure
    // the observation is used as the first member in the composite key:
    typedef boost::multi_index_container<ObservationMatch, boost::multi_index::indexed_by<boost::multi_index::ordered_unique<boost::multi_index::composite_key<
                                                             ObservationMatch, boost::multi_index::member<ObservationMatch, ObservationRef, &ObservationMatch::observation_ref>,
                                                             boost::multi_index::member<ObservationMatch, IdentifiedMolecule, &ObservationMatch::identified_molecule_var>,
                                                             boost::multi_index::member<ObservationMatch, AdductOpt, &ObservationMatch::adduct_opt>>>>>
      ObsMatches;

    struct ObservationMatches : public ObsMatches {
      ObservationMatches() : ObsMatches()
      {
      }
      ObservationMatches(const ObservationMatches& other) : ObsMatches(other)
      {
      }
      ObservationMatches(const ObsMatches& other) : ObsMatches(other)
      {
      }
    };

    typedef IteratorWrapper<ObservationMatches::iterator, ObservationMatch> ObsMatref;

    struct ObservationMatchRef : public ObsMatref {
      ObservationMatchRef() : ObsMatref()
      {
      }
      ObservationMatchRef(const ObservationMatchRef& other) : ObsMatref(other)
      {
      }
      ObservationMatchRef(const ObsMatref& other) : ObsMatref(other)
      {
      }
      ObservationMatchRef(
        const boost::multi_index::detail::bidir_node_iterator<boost::multi_index::detail::ordered_index_node<
          boost::multi_index::detail::null_augment_policy,
          boost::multi_index::detail::index_node_base<OpenMS::IdentificationDataInternal::ObservationMatch, std::allocator<OpenMS::IdentificationDataInternal::ObservationMatch>>>>& other) :
          ObsMatref(other)
      {
      }
      ObservationMatchRef operator=(const ObservationMatchRef& other)
      {
        return ObsMatref::operator=(other);
      }
    };

  } // namespace IdentificationDataInternal
} // namespace OpenMS
