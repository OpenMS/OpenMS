// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/ObservationMatch.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /** @brief: Group of related (co-identified) input matches

      E.g. for cross-linking data or multiplexed spectra.
    */
    struct ObservationMatchGroup: public ScoredProcessingResult
    {
      std::set<ObservationMatchRef> observation_match_refs;

      bool allSameMolecule() const
      {
        // @TODO: return true or false for the empty set?
        if (observation_match_refs.size() <= 1) return true;
        const IdentifiedMolecule var =
          (*observation_match_refs.begin())->identified_molecule_var;
        for (auto it = ++observation_match_refs.begin();
             it != observation_match_refs.end(); ++it)
        {
          if (!((*it)->identified_molecule_var == var)) return false;
        }
        return true;
      }

      bool allSameQuery() const
      {
        // @TODO: return true or false for the empty set?
        if (observation_match_refs.size() <= 1) return true;
        ObservationRef ref = (*observation_match_refs.begin())->observation_ref;
        for (auto it = ++observation_match_refs.begin();
             it != observation_match_refs.end(); ++it)
        {
          if ((*it)->observation_ref != ref) return false;
        }
        return true;
      }

      bool operator==(const ObservationMatchGroup& rhs) const
      {
        return ((rhs.observation_match_refs == observation_match_refs) &&
                (rhs.steps_and_scores == steps_and_scores));
      }

      bool operator!=(const ObservationMatchGroup& rhs) const
      {
        return !operator==(rhs);
      }
    };

    typedef boost::multi_index_container<
      ObservationMatchGroup,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
          boost::multi_index::member<ObservationMatchGroup, std::set<ObservationMatchRef>,
                                     &ObservationMatchGroup::observation_match_refs>>>
      > ObservationMatchGroups;
    typedef IteratorWrapper<ObservationMatchGroups::iterator> MatchGroupRef;
  }
}
