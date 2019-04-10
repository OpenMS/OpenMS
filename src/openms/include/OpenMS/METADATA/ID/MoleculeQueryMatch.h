// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/METADATA/ID/MetaData.h>
#include <OpenMS/METADATA/ID/IdentifiedCompound.h>
#include <OpenMS/METADATA/ID/IdentifiedSequence.h>
#include <OpenMS/METADATA/PeptideHit.h> // for "PeakAnnotation"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/variant.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    // @TODO: move "PeakAnnotation" out of "PeptideHit"
    typedef std::vector<PeptideHit::PeakAnnotation> PeakAnnotations;
    typedef std::map<boost::optional<ProcessingStepRef>,
                     PeakAnnotations> PeakAnnotationSteps;

    typedef boost::variant<IdentifiedPeptideRef, IdentifiedCompoundRef,
                           IdentifiedOligoRef> IdentifiedMoleculeRef;

    /** @brief Meta data for a search hit (e.g. peptide-spectrum match).
    */
    struct MoleculeQueryMatch: public ScoredProcessingResult
    {
      IdentifiedMoleculeRef identified_molecule_ref;

      DataQueryRef data_query_ref;

      Int charge;

      // peak annotations (fragment ion matches), potentially from different
      // data processing steps:
      PeakAnnotationSteps peak_annotations;

      explicit MoleculeQueryMatch(
        IdentifiedMoleculeRef identified_molecule_ref,
        DataQueryRef data_query_ref, Int charge = 0,
        const AppliedProcessingSteps& steps_and_scores =
        AppliedProcessingSteps(),
        const PeakAnnotationSteps& peak_annotations = PeakAnnotationSteps()):
        ScoredProcessingResult(steps_and_scores),
        identified_molecule_ref(identified_molecule_ref),
        data_query_ref(data_query_ref), charge(charge),
        peak_annotations(peak_annotations)
      {
      }

      MoleculeQueryMatch(const MoleculeQueryMatch& other) = default;

      MoleculeType getMoleculeType() const
      {
        if (boost::get<IdentifiedPeptideRef>(&identified_molecule_ref))
        {
          return MoleculeType::PROTEIN;
        }
        if (boost::get<IdentifiedCompoundRef>(&identified_molecule_ref))
        {
          return MoleculeType::COMPOUND;
        }
        if (boost::get<IdentifiedOligoRef>(&identified_molecule_ref))
        {
          return MoleculeType::RNA;
        }
        return MoleculeType::SIZE_OF_MOLECULETYPE; // this shouldn't happen
      }

      IdentifiedPeptideRef getIdentifiedPeptideRef() const
      {
        if (const IdentifiedPeptideRef* ref_ptr =
            boost::get<IdentifiedPeptideRef>(&identified_molecule_ref))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not a peptide";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      IdentifiedCompoundRef getIdentifiedCompoundRef() const
      {
        if (const IdentifiedCompoundRef* ref_ptr =
            boost::get<IdentifiedCompoundRef>(&identified_molecule_ref))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not a compound";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      IdentifiedOligoRef getIdentifiedOligoRef() const
      {
        if (const IdentifiedOligoRef* ref_ptr =
            boost::get<IdentifiedOligoRef>(&identified_molecule_ref))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not an oligonucleotide";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      MoleculeQueryMatch& operator+=(const MoleculeQueryMatch& other)
      {
        ScoredProcessingResult::operator+=(other);
        if (charge == 0) charge = other.charge;
        peak_annotations.insert(other.peak_annotations.begin(),
                                other.peak_annotations.end());
        return *this;
      }
    };

    // all matches for the same data query should be consecutive!
    typedef boost::multi_index_container<
      MoleculeQueryMatch,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
          boost::multi_index::composite_key<
            MoleculeQueryMatch,
            boost::multi_index::member<MoleculeQueryMatch, DataQueryRef,
                                       &MoleculeQueryMatch::data_query_ref>,
            boost::multi_index::member<
              MoleculeQueryMatch, IdentifiedMoleculeRef,
              &MoleculeQueryMatch::identified_molecule_ref>>>>
      > MoleculeQueryMatches;
    typedef IteratorWrapper<MoleculeQueryMatches::iterator> QueryMatchRef;

  }
}
