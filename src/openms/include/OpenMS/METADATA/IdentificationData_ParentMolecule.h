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
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IDENTIFICATIONDATA_PARENTMOLECULE_H
#define OPENMS_METADATA_IDENTIFICATIONDATA_PARENTMOLECULE_H

#include <OpenMS/METADATA/IdentificationData_MetaData.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Base class for data with scores and processing steps (and meta info)
    struct ScoredProcessingResult: public MetaInfoInterface
    {
      ScoreList scores;

      // @TODO: use a "boost::multi_index_container" here for efficient look-up?
      std::vector<ProcessingStepRef> processing_step_refs;

      ScoredProcessingResult& operator+=(const ScoredProcessingResult& other)
      {
        // merge processing steps:
        for (auto step_ref : other.processing_step_refs)
        {
          if (std::find(processing_step_refs.begin(),
                        processing_step_refs.end(), step_ref) ==
              processing_step_refs.end())
          {
            processing_step_refs.push_back(step_ref);
          }
        }
        // merge scores:
        for (auto score_pair : other.scores)
        {
          // @TODO: should we overwrite scores?
          if (std::find(scores.begin(), scores.end(), score_pair) ==
              scores.end())
          {
            scores.push_back(score_pair);
          }
        }
        // merge meta info:
        std::vector<UInt> keys;
        other.getKeys(keys);
        for (const UInt key : keys)
        {
          // @TODO: should we overwrite meta values?
          if (!metaValueExists(key))
          {
            setMetaValue(key, other.getMetaValue(key));
          }
        }

        return *this;
      }

      std::pair<double, bool> getScore(ScoreTypeRef score_ref) const
      {
        // give priority to "later" scores in the list:
        for (ScoreList::const_reverse_iterator it = scores.rbegin();
             it != scores.rend(); ++it)
        {
          if (it->first == score_ref) return std::make_pair(it->second, true);
        }
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), false);
      }

    protected:
      explicit ScoredProcessingResult(
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepRef>& processing_step_refs =
        std::vector<ProcessingStepRef>()):
        scores(scores), processing_step_refs(processing_step_refs)
      {
      }

      ScoredProcessingResult(const ScoredProcessingResult& other) = default;
    };


    /*!
      Representation of a parent molecule that is identified only indirectly (e.g. a protein).
    */
    struct ParentMolecule: public ScoredProcessingResult
    {
      String accession;

      enum MoleculeType molecule_type;

      // @TODO: if there are modifications in the sequence, "sequence.size()"
      // etc. will be misleading!
      String sequence;

      String description;

      double coverage; //< sequence coverage as a fraction between 0 and 1

      bool is_decoy;

      explicit ParentMolecule(
        const String& accession,
        MoleculeType molecule_type = MoleculeType::PROTEIN,
        const String& sequence = "", const String& description = "",
        double coverage = 0.0, bool is_decoy = false,
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepRef>& processing_step_refs =
        std::vector<ProcessingStepRef>()):
        ScoredProcessingResult(scores, processing_step_refs),
        accession(accession), molecule_type(molecule_type), sequence(sequence),
        description(description), coverage(coverage), is_decoy(is_decoy)
      {
      }

      ParentMolecule(const ParentMolecule& other) = default;

      ParentMolecule& operator+=(const ParentMolecule& other)
      {
        ScoredProcessingResult::operator+=(other);
        if (sequence.empty()) sequence = other.sequence;
        if (description.empty()) description = other.description;
        if (!is_decoy) is_decoy = other.is_decoy; // believe it when it's set
        // @TODO: what about coverage? (not reliable if we're merging data)

        return *this;
      }
    };

    // parent molecules indexed by their accessions:
    // @TODO: allow querying/iterating over proteins and RNAs separately
    typedef boost::multi_index_container<
      ParentMolecule,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          ParentMolecule, String, &ParentMolecule::accession>>>
      > ParentMolecules;
    typedef IteratorWrapper<ParentMolecules::iterator> ParentMoleculeRef;


    /*!
      Group of ambiguously identified parent molecules
    */
    struct ParentMoleculeGroup
    {
      ScoreList scores;
      // @TODO: does this need a "leader" or some such?
      std::set<ParentMoleculeRef> parent_molecule_refs;
    };

    typedef boost::multi_index_container<
      ParentMoleculeGroup,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
        boost::multi_index::member<
          ParentMoleculeGroup, std::set<ParentMoleculeRef>,
          &ParentMoleculeGroup::parent_molecule_refs>>>
      > ParentMoleculeGroups;
    typedef IteratorWrapper<ParentMoleculeGroups::iterator> ParentGroupRef;

    /*!
      Set of groups of ambiguously identified parent molecules (e.g. results of running a protein inference algorithm)
    */
    struct ParentMoleculeGrouping
    {
      String label; // @TODO: use "label" as a uniqueness constraint?
      std::vector<ProcessingStepRef> processing_step_refs;
      ParentMoleculeGroups groups;
    };

    typedef std::vector<ParentMoleculeGrouping> ParentMoleculeGroupings;


    /*!
      Meta data for the association between an identified molecule (e.g. peptide) and a parent molecule (e.g. protein).
    */
    struct MoleculeParentMatch: public MetaInfoInterface
    {
      // in extraordinary cases (e.g. database searches that allow insertions/
      // deletions), the length of the identified molecule may differ from the
      // length of the subsequence in the parent; therefore, store "end_pos":
      Size start_pos, end_pos;

      // String instead of char so modified residues can be represented:
      String left_neighbor, right_neighbor; // neighboring sequence elements

      static const Size UNKNOWN_POSITION; // = Size(-1)
      static const char UNKNOWN_NEIGHBOR; // = 'X'
      static const char LEFT_TERMINUS; // = '['
      static const char RIGHT_TERMINUS; // = ']'

      explicit MoleculeParentMatch(Size start_pos = UNKNOWN_POSITION,
                                   Size end_pos = UNKNOWN_POSITION,
                                   String left_neighbor = UNKNOWN_NEIGHBOR,
                                   String right_neighbor = UNKNOWN_NEIGHBOR):
        start_pos(start_pos), end_pos(end_pos), left_neighbor(left_neighbor),
        right_neighbor(right_neighbor)
      {
      }

      bool operator<(const MoleculeParentMatch& other) const
      {
        // positions determine neighbors - no need to compare those:
        return (std::tie(start_pos, end_pos) <
                std::tie(other.start_pos, other.end_pos));
      }

      bool operator==(const MoleculeParentMatch& other) const
      {
        // positions determine neighbors - no need to compare those:
        return (std::tie(start_pos, end_pos) ==
                std::tie(other.start_pos, other.end_pos));
      }

      bool hasValidPositions(Size molecule_length = 0, Size parent_length = 0) const
      {
        if ((start_pos == UNKNOWN_POSITION) || (end_pos == UNKNOWN_POSITION))
        {
          return false;
        }
        if (end_pos < start_pos) return false;
        if (molecule_length && (end_pos - start_pos + 1 != molecule_length))
        {
          return false;
        }
        if (parent_length && (end_pos >= parent_length)) return false;
        return true;
      }
    };

    // mapping: parent molecule -> match information
    typedef std::map<ParentMoleculeRef,
                     std::set<MoleculeParentMatch>> ParentMatches;
  }
}

#endif
