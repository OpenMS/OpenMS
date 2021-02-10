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
        DataQueryRef data_query_ref, Int m_charge = 0,
        const AppliedProcessingSteps& steps_and_scores = AppliedProcessingSteps(),
        const PeakAnnotationSteps& peak_annotations = PeakAnnotationSteps()
      )
        : ScoredProcessingResult(steps_and_scores),
          identified_molecule_ref(identified_molecule_ref),
          data_query_ref(data_query_ref), charge(m_charge),
          peak_annotations(peak_annotations)
      {
      }

      MoleculeQueryMatch(const MoleculeQueryMatch&) = default;

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

    /** @brief Stable reference to a query match.
     *  This reference stays valid (as opposed to using e.g., iterators and memory addresses)
     *  if an ID datastructure is copied or stored/loaded back from disc.       
    */
    struct StableQueryMatchRef
    {
      std::string basename; // filename with extension (without path)
      std::string native_id; // the spectrum or consensus/feature native id
      std::set<std::pair<MoleculeType, std::string>> identified_molecules; // identified molecule (peptide/RNA/compound) and identifier (sequence or database id)

      // Construct from string representation e.g. "myfile1.mzML|spectrum=123|P_DEPIANGER|P_TESTPEPTIDER"
      StableQueryMatchRef(std::string s)
      {
        size_t a = s.find("|");
        if (a == string::npos) 
        {
          std::string msg = "Invalid string. Conversion to stable reference not possible.";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
        }
        basename = s.substr(0, a)

        size_t b = s.find("|", a + 1);
        if (b == string::npos)
        {
          std::string msg = "Invalid string. Conversion to stable reference not possible.";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);

        }
        native_id = s.substr(a + 1, b - a - 1);

        if (b == s.size() - 1) // end of native id was at end of string
        {
          // we have a valid query match to reference (no molecules have been identified->empty set)
          return;
        }

        // parse set of molecules
        size_t c = s.find("|", b + 1); // let c point to next | (if present)
        while (c <= std::string::npos)
        {
          if (s[b + 1] == 'P')
          {
            identified_molecules.insert({MoleculeType::PROTEIN, s.substr(b + 3, c - b - 3)});
          }
          else if (s[b + 1] == 'C')
          {
            identified_molecules.insert({MoleculeType::COMPOUND, s.substr(b + 3, c - b - 3)});
          }
          else if (s[b + 1] == 'O')
          {
            identified_molecules.insert({MoleculeType::OLIGO, s.substr(b + 3, c - b - 3)});
          }
          else 
          {
            std::string msg = "Invalid character for molecule encoding found. Conversion to stable reference not possible.";
            throw Exception::IllegalArgument(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, msg);           
          }

          if (c == std::string::npos) return;

          b = c;
    	    c = s.find("|", b + 1);
        }
      }

      // Convert to string representation e.g. "myfile2.mzML|spectrum=4|C_HMDB:23433"
      // format: "basename|native_id|" is mandatory, 
      //         followed by one or more molecule type (encoded by "P_", "C_", "O_" and molecule identifier.
      //         Molecule identifier are sequences for Proteins and Oligos or database identifier for Compounds
      std::string toStdString() const
      {
        std::string s = basename + '|' + native_id + '|';
        for (const auto& m : identified_molecules)
        {
          switch (m.first)
          {
            case MoleculeType::PROTEIN:
              s += "P_";
            break;
            case MoleculeType::COMPOUND:
              s += "C_";
            break;
            case MoleculeType::OLIGO:
              s += "O_";
            break;            
          }

          if (&m != &identified_molecules.back()) s += m.second + "|";
        }
        return s;
      }      
    };

  }
}
