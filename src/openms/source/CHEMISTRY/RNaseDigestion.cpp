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
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>

using namespace std;

namespace OpenMS
{
  void RNaseDigestion::setEnzyme(const DigestionEnzyme* enzyme)
  {
    EnzymaticDigestion::setEnzyme(enzyme);

    const DigestionEnzymeRNA* rnase =
      dynamic_cast<const DigestionEnzymeRNA*>(enzyme_);
    String five_prime_code = rnase->getFivePrimeGain();
    if (five_prime_code == "p") five_prime_code = "5'-p";
    String three_prime_code = rnase->getThreePrimeGain();
    if (three_prime_code == "p") three_prime_code = "3'-p";

    static RibonucleotideDB* ribo_db = RibonucleotideDB::getInstance();
    five_prime_gain_ = five_prime_code.empty() ?
      nullptr : ribo_db->getRibonucleotide(five_prime_code);
    three_prime_gain_ = three_prime_code.empty() ?
      nullptr : ribo_db->getRibonucleotide(three_prime_code);

    cuts_after_regex_.assign(rnase->getCutsAfterRegEx());
    cuts_before_regex_.assign(rnase->getCutsBeforeRegEx());
  }


  void RNaseDigestion::setEnzyme(const String& enzyme_name)
  {
    setEnzyme(RNaseDB::getInstance()->getEnzyme(enzyme_name));
  }


  vector<pair<Size, Size>> RNaseDigestion::getFragmentPositions_(
    const NASequence& rna, Size min_length, Size max_length) const
  {
    if (min_length == 0) min_length = 1;
    if ((max_length == 0) || (max_length > rna.size()))
    {
      max_length = rna.size();
    }

    vector<pair<Size, Size>> result;
    if (enzyme_->getName() == NoCleavage) // no cleavage
    {
      Size length = rna.size();
      if ((length >= min_length) && (length <= max_length))
      {
        result.emplace_back(0, length);
      }
    }
    else if (enzyme_->getName() == UnspecificCleavage) // unspecific cleavage
    {
      result.reserve(rna.size() * (max_length - min_length + 1));
      for (Size i = 0; i <= rna.size() - min_length; ++i)
      {
        const Size right = std::min(i + max_length, rna.size());
        for (Size j = i + min_length; j <= right; ++j)
        {
          result.emplace_back(i, j - i);
        }
      }
    }
    else // proper enzyme cleavage
    {
      vector<Size> fragment_pos(1, 0);
      for (Size i = 1; i < rna.size(); ++i)
      {
        if (boost::regex_search(rna[i - 1]->getCode(), cuts_after_regex_) &&
            boost::regex_search(rna[i]->getCode(), cuts_before_regex_))
        {
          fragment_pos.push_back(i);
        }
      }
      fragment_pos.push_back(rna.size());

      // "fragment_pos" has at least two elements (zero and "rna.size()"):
      for (Size start_it = 0; start_it < fragment_pos.size() - 1; ++start_it)
      {
        Size start_pos = fragment_pos[start_it];
        for (Size offset = 0; offset <= missed_cleavages_; ++offset)
        {
          Size end_it = start_it + offset + 1;
          if (end_it >= fragment_pos.size()) break;
          Size end_pos = fragment_pos[end_it];

          Size length = end_pos - start_pos;
          if ((length >= min_length) && (length <= max_length))
          {
            result.emplace_back(start_pos, length);
          }
        }
      }
    }

    return result;
  }

  void RNaseDigestion::digest(const NASequence& rna, vector<NASequence>& output,
                              Size min_length, Size max_length) const
  {
    output.clear();
    if (rna.empty()) return;

    vector<pair<Size, Size>> positions = getFragmentPositions_(rna, min_length,
                                                               max_length);

    for (const auto& pos : positions)
    {
      NASequence fragment = rna.getSubsequence(pos.first, pos.second);
      if (pos.first > 0)
      {
        fragment.setFivePrimeMod(five_prime_gain_);
      }
      if (pos.first + pos.second < rna.size())
      {
        fragment.setThreePrimeMod(three_prime_gain_);
      }
      output.push_back(fragment);
    }
  }


  void RNaseDigestion::digest(IdentificationData& id_data, Size min_length,
                              Size max_length) const
  {
    for (IdentificationData::ParentMoleculeRef parent_ref =
           id_data.getParentMolecules().begin(); parent_ref !=
           id_data.getParentMolecules().end(); ++parent_ref)
    {
      if (parent_ref->molecule_type != IdentificationData::MoleculeType::RNA)
      {
        continue;
      }

      NASequence rna = NASequence::fromString(parent_ref->sequence);
      vector<pair<Size, Size>> positions =
        getFragmentPositions_(rna, min_length, max_length);

      for (const auto& pos : positions)
      {
        NASequence fragment = rna.getSubsequence(pos.first, pos.second);
        if (pos.first > 0)
        {
          fragment.setFivePrimeMod(five_prime_gain_);
        }
        if (pos.first + pos.second < rna.size())
        {
          fragment.setThreePrimeMod(three_prime_gain_);
        }
        IdentificationData::IdentifiedOligo oligo(fragment);
        Size end_pos = pos.first + pos.second; // past-the-end position!
        IdentificationData::MoleculeParentMatch match(pos.first, end_pos - 1);
        match.left_neighbor = (pos.first > 0) ? rna[pos.first - 1]->getCode() :
          IdentificationData::MoleculeParentMatch::LEFT_TERMINUS;
        match.right_neighbor = (end_pos < rna.size()) ?
          rna[end_pos]->getCode() :
          IdentificationData::MoleculeParentMatch::RIGHT_TERMINUS;
        oligo.parent_matches[parent_ref].insert(match);
        id_data.registerIdentifiedOligo(oligo);
      }
    }
  }

} //namespace
