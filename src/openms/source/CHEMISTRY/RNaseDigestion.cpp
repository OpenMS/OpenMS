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
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>

using namespace std;

namespace OpenMS
{
  void RNaseDigestion::setEnzyme(const String& enzyme_name)
  {
    enzyme_ = RNaseDB::getInstance()->getEnzyme(enzyme_name);
  }

  void RNaseDigestion::digest(const String& rna, vector<String>& output,
                              Size min_length, Size max_length) const
  {
    // initialization
    output.clear();
    if (rna.empty() || ((rna.size() == 1) && (rna[0] == 'p'))) return;

    // handle terminal phosphates in input:
    bool has_5prime_p = (rna[0] == 'p');
    bool has_3prime_p = (rna[rna.size() - 1] == 'p');
    // mark sequence ends:
    String temp_rna = "^" + rna.substr(has_5prime_p, rna.size() - has_5prime_p -
                                       has_3prime_p) + "$";

    // beginning positions of "naive" fragments:
    vector<int> fragment_pos = tokenize_(temp_rna);
    // after "^" or before "$" aren't valid cleavages:
    if (fragment_pos.size() > 1)
    {
      if (fragment_pos[1] == 1)
      {
        fragment_pos.erase(fragment_pos.begin() + 1);
      }
      if (fragment_pos.back() == (int)temp_rna.size() - 1)
      {
        fragment_pos.resize(fragment_pos.size() - 1);
      }
    }

    vector<StringView> unmod_output;
    // don't apply length filters yet, because we modified the original string:
    digestAfterTokenize_(fragment_pos, temp_rna, unmod_output);

    String three_prime_gain =
      dynamic_cast<const DigestionEnzymeRNA*>(enzyme_)->getThreePrimeGain();
    String five_prime_gain =
      dynamic_cast<const DigestionEnzymeRNA*>(enzyme_)->getFivePrimeGain();

    for (vector<StringView>::iterator it = unmod_output.begin();
         it != unmod_output.end(); ++it)
    {
      String fragment = it->getString();
      Size actual_length = fragment.size();
      bool is_5prime_end = (fragment[0] == '^');
      bool is_3prime_end = (fragment[fragment.size() - 1] == '$');
      if (is_5prime_end) // original 5' end -> no 5' enzyme mod
      {
        actual_length--; // don't count the "^"
        if (has_5prime_p)
        {
          fragment[0] = 'p';
        }
        else
        {
          fragment = fragment.substr(1);
        }
      }
      else
      {
        fragment = five_prime_gain + fragment;
      }
      if (is_3prime_end) // original 3' end -> no 3' enzyme mod
      {
        actual_length--; // don't count the "$"
        if (has_3prime_p)
        {
          fragment[fragment.size() - 1] = 'p';
        }
        else
        {
          fragment = fragment.substr(0, fragment.size() - 1);
        }
      }
      else
      {
        fragment += three_prime_gain;
      }

      if (((min_length == 0) || (actual_length >= min_length)) &&
          ((max_length == 0) || (actual_length <= max_length)))
      {
        output.push_back(fragment);
      }
    }
  }


  vector<pair<Size, Size>> RNaseDigestion::getFragmentPositions_(
    const NASequence& rna, Size min_length, Size max_length,
    const boost::regex& cuts_after_regex, const boost::regex& cuts_before_regex) const
  {
    vector<pair<Size, Size>> result;
    vector<Size> fragment_pos(1, 0);
    if (!cuts_after_regex.empty() && !cuts_before_regex.empty())
    {
      for (Size i = 1; i < rna.size(); ++i)
      {
        if (boost::regex_search(rna[i - 1]->getCode(), cuts_after_regex) &&
            boost::regex_search(rna[i]->getCode(), cuts_before_regex))
        {
          fragment_pos.push_back(i);
        }
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
        if (((min_length == 0) || (length >= min_length)) &&
            ((max_length == 0) || (length <= max_length)))
        {
          result.push_back(make_pair(start_pos, length));
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

    const DigestionEnzymeRNA* rnase =
      dynamic_cast<const DigestionEnzymeRNA*>(enzyme_);

    String five_prime_code = rnase->getFivePrimeGain();
    if (five_prime_code == "p") five_prime_code = "5'-p";
    String three_prime_code = rnase->getThreePrimeGain();
    if (three_prime_code == "p") three_prime_code = "3'-p";

    static RibonucleotideDB* ribo_db = RibonucleotideDB::getInstance();
    const Ribonucleotide* five_prime_gain = five_prime_code.empty() ?
      nullptr : ribo_db->getRibonucleotide(five_prime_code);
    const Ribonucleotide* three_prime_gain = three_prime_code.empty() ?
      nullptr : ribo_db->getRibonucleotide(three_prime_code);

    boost::regex cuts_after_regex(rnase->getCutsAfterRegEx());
    boost::regex cuts_before_regex(rnase->getCutsBeforeRegEx());
    // @TODO: add special case for unspecific cleavage?

    vector<pair<Size, Size>> positions = getFragmentPositions_(
      rna, min_length, max_length, cuts_after_regex, cuts_before_regex);

    for (const auto& pos : positions)
    {
      NASequence fragment = rna.getSubsequence(pos.first, pos.second);
      if (pos.first > 0)
      {
        fragment.setFivePrimeMod(five_prime_gain);
      }
      if (pos.first + pos.second < rna.size())
      {
        fragment.setThreePrimeMod(three_prime_gain);
      }
      output.push_back(fragment);
    }
  }


  void RNaseDigestion::digest(IdentificationData& id_data, Size min_length,
                              Size max_length) const
  {
    const DigestionEnzymeRNA* rnase =
      dynamic_cast<const DigestionEnzymeRNA*>(enzyme_);

    String five_prime_code = rnase->getFivePrimeGain();
    if (five_prime_code == "p") five_prime_code = "5'-p";
    String three_prime_code = rnase->getThreePrimeGain();
    if (three_prime_code == "p") three_prime_code = "3'-p";

    static RibonucleotideDB* ribo_db = RibonucleotideDB::getInstance();
    const Ribonucleotide* five_prime_gain = five_prime_code.empty() ?
      nullptr : ribo_db->getRibonucleotide(five_prime_code);
    const Ribonucleotide* three_prime_gain = three_prime_code.empty() ?
      nullptr : ribo_db->getRibonucleotide(three_prime_code);

    boost::regex cuts_after_regex(rnase->getCutsAfterRegEx());
    boost::regex cuts_before_regex(rnase->getCutsBeforeRegEx());
    // @TODO: add special case for unspecific cleavage?

    for (IdentificationData::ParentMoleculeRef parent_ref =
           id_data.getParentMolecules().begin(); parent_ref !=
           id_data.getParentMolecules().end(); ++parent_ref)
    {
      if (parent_ref->molecule_type != IdentificationData::MoleculeType::RNA)
      {
        continue;
      }

      NASequence rna = NASequence::fromString(parent_ref->sequence);
      vector<pair<Size, Size>> positions = getFragmentPositions_(
        rna, min_length, max_length, cuts_after_regex, cuts_before_regex);

      for (const auto& pos : positions)
      {
        NASequence fragment = rna.getSubsequence(pos.first, pos.second);
        if (pos.first > 0)
        {
          fragment.setFivePrimeMod(five_prime_gain);
        }
        if (pos.first + pos.second < rna.size())
        {
          fragment.setThreePrimeMod(three_prime_gain);
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
