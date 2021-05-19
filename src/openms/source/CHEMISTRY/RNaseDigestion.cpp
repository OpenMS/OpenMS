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
  RNaseDigestion::RNaseDigestion():
    variable_inosine_(0)
  {
    setEnzyme("RNase_T1");
  }


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


  void RNaseDigestion::setVariableInosine(Size max_per_oligo)
  {
    variable_inosine_ = max_per_oligo;
  }


  Size RNaseDigestion::getVariableInosine()
  {
    return variable_inosine_;
  }


  void RNaseDigestion::digest(const NASequence& rna, vector<NASequence>& output,
                              Size min_length, Size max_length) const
  {
    output.clear();
    if (rna.empty()) return;

    vector<DigestionFragment> fragments = getFragmentPositions_(rna, min_length,
                                                                max_length);

    for (const DigestionFragment& fragment : fragments)
    {
      NASequence seq = rna.getSubsequence(fragment.start, fragment.length);
      if (fragment.start > 0)
      {
        seq.setFivePrimeMod(five_prime_gain_);
      }
      if (fragment.start + fragment.length < rna.size())
      {
        seq.setThreePrimeMod(three_prime_gain_);
      }
      output.push_back(seq);
    }
  }


  vector<RNaseDigestion::DigestionFragment> RNaseDigestion::getFragmentPositions_(
    const NASequence& rna, Size min_length, Size max_length) const
  {
    if (min_length == 0) min_length = 1;
    if ((max_length == 0) || (max_length > rna.size()))
    {
      max_length = rna.size();
    }

    vector<DigestionFragment> result;
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
      // with variable A->I, would enzyme cut after/before I, but not unmod. A?
      bool cuts_after_I = (variable_inosine_ &&
                           boost::regex_search("I", cuts_after_regex_) &&
                           !boost::regex_search("A", cuts_after_regex_));
      bool cuts_before_I = (variable_inosine_ &&
                            boost::regex_search("I", cuts_before_regex_) &&
                            !boost::regex_search("A", cuts_before_regex_));
      // list of fragment (start) positions, and whether cleavage there would
      // require a variable modification (A->I) before or after the site:
      vector<tuple<Size, bool, bool>> fragment_pos;
      fragment_pos.emplace_back(0, false, false);
      for (Size i = 1; i < rna.size(); ++i)
      {
        // standard case:
        bool cuts_after_regular =
          boost::regex_search(rna[i - 1]->getCode(), cuts_after_regex_);
        bool cuts_before_regular =
          boost::regex_search(rna[i]->getCode(), cuts_before_regex_);
        if (cuts_after_regular && cuts_before_regular)
        {
          fragment_pos.emplace_back(i, false, false);
          continue;
        }
        // cases involving variable modifications:
        bool cuts_after_special = (!cuts_after_regular && cuts_after_I &&
                                   (rna[i - 1]->getCode() == "A"));
        if (cuts_after_special && cuts_before_regular)
        {
          fragment_pos.emplace_back(i, true, false); // var. mod. required before cutting site
          continue;
        }
        bool cuts_before_special = (!cuts_before_regular && cuts_before_I &&
                                    (rna[i]->getCode() == "A"));

        if (cuts_after_regular && cuts_before_special)
        {
          fragment_pos.emplace_back(i, false, true); // var. mod. required after cutting site
          continue;
        }
        if (cuts_after_special && cuts_before_special)
        {
          fragment_pos.emplace_back(i, true, true); // var. mod. required before and after cutting site
          continue;
        }
      }
      fragment_pos.emplace_back(rna.size(), false, false); // no var. mods. required

      // "fragment_pos" has at least two elements (zero and "rna.size()"):
      for (Size start_it = 0; start_it < fragment_pos.size() - 1; ++start_it)
      {
        Size start_pos = get<0>(fragment_pos[start_it]);
        Size missed_cleavages = 0;
        for (Size end_it = start_it + 1;
             (missed_cleavages <= missed_cleavages_) &&
             (end_it < fragment_pos.size()); ++end_it)
        {
          Size end_pos = get<0>(fragment_pos[end_it]);

          Size length = end_pos - start_pos;
          if ((length >= min_length) && (length <= max_length))
          {
            DigestionFragment fragment(start_pos, length, missed_cleavages);
            fragment.req_mod_first = get<2>(fragment_pos[start_it]);
            fragment.req_mod_last = get<1>(fragment_pos[end_it]);
            if (fragment.req_mod_first + fragment.req_mod_last <=
                variable_inosine_)
            {
              result.push_back(fragment);
            }
          }
          // missed cleavage only occurs if no var. mod. is required:
          if (!get<1>(fragment_pos[end_it]) && !get<2>(fragment_pos[end_it]))
          {
            missed_cleavages++;
          }
        }
      }
    }

    return result;
  }

  void RNaseDigestion::digest(IdentificationData& id_data, Size min_length,
                              Size max_length) const
  {
    const Ribonucleotide* inosine =
      (variable_inosine_ ?
       RibonucleotideDB::getInstance()->getRibonucleotide("I") : 0);

    for (IdentificationData::ParentMoleculeRef parent_ref =
           id_data.getParentMolecules().begin(); parent_ref !=
           id_data.getParentMolecules().end(); ++parent_ref)
    {
      if (parent_ref->molecule_type != IdentificationData::MoleculeType::RNA)
      {
        continue;
      }

      NASequence rna = NASequence::fromString(parent_ref->sequence);
      vector<DigestionFragment> fragments =
        getFragmentPositions_(rna, min_length, max_length);

      for (const DigestionFragment& fragment : fragments)
      {
        NASequence seq = rna.getSubsequence(fragment.start, fragment.length);
        // if the fragment is only valid with mods, apply them:
        Size n_var_mods = 0;
        if (fragment.req_mod_first)
        {
          seq[0] = inosine;
          ++n_var_mods;
        }
        if (fragment.req_mod_last)
        {
          seq[fragment.length - 1] = inosine;
          ++n_var_mods;
        }
        if (fragment.start > 0)
        {
          seq.setFivePrimeMod(five_prime_gain_);
        }
        Size end_pos = fragment.start + fragment.length; // past-the-end pos.!
        if (end_pos < rna.size())
        {
          seq.setThreePrimeMod(three_prime_gain_);
        }
        IdentificationData::IdentifiedOligo oligo(seq);
        IdentificationData::MoleculeParentMatch match(fragment.start,
                                                      end_pos - 1);
        match.left_neighbor = (fragment.start > 0) ?
          rna[fragment.start - 1]->getCode() :
          IdentificationData::MoleculeParentMatch::LEFT_TERMINUS;
        match.right_neighbor = (end_pos < rna.size()) ?
          rna[end_pos]->getCode() :
          IdentificationData::MoleculeParentMatch::RIGHT_TERMINUS;
        oligo.parent_matches[parent_ref].insert(match);
        if (n_var_mods > 0)
        {
          oligo.setMetaValue("variable_mods", n_var_mods);
        }
        oligo.setMetaValue("missed_cleavages", fragment.missed_cleavages);
        id_data.registerIdentifiedOligo(oligo);
      }
    }
  }

} //namespace
