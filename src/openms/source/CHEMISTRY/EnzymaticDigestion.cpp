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
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <boost/regex.hpp>

#include <iostream>
#include <limits>

using namespace std;

namespace OpenMS
{
  const std::string EnzymaticDigestion::NamesOfSpecificity[] = {"full", "semi", "none"};
  const std::string EnzymaticDigestion::UnspecificCleavage = "unspecific cleavage";

  EnzymaticDigestion::EnzymaticDigestion() :
    missed_cleavages_(0),
    enzyme_(ProteaseDB::getInstance()->getEnzyme("Trypsin")), // @TODO: keep trypsin as default?
    specificity_(SPEC_FULL)
  {
  }

  EnzymaticDigestion::~EnzymaticDigestion()
  {
  }

  Size EnzymaticDigestion::getMissedCleavages() const
  {
    return missed_cleavages_;
  }

  void EnzymaticDigestion::setMissedCleavages(Size missed_cleavages)
  {
    missed_cleavages_ = missed_cleavages;
  }

  void EnzymaticDigestion::setEnzyme(const DigestionEnzyme* enzyme)
  {
    enzyme_ = enzyme;
  }

  String EnzymaticDigestion::getEnzymeName() const
  {
    return enzyme_->getName();
  }

  EnzymaticDigestion::Specificity EnzymaticDigestion::getSpecificityByName(const String& name)
  {
    for (Size i = 0; i < SIZE_OF_SPECIFICITY; ++i)
    {
      if (name == NamesOfSpecificity[i]) return Specificity(i);
    }
    return SIZE_OF_SPECIFICITY;
  }

  EnzymaticDigestion::Specificity EnzymaticDigestion::getSpecificity() const
  {
    return specificity_;
  }

  void EnzymaticDigestion::setSpecificity(Specificity spec)
  {
    specificity_ = spec;
  }

  std::vector<Size> EnzymaticDigestion::tokenize_(const String& sequence) const
  {
    std::vector<Size> positions;
    Size pos = 0;
    if (enzyme_->getRegEx() != "()") // if it's not "no cleavage"
    {
      boost::regex re(enzyme_->getRegEx());
      boost::sregex_token_iterator i(sequence.begin(), sequence.end(), re, -1);
      boost::sregex_token_iterator j;
      while (i != j)
      {
        positions.push_back(pos); // 1st push 0, then all the real cleavage sites (instead of all cleavage sites and end-of-string)
        pos += i->length();
        ++i;
      }
    }
    else
    {
      positions.push_back(pos);
    }
    return positions;
  }

  bool EnzymaticDigestion::isValidProduct(const String& sequence,
                                          Size pos,
                                          Size length,
                                          bool ignore_missed_cleavages) const
  {
    if (pos >= sequence.size())
    {
      LOG_WARN << "Error: start of fragment (" << pos << ") is beyond end of sequence '" << sequence << "'!" << endl;
      return false;
    }
    if (pos + length > sequence.size())
    {
      LOG_WARN << "Error: end of fragment (" << (pos + length) << ") is beyond end of sequence '" << sequence << "'!" << endl;
      return false;
    }
    if (length == 0 || sequence.empty())
    {
      LOG_WARN << "Error: fragment and sequence must not be empty!" << endl;
      return false;
    }

    // ignore specificity and missed cleavage settings for unspecific cleavage
    if (enzyme_->getName() == UnspecificCleavage) { return true; }

    const Size end = pos + length; // past-the-end index into sequence of last fragment position

    if (specificity_ == SPEC_NONE)
    { // we don't care about terminal ends
      if (ignore_missed_cleavages) return true;
      const std::vector<Size> cleavage_positions = tokenize_(sequence); // has '0' as first site
      return (countMissedCleavages_(cleavage_positions, pos, end) <= missed_cleavages_);
    }
    else // either SPEC_SEMI or SPEC_FULL
    {
      bool spec_c = false, spec_n = false;
      const std::vector<Size> cleavage_positions = tokenize_(sequence); // has '0' as first site
      //
      // test each terminal end of the fragment
      //
      // left end (N-term for peptides):
      if (std::find(cleavage_positions.begin(), cleavage_positions.end(), pos) != cleavage_positions.end())
      { // '0' is included in cleavage_positions, so starting fragments will be found as well
        spec_n = true;
      }
      // right end (C-term for peptides):
      if (end == sequence.size())
      { // full length match (end of sequence is not in cleavage_positions)
        spec_c = true;
      }
      else if (std::find(cleavage_positions.begin(), cleavage_positions.end(), end) != cleavage_positions.end())
      {
        spec_c = true;
      }

      if ((spec_n && spec_c) || // full spec
          ((specificity_ == SPEC_SEMI) && (spec_n || spec_c))) // semi spec
      {
        if (ignore_missed_cleavages) return true;
        return (countMissedCleavages_(cleavage_positions, pos, end) <= missed_cleavages_);
      }
      return false;
    }
  }

  bool EnzymaticDigestion::filterByMissedCleavages(const String& sequence, std::function<bool(Int)> filter) const
  {
    return filter(tokenize_(sequence).size() - 1);
  }

  Size EnzymaticDigestion::countMissedCleavages_(const std::vector<Size>& cleavage_positions, Size seq_start, Size seq_end) const
  {
    Size count(0);
    for (std::vector<Size>::const_iterator it = cleavage_positions.begin(); it != cleavage_positions.end(); ++it)
    { // count MCs within fragment borders
      if ((seq_start < *it) && (*it < seq_end)) ++count;
    }
    return count;
  }

  void EnzymaticDigestion::digestAfterTokenize_(const std::vector<Size>& fragment_positions, const StringView& sequence, std::vector<StringView>& output, Size min_length, Size max_length) const
  {
    Size count = fragment_positions.size();

    // no cleavage sites? return full string
    if (count == 0)
    {
      if (sequence.size() >= min_length && sequence.size() <= max_length)
      {
        output.push_back(sequence);
      }
      return;
    }

    for (Size i = 1; i != count; ++i)
    {
      // add if cleavage product larger then min length
      Size l = fragment_positions[i] - fragment_positions[i - 1];
      if (l >= min_length && l <= max_length)
      {
        output.push_back(sequence.substr(fragment_positions[i - 1], l));
      }
    }

    // add last cleavage product (need to add because end is not a cleavage site) if larger then min length
    Size l = sequence.size() - fragment_positions[count - 1];
    if (l >= min_length && l <= max_length)
    {
      output.push_back(sequence.substr(fragment_positions[count - 1], l));
    }

    // generate fragments with missed cleavages
    for (Size i = 1; ((i <= missed_cleavages_) && (i < count)); ++i)
    {
      for (Size j = 1; j < count - i; ++j)
      {
        Size l = fragment_positions[j + i] - fragment_positions[j - 1];
        if (l >= min_length && l <= max_length)
        {
          output.push_back(sequence.substr(fragment_positions[j - 1], l));
        }
      }

      // add last cleavage product (need to add because end is not a cleavage site)
      Size l = sequence.size() - fragment_positions[count - i - 1];
      if (l >= min_length && l <= max_length)
      {
        output.push_back(sequence.substr(fragment_positions[count - i - 1], l));
      }
    }
  }

  void EnzymaticDigestion::digestUnmodified(const StringView& sequence, std::vector<StringView>& output, Size min_length, Size max_length) const
  {
    // initialization
    output.clear();

    // disable max length filter by setting to maximum length
    if (max_length == 0 || max_length > sequence.size())
    {
      max_length = sequence.size();
    }

    // Unspecific cleavage:
    // For unspecific cleavage every site is a cutting position.
    // All substrings of length min_size..max_size are generated.
    if (enzyme_->getName() == UnspecificCleavage)
    {
      output.reserve(sequence.size() * (max_length - min_length + 1));
      for (Size i = 0; i <= sequence.size() - min_length; ++i)
      {
        const Size right = std::min(i + max_length, sequence.size());
        for (Size j = i + min_length; j <= right; ++j)
        {
          output.push_back(sequence.substr(i, j - 1));
        }
      }
      return;
    }

    // naive cleavage sites
    std::vector<Size> fragment_positions = tokenize_(sequence.getString());
    digestAfterTokenize_(fragment_positions, sequence, output, min_length, max_length);
  }
} //namespace
