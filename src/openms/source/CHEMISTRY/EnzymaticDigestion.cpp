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
// $Authors: Marc Sturm, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{
  const std::string EnzymaticDigestion::NamesOfSpecificity[] = {"full", "semi", "none"};
  const std::string EnzymaticDigestion::NoCleavage = "no cleavage";
  const std::string EnzymaticDigestion::UnspecificCleavage = "unspecific cleavage";

  EnzymaticDigestion::EnzymaticDigestion() :
    missed_cleavages_(0),
    enzyme_(ProteaseDB::getInstance()->getEnzyme("Trypsin")), // @TODO: keep trypsin as default?
    re_(enzyme_->getRegEx()),
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
    re_ = boost::regex(enzyme_->getRegEx());
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

  std::vector<int> EnzymaticDigestion::tokenize_(const String& sequence, int start, int end) const
  {
    std::vector<int> positions;
    // set proper boundaries
    start = std::max(0, start);
    if (end < 0 || end > (int)sequence.size()) end = (int)sequence.size();

    if (enzyme_->getRegEx() != "()") // if it's not "no cleavage"
    {
      boost::sregex_token_iterator i(sequence.begin() + start, sequence.begin() + end, re_, -1);
      boost::sregex_token_iterator j;
      while (i != j)
      {
        positions.push_back(start); // first push 'start' (usually 0), then all the real cleavage sites
        start += (int)i->length();
        ++i;
      }
    }
    else
    {
      positions.push_back(start);
    }
    return positions;
  }

  bool EnzymaticDigestion::isValidProduct(const String& sequence,
                                          int pos,
                                          int length,
    bool ignore_missed_cleavages) const
  {
    return isValidProduct_(sequence, pos, length, ignore_missed_cleavages, false, false);
  }

  bool EnzymaticDigestion::filterByMissedCleavages(const String& sequence, std::function<bool(Int)> filter) const
  {
    return filter(Int(tokenize_(sequence).size() - 1));
  }

  bool EnzymaticDigestion::isValidProduct_(const String& sequence,
                                           int pos,
                                           int length,
                                           bool ignore_missed_cleavages,
                                           bool allow_nterm_protein_cleavage,
                                           bool allow_random_asp_pro_cleavage) const
    {
    // for XTandem specific rules (see https://github.com/OpenMS/OpenMS/issues/2497)
    // M or MX at the N-terminus might have been cleaved.
    if (allow_nterm_protein_cleavage && (pos <= 2) && (sequence[0] == 'M'))
    {
      // check the N-terminal peptide for a C-terminal cleavage site:
      length += pos;
      pos = 0;
    }

    if (pos >= (int)sequence.size())
    {
      OPENMS_LOG_WARN << "Error: start of fragment (" << pos << ") is beyond end of sequence '" << sequence << "'!" << endl;
      return false;
    }
    if (pos + length > (int)sequence.size())
    {
      OPENMS_LOG_WARN << "Error: end of fragment (" << (pos + length) << ") is beyond end of sequence '" << sequence << "'!" << endl;
      return false;
    }
    if (length == 0 || sequence.empty())
    {
      OPENMS_LOG_WARN << "Error: fragment and sequence must not be empty!" << endl;
      return false;
    }

    // ignore specificity and missed cleavage settings for unspecific cleavage
    if (enzyme_->getName() == UnspecificCleavage) { return true; }

    const int end = pos + length; // past-the-end index into sequence of last fragment position

    if (specificity_ == SPEC_NONE)
    { // we don't care about terminal ends
      if (ignore_missed_cleavages) return true;
      // tokenize_ is really slow, so reduce work by working on substring:
      const std::vector<int> cleavage_positions = tokenize_(sequence, pos, end); // has 'pos' as first site
      return (cleavage_positions.size() - 1) <= missed_cleavages_;
    }
    else // either SPEC_SEMI or SPEC_FULL
    {
      bool spec_c = false, spec_n = false;
      // tokenize_ is really slow, so reduce work by working on substring with +-2 chars margin:
      const std::vector<int> cleavage_positions = tokenize_(sequence, pos - 2, end + 2); // has max(0,pos-2) as first site 
      
      //
      // test each terminal end of the fragment
      //
      // left end (N-term for peptides):
      if (std::find(cleavage_positions.begin(), cleavage_positions.end(), pos) != cleavage_positions.end())
      { // '0' is included in cleavage_positions, so starting fragments will be found as well
        spec_n = true;
      }
      // pos is > 0 at this point, so [pos-1] is valid
      else if (allow_random_asp_pro_cleavage && (sequence[pos - 1] == 'D') && (sequence[pos] == 'P'))
      {
        spec_n = true;
      }

      // right end (C-term for peptides):
      if (end == (int)sequence.size())
      { // full length match (end of sequence is not in cleavage_positions)
        spec_c = true;
      }
      else if (std::find(cleavage_positions.rbegin(), cleavage_positions.rend(), end) != cleavage_positions.rend())
      { // use rbegin() since we expect this to be the correct hit
        spec_c = true;
      }
      else if (allow_random_asp_pro_cleavage && (sequence[end - 1] == 'D') && (sequence[end] == 'P'))
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

  Size EnzymaticDigestion::countMissedCleavages_(const std::vector<int>& cleavage_positions, Size seq_start, Size seq_end) const
  {
    Size count(0);
    for (auto it = cleavage_positions.begin(); it != cleavage_positions.end(); ++it)
    { // count MCs within fragment borders
      if (((int)seq_start < *it) && (*it < (int)seq_end)) ++count;
    }
    return count;
  }

  Size EnzymaticDigestion::digestAfterTokenize_(const std::vector<int>& fragment_positions, const StringView& sequence, std::vector<std::pair<Size,Size>>& output, Size min_length, Size max_length) const
  {
    Size count = fragment_positions.size();
    Size wrong_size(0);
    Size l(0); //length

    // no cleavage sites? return full string
    if (count == 0)
    {
      if (sequence.size() >= min_length && sequence.size() <= max_length)
      {
        output.emplace_back(0, sequence.size() - 1);
      }
      return wrong_size;
    }

    for (Size i = 1; i != count; ++i)
    {
      // add if cleavage product larger than min length
      l = fragment_positions[i] - fragment_positions[i - 1];
      if (l >= min_length && l <= max_length)
      {
        output.emplace_back(fragment_positions[i - 1], l);
      }
      else ++wrong_size;
    }

    // add last cleavage product (need to add because end is not a cleavage site) if larger than min length
    l = sequence.size() - fragment_positions[count - 1];
    if (l >= min_length && l <= max_length)
    {
      output.emplace_back(fragment_positions[count - 1], l);
    }
    else ++wrong_size;

    // generate fragments with missed cleavages
    for (Size i = 1; ((i <= missed_cleavages_) && (i < count)); ++i)
    {
      for (Size j = 1; j < count - i; ++j)
      {
        l = fragment_positions[j + i] - fragment_positions[j - 1];
        if (l >= min_length && l <= max_length)
        {
          output.emplace_back(fragment_positions[j - 1], l);
        }
        else ++wrong_size;
      }

      // add last cleavage product (need to add because end is not a cleavage site)
      l = sequence.size() - fragment_positions[count - i - 1];
      if (l >= min_length && l <= max_length)
      {
        output.emplace_back(fragment_positions[count - i - 1], l);
      }
      else ++wrong_size;
    }
    return wrong_size;
  }

  Size EnzymaticDigestion::digestAfterTokenize_(const std::vector<int>& fragment_positions, const StringView& sequence, std::vector<StringView>& output, Size min_length, Size max_length) const
  {
    Size count = fragment_positions.size();
    Size wrong_size(0);

    // no cleavage sites? return full string
    if (count == 0)
    {
      if (sequence.size() >= min_length && sequence.size() <= max_length)
      {
        output.push_back(sequence);
      }
      return wrong_size;
    }

    for (Size i = 1; i != count; ++i)
    {
      // add if cleavage product larger than min length
      Size l = fragment_positions[i] - fragment_positions[i - 1];
      if (l >= min_length && l <= max_length)
      {
        output.push_back(sequence.substr(fragment_positions[i - 1], l));
      }
      else ++wrong_size;
    }

    // add last cleavage product (need to add because end is not a cleavage site) if larger than min length
    Size l = sequence.size() - fragment_positions[count - 1];
    if (l >= min_length && l <= max_length)
    {
      output.push_back(sequence.substr(fragment_positions[count - 1], l));
    }
    else ++wrong_size;

    // generate fragments with missed cleavages
    for (Size i = 1; ((i <= missed_cleavages_) && (i < count)); ++i)
    {
      for (Size j = 1; j < count - i; ++j)
      {
        Size m = fragment_positions[j + i] - fragment_positions[j - 1];
        if (m >= min_length && m <= max_length)
        {
          output.push_back(sequence.substr(fragment_positions[j - 1], m));
        }
        else ++wrong_size;
      }

      // add last cleavage product (need to add because end is not a cleavage site)
      Size n = sequence.size() - fragment_positions[count - i - 1];
      if (n >= min_length && n <= max_length)
      {
        output.push_back(sequence.substr(fragment_positions[count - i - 1], n));
      }
      else ++wrong_size;
    }
    return wrong_size;
  }

  Size EnzymaticDigestion::digestUnmodified(const StringView& sequence, std::vector<StringView>& output, Size min_length, Size max_length) const
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
          output.emplace_back(sequence.substr(i, j - i));
        }
      }
      return 0;
    }

    // naive cleavage sites
    std::vector<int> fragment_positions = tokenize_(sequence.getString());
    return digestAfterTokenize_(fragment_positions, sequence, output, min_length, max_length);
  }

  Size EnzymaticDigestion::digestUnmodified(const StringView& sequence, std::vector<std::pair<Size,Size>>& output, Size min_length, Size max_length) const
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
          output.emplace_back(i, j - i);
        }
      }
      return 0;
    }

    // naive cleavage sites
    std::vector<int> fragment_positions = tokenize_(sequence.getString());
    return digestAfterTokenize_(fragment_positions, sequence, output, min_length, max_length);
  }

} //namespace
