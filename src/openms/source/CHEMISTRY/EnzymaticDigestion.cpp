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
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <boost/regex.hpp>

#include <iostream>
#include <limits>

using namespace std;

namespace OpenMS
{
  const std::string EnzymaticDigestion::NamesOfSpecificity[] = {"full", "semi", "none"};


  EnzymaticDigestion::EnzymaticDigestion() :
    missed_cleavages_(0),
    enzyme_(*EnzymesDB::getInstance()->getEnzyme("Trypsin")),
    specificity_(SPEC_FULL)
  {
  }

  EnzymaticDigestion::EnzymaticDigestion(const EnzymaticDigestion& rhs) :
    missed_cleavages_(rhs.missed_cleavages_),
    enzyme_(rhs.enzyme_),
    specificity_(rhs.specificity_)
  {
  }

  /// Assignment operator
  EnzymaticDigestion& EnzymaticDigestion::operator=(const EnzymaticDigestion& rhs)
  {
    if (this != &rhs)
    {
      missed_cleavages_ = rhs.missed_cleavages_;
      enzyme_ = rhs.enzyme_;
      specificity_ = rhs.specificity_;
    }
    return *this;
  }

  Size EnzymaticDigestion::getMissedCleavages() const
  {
    return missed_cleavages_;
  }

  void EnzymaticDigestion::setMissedCleavages(Size missed_cleavages)
  {
    missed_cleavages_ = missed_cleavages;
  }

  void EnzymaticDigestion::setEnzyme(const String enzyme_name)
  {
    enzyme_ = *EnzymesDB::getInstance()->getEnzyme(enzyme_name);
  }

  String EnzymaticDigestion::getEnzymeName() const
  {
    return enzyme_.getName();
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

  std::vector<Size> EnzymaticDigestion::tokenize_(const String& protein) const
  {
    std::vector<Size> pep_positions;
    Size pos = 0;
    if (enzyme_.getRegEx() != "()") // if it's not "no cleavage"
    {
      boost::regex re(enzyme_.getRegEx());
      boost::sregex_token_iterator i(protein.begin(), protein.end(), re, -1);
      boost::sregex_token_iterator j;
      while (i != j)
      {
        pep_positions.push_back(pos); // 1st push 0, then all the real cleavage sites (instead of all cleavage sites and end-of-string)
        pos += i->length();
        ++i;
      }
    }
    else
    {
      pep_positions.push_back(pos);
    }
    return pep_positions;
  }


  bool EnzymaticDigestion::isValidProduct(const AASequence& protein,
    Size pep_pos,
    Size pep_length,
    bool methionine_cleavage,
    bool ignore_missed_cleavages) const
  {
    return isValidProduct(protein.toUnmodifiedString(), pep_pos, pep_length, methionine_cleavage, ignore_missed_cleavages);
  }

  bool EnzymaticDigestion::isValidProduct(const String& protein,
                                          Size pep_pos,
                                          Size pep_length,
                                          bool methionine_cleavage,
                                          bool ignore_missed_cleavages) const
  {
    if (pep_pos >= protein.size())
    {
      LOG_WARN << "Error: start of peptide (" << pep_pos << ") is beyond end of protein '" << protein << "'!" << endl;
      return false;
    }
    else if (pep_pos + pep_length > protein.size())
    {
      LOG_WARN << "Error: end of peptide (" << (pep_pos + pep_length) << ") is beyond end of protein '" << protein << "'!" << endl;
      return false;
    }
    else if (pep_length == 0 || protein.size() == 0)
    {
      LOG_WARN << "Error: peptide or protein must not be empty!" << endl;
      return false;
    }
    
    // ignore specificity and missed cleavage settings for unspecific cleavage
    if (enzyme_.getName() == "unspecific cleavage") { return true; }
    
    const Size pep_end = pep_pos + pep_length; // past-the-end index into protein of last peptide amino acid

    if (specificity_ == SPEC_NONE)
    { // we don't care about terminal ends
      if (ignore_missed_cleavages) { return true; }
      const std::vector<Size> cleavage_positions = tokenize_(protein); // has '0' as first site
      return (countMissedCleavages_(cleavage_positions, pep_pos, pep_end) <= missed_cleavages_);
    }
    else // either SPEC_SEMI or SPEC_FULL
    {
      bool spec_c = false, spec_n = false;
      const std::vector<Size> cleavage_positions = tokenize_(protein); // has '0' as first site
      //
      // test each terminal end of peptide
      //
      // N-term:
      if (std::find(cleavage_positions.begin(), cleavage_positions.end(), pep_pos) != cleavage_positions.end())
      { // '0' is included in cleavage_positions, so N-terminal peptides will be found as well
        spec_n = true;
      }
      else if (pep_pos == 1 && methionine_cleavage && protein[0] == 'M')
      { // if allow methionine cleavage at the protein start position
        // if there were a real cleavage site at pos '1' it would be found in the first if-block, so we're not overlooking a MC here
        spec_n = true;
      }
      
      // C-term:
      if (pep_end == protein.size())
      { // full length match (protein C-term is not in cleavage_positions)
        spec_c = true;
      } 
      else if (std::find(cleavage_positions.begin(), cleavage_positions.end(), pep_end) != cleavage_positions.end())
      {
        spec_c = true;
      }

      if ( (spec_n && spec_c) || // full spec
           ((specificity_ == SPEC_SEMI) && (spec_n || spec_c))) // semi spec
      {
        if (ignore_missed_cleavages)
        {
          return true;
        }
        return (countMissedCleavages_(cleavage_positions, pep_pos, pep_end) <= missed_cleavages_);
      }
      return false;
    }
  }

  Size EnzymaticDigestion::countMissedCleavages_(const std::vector<Size>& cleavage_positions, Size pep_start, Size pep_end) const
  {
    Size count(0);
    for (std::vector<Size>::const_iterator it = cleavage_positions.begin(); it != cleavage_positions.end(); ++it)
    { // count MCs within peptide borders
      if ((pep_start < *it) && (*it < pep_end)) ++count;
    }
    return count;
  }

  Size EnzymaticDigestion::peptideCount(const AASequence& protein)
  {
    // For unspecific cleavage every cutting position may be skipped. Thus, we get (n + 1) \choose 2 products.
    if (enzyme_.getName() == "unspecific cleavage") { return (protein.size() + 1) * (protein.size()) / 2; };
    
    std::vector<Size> pep_positions = tokenize_(protein.toUnmodifiedString());
    Size count = pep_positions.size();
    // missed cleavages
    Size sum = count;
    for (Size i = 1; i < count; ++i)
    {
      if (i > missed_cleavages_) break;
      sum += count - i;
    }
    return sum;
  }

  void EnzymaticDigestion::digest(const AASequence& protein, vector<AASequence>& output) const
  {
    // initialization
    output.clear();
    
    Size mc = (enzyme_.getName() == "unspecific cleavage") ? std::numeric_limits<Size>::max() : missed_cleavages_;
    
    // naive cleavage sites
    std::vector<Size> pep_positions = tokenize_(protein.toUnmodifiedString());
    Size count = pep_positions.size();
    Size begin = pep_positions[0];
    for (Size i = 1; i < count; ++i)
    {
      output.push_back(protein.getSubsequence(begin, pep_positions[i] - begin));
      begin = pep_positions[i];
    }
    output.push_back(protein.getSubsequence(begin, protein.size() - begin));

    // missed cleavages
    if (pep_positions.size() > 0 && mc != 0) // there is at least one cleavage site!
    {
      // generate fragments with missed cleavages
      for (Size i = 1; ((i <= mc) && (count > i)); ++i)
      {
        begin = pep_positions[0];
        for (Size j = 1; j < count - i; ++j)
        {
          output.push_back(protein.getSubsequence(begin, pep_positions[j + i] - begin));
          begin = pep_positions[j];
        }
        output.push_back(protein.getSubsequence(begin, protein.size() - begin));
      }
    }
  }

  void EnzymaticDigestion::digestUnmodifiedString(const StringView sequence, std::vector<StringView>& output, Size min_length, Size max_length) const
  {
    // initialization
    output.clear();

    // disable max length filter by setting to maximum length
    if (max_length == 0)
    {
      max_length = sequence.size();
    }

    // Unspecific cleavage:
    // For unspecific cleavages every amino acid is a cutting position.
    // And all substrings of legnth minx_size to max_size are generated.
    if (enzyme_.getName() == "unspecific cleavage")
    {
      for (Size i = 0; i <= sequence.size() - min_length; ++i)
      {
        for (Size j = i + min_length; j <= i + max_length; ++j)
        {
          if (j > sequence.size()) { continue; }
          output.push_back(sequence.substr(i, j - 1));
        }
      }
      return;
    }
    
    // naive cleavage sites
    std::vector<Size> pep_positions = tokenize_(sequence.getString());
    Size count = pep_positions.size();

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
      Size l = pep_positions[i] - pep_positions[i - 1];
      if (l >= min_length && l <= max_length)
      {
        output.push_back(sequence.substr(pep_positions[i - 1], pep_positions[i] - 1));
      }
    }

    // add last cleavage product (need to add because end is not a cleavage site) if larger then min length
    Size l = sequence.size() - pep_positions[count - 1];
    if (l >= min_length && l <= max_length)
    {
      output.push_back(sequence.substr(pep_positions[count - 1], sequence.size() - 1));
    }

    // generate fragments with missed cleavages
    for (Size i = 1; ((i <= missed_cleavages_) && (i < count)); ++i)
    {
      for (Size j = 1; j < count - i; ++j)
      {
        Size l = pep_positions[j + i] - pep_positions[j - 1];
        if (l >= min_length && l <= max_length)
        {
          output.push_back(sequence.substr(pep_positions[j - 1], pep_positions[j + i] - 1));
        }
      }

      // add last cleavage product (need to add because end is not a cleavage site)
      Size l = sequence.size() - pep_positions[count - i - 1];
      if (l >= min_length && l <= max_length)
      {
        output.push_back(sequence.substr(pep_positions[count - i - 1], sequence.size() - 1 ));
      }
    }
  }
} //namespace
