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

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <boost/regex.hpp>

#include <iostream>
#include <limits>

using namespace std;

namespace OpenMS
{
  void ProteaseDigestion::setEnzyme(const String& enzyme_name)
  {
    enzyme_ = ProteaseDB::getInstance()->getEnzyme(enzyme_name);
    re_ = boost::regex(enzyme_->getRegEx());
  }

  bool ProteaseDigestion::isValidProduct(const String& protein,
                                          int pos,
                                          int length,
                                          bool ignore_missed_cleavages,
                                          bool allow_nterm_protein_cleavage,
                                          bool allow_random_asp_pro_cleavage) const
  {
    return isValidProduct_(protein, pos, length, ignore_missed_cleavages, allow_nterm_protein_cleavage, allow_random_asp_pro_cleavage);
  }

  bool ProteaseDigestion::isValidProduct(const AASequence& protein,
                                         int pep_pos,
                                         int pep_length,
                                         bool ignore_missed_cleavages,
                                         bool allow_nterm_protein_cleavage,
                                         bool allow_random_asp_pro_cleavage) const
  {
    String seq = protein.toUnmodifiedString();
    return isValidProduct_(seq, pep_pos, pep_length, ignore_missed_cleavages, allow_nterm_protein_cleavage, allow_random_asp_pro_cleavage);
  }

  Size ProteaseDigestion::peptideCount(const AASequence& protein)
  {
    // For unspecific cleavage every cutting position may be skipped. Thus, we get (n + 1) \choose 2 products.
    if (enzyme_->getName() == UnspecificCleavage) 
    {
      return (protein.size() + 1) * protein.size() / 2;
    };

    std::vector<int> pep_positions = tokenize_(protein.toUnmodifiedString());
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

  Size ProteaseDigestion::digest(const AASequence& protein, vector<AASequence>& output, Size min_length, Size max_length) const
  {
    // initialization
    output.clear();

    // disable max length filter by setting to maximum length
    if (max_length == 0 || max_length > protein.size())
    {
      max_length = protein.size();
    }

    Size mc = (enzyme_->getName() == UnspecificCleavage) ? std::numeric_limits<Size>::max() : missed_cleavages_;
    Size wrong_size(0);

    // naive cleavage sites
    std::vector<int> pep_positions = tokenize_(protein.toUnmodifiedString());
    pep_positions.push_back(protein.size()); // positions now contains 0, x1, ... xn, end
    Size count = pep_positions.size();
    Size begin = pep_positions[0];
    for (Size i = 1; i < count; ++i)
    {
      Size l = pep_positions[i] - begin;
      if (l >= min_length && l <= max_length) output.push_back(protein.getSubsequence(begin, l));
      else ++wrong_size;
      begin = pep_positions[i];
    }

    // missed cleavages
    if (pep_positions.size() > 1 && mc != 0) // there is at least one cleavage site (in addition to last position)!
    {
      // generate fragments with missed cleavages
      for (Size mcs = 1; ((mcs <= mc) && (mcs < count - 1)); ++mcs)
      {
        begin = pep_positions[0];
        for (Size j = 1; j < count - mcs; ++j)
        {
          Size l = pep_positions[j + mcs] - begin;
          if (l >= min_length && l <= max_length) output.push_back(protein.getSubsequence(begin, l));
          else ++wrong_size;
          begin = pep_positions[j];
        }
      }
    }
    return wrong_size;
  }

} //namespace OpenMS

