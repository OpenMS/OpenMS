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
  }

  bool ProteaseDigestion::isValidProduct(const String& protein,
                                          Size pep_pos,
                                          Size pep_length,
                                          bool ignore_missed_cleavages,
                                          bool methionine_cleavage) const
  {
    if (methionine_cleavage && (pep_pos == 1) && (protein[0] == 'M'))
    {
      // check the N-terminal peptide for a C-terminal cleavage site:
      pep_pos = 0;
      pep_length++;
    }
    return EnzymaticDigestion::isValidProduct(protein, pep_pos, pep_length, ignore_missed_cleavages);
  }

  bool ProteaseDigestion::isValidProduct(const AASequence& protein,
                                         Size pep_pos,
                                         Size pep_length,
                                         bool ignore_missed_cleavages,
                                         bool methionine_cleavage) const
  {
    String seq = protein.toUnmodifiedString();
    return isValidProduct(seq, pep_pos, pep_length, ignore_missed_cleavages, methionine_cleavage);
  }

  Size ProteaseDigestion::peptideCount(const AASequence& protein)
  {
    // For unspecific cleavage every cutting position may be skipped. Thus, we get (n + 1) \choose 2 products.
    if (enzyme_->getName() == UnspecificCleavage) { return (protein.size() + 1) * protein.size() / 2; };

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

  void ProteaseDigestion::digest(const AASequence& protein, vector<AASequence>& output) const
  {
    // initialization
    output.clear();

    Size mc = (enzyme_->getName() == UnspecificCleavage) ? std::numeric_limits<Size>::max() : missed_cleavages_;

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

} //namespace
