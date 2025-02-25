// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <boost/regex.hpp>

#include <limits>

using namespace std;

namespace OpenMS
{
  void ProteaseDigestion::setEnzyme(const String& enzyme_name)
  {
    enzyme_ = ProteaseDB::getInstance()->getEnzyme(enzyme_name);
    re_.reset(new boost::regex(enzyme_->getRegEx()));
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
      if (i > missed_cleavages_)
      {
        break;
      }
      sum += count - i;
    }
    return sum;
  }

  Size ProteaseDigestion::digest(const AASequence& protein, vector<AASequence>& output, Size min_length, Size max_length) const
  {
    // initialization
    output.clear();
    std::vector<std::pair<size_t,size_t>> idcs; // small overhead filling intermediate vector first and iterating again
    Size wrong_size = digest(protein, idcs, min_length, max_length);
    output.reserve(idcs.size());
    std::transform(idcs.begin(), idcs.end(), std::back_inserter(output),
      [&protein](std::pair<size_t, size_t>& start_end)
      {
        return protein.getSubsequence(start_end.first, UInt(start_end.second - start_end.first));
      }
    );
    return wrong_size;
  }

  Size ProteaseDigestion::digest(const AASequence& protein, vector<std::pair<size_t,size_t>>& output, Size min_length, Size max_length) const
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
      if (l >= min_length && l <= max_length)
      {
        output.emplace_back(begin, pep_positions[i]);
      }
      else
      {
        ++wrong_size;
      }
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
          if (l >= min_length && l <= max_length)
          {
            output.emplace_back(begin, pep_positions[j + mcs]);
          }
          else
          {
            ++wrong_size;
          }
          begin = pep_positions[j];
        }
      }
    }
    return wrong_size;
  }

} //namespace OpenMS

