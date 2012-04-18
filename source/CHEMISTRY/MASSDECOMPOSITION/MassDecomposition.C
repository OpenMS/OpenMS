// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>

using namespace std;

namespace OpenMS
{
  MassDecomposition::MassDecomposition() :
    number_of_max_aa_(0)
  {}

  MassDecomposition::MassDecomposition(const String & deco) :
    number_of_max_aa_(0)
  {
    String tmp = deco;
    vector<String> split;

    // some more info per line
    if (deco.has('('))
    {
      Size pos = tmp.find('(', 0);
      tmp = tmp.substr(0, pos);
      tmp.trim();
    }

    tmp.split(' ', split);
    number_of_max_aa_ = 0;
    // only one aa type?
    if (!split.empty())
    {
      Size sum = 0;

      for (Size i = 0; i != split.size(); ++i)
      {
        char aa = split[i][0];
        String s = split[i];
        s.erase(0, 1);
        Size n = (Size)s.toInt();
        if (number_of_max_aa_ < n)
        {
          number_of_max_aa_ = n;
        }
        sum += n;
        decomp_[aa] = n;
      }
    }
  }

  MassDecomposition::MassDecomposition(const MassDecomposition & rhs) :
    decomp_(rhs.decomp_),
    number_of_max_aa_(rhs.number_of_max_aa_)
  {}

  MassDecomposition & MassDecomposition::operator = (const MassDecomposition &rhs)
  {
    if (&rhs != this)
    {
      decomp_ = rhs.decomp_;
      number_of_max_aa_ = rhs.number_of_max_aa_;
    }
    return *this;
  }

  MassDecomposition & MassDecomposition::operator += (const MassDecomposition &d)
  {
    for (Map<char, Size>::const_iterator it = d.decomp_.begin(); it != d.decomp_.end(); ++it)
    {
      if (!decomp_.has(it->first))
      {
        decomp_.insert(*it);
        if (it->second > number_of_max_aa_)
        {
          number_of_max_aa_ = it->second;
        }
      }
      else
      {
        decomp_[it->first] += it->second;
        if (decomp_[it->first] > number_of_max_aa_)
        {
          number_of_max_aa_ = decomp_[it->first];
        }
      }
    }

    return *this;
  }

  bool MassDecomposition::operator < (const MassDecomposition &rhs) const
  {
    return decomp_ < rhs.decomp_;
  }

  bool MassDecomposition::operator == (const String &deco) const
  {
    MassDecomposition md(deco);

    return decomp_ == md.decomp_ && number_of_max_aa_ == md.number_of_max_aa_;
  }

  String MassDecomposition::toString() const
  {
    String s;
    for (Map<char, Size>::const_iterator it = decomp_.begin(); it != decomp_.end(); ++it)
    {
      s += it->first + String(it->second) + String(" ");
    }
    return s.trim();
  }

  String MassDecomposition::toExpandedString() const
  {
    String s;
    for (Map<char, Size>::const_iterator it = decomp_.begin(); it != decomp_.end(); ++it)
    {
      s += String(it->second, it->first);
    }
    return s;
  }

  bool MassDecomposition::containsTag(const String & tag) const
  {
    Map<char, Size> tmp;
    for (String::ConstIterator it = tag.begin(); it != tag.end(); ++it)
    {
      char aa = *it;
      if (!decomp_.has(aa))
      {
        return false;
      }
      if (tmp.has(aa))
      {
        tmp[aa]++;
      }
      else
      {
        tmp[aa] = 1;
      }
    }

    // check if tag decomp_ is compatible with decomp_
    for (Map<char, Size>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
    {
      if (decomp_.find(it->first)->second < it->second)
      {
        return false;
      }
    }

    return true;
  }

  bool MassDecomposition::compatible(const MassDecomposition & deco) const
  {
    for (Map<char, Size>::const_iterator it = deco.decomp_.begin(); it != deco.decomp_.end(); ++it)
    {
      if (!decomp_.has(it->first) || decomp_.find(it->first)->second < it->second)
      {
        cerr << it->first << " " << it->second << endl;
        return false;
      }
    }
    return true;
  }

  MassDecomposition MassDecomposition::operator + (const MassDecomposition &rhs) const
  {
    MassDecomposition d(* this);
    for (Map<char, Size>::const_iterator it = rhs.decomp_.begin(); it != rhs.decomp_.end(); ++it)
    {
      if (!d.decomp_.has(it->first))
      {
        d.decomp_.insert(*it);
        if (it->second > number_of_max_aa_)
        {
          d.number_of_max_aa_ = it->second;
        }
      }
      else
      {
        d.decomp_[it->first] += it->second;
        if (d.decomp_[it->first] > d.number_of_max_aa_)
        {
          d.number_of_max_aa_ = d.decomp_[it->first];
        }
      }
    }
    return d;
  }

  Size MassDecomposition::getNumberOfMaxAA() const
  {
    return number_of_max_aa_;
  }

} // namespace OpenMS
