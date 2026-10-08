// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
  MassDecomposition::MassDecomposition() :
    number_of_max_aa_(0)
  {
  }

  MassDecomposition::MassDecomposition(const std::string& deco) :
    number_of_max_aa_(0)
  {
    std::string tmp = deco;
    vector<std::string> split;

    // some more info per line
    if (StringUtils::has(deco, '('))
    {
      Size pos = tmp.find('(', 0);
      tmp = StringUtils::substr(tmp, 0, pos);
      StringUtils::trim(tmp);
    }

    StringUtils::split(tmp, ' ', split);
    number_of_max_aa_ = 0;
    // only one aa type?
    if (!split.empty())
    {
      for (Size i = 0; i != split.size(); ++i)
      {
        char aa = split[i][0];
        std::string s = split[i];
        s.erase(0, 1);
        Size n = (Size)StringUtils::toInt32(s);
        if (number_of_max_aa_ < n)
        {
          number_of_max_aa_ = n;
        }
        decomp_[aa] = n;
      }
    }
  }

  MassDecomposition::MassDecomposition(const MassDecomposition& rhs) = default;

  MassDecomposition& MassDecomposition::operator=(const MassDecomposition& rhs)
  {
    if (&rhs != this)
    {
      decomp_ = rhs.decomp_;
      number_of_max_aa_ = rhs.number_of_max_aa_;
    }
    return *this;
  }

  MassDecomposition& MassDecomposition::operator+=(const MassDecomposition& d)
  {
    for (map<char, Size>::const_iterator it = d.decomp_.begin(); it != d.decomp_.end(); ++it)
    {
      map<char, Size>::iterator it2 = decomp_.find(it->first);
      if (it2 == decomp_.end())
      {
        decomp_.insert(*it);
        if (it->second > number_of_max_aa_)
        {
          number_of_max_aa_ = it->second;
        }
      }
      else
      {
        it2->second += it->second;
        if (it2->second > number_of_max_aa_)
        {
          number_of_max_aa_ = it2->second;
        }
      }
    }

    return *this;
  }

  bool MassDecomposition::operator<(const MassDecomposition& rhs) const
  {
    return decomp_ < rhs.decomp_;
  }

  bool MassDecomposition::operator==(const std::string& deco) const
  {
    MassDecomposition md(deco);

    return decomp_ == md.decomp_ && number_of_max_aa_ == md.number_of_max_aa_;
  }

  std::string MassDecomposition::toString() const
  {
    std::string s;
    for (map<char, Size>::const_iterator it = decomp_.begin(); it != decomp_.end(); ++it)
    {
      s += it->first + StringUtils::toStr(it->second) + std::string(" ");
    }
    return StringUtils::trim(s);
  }

  std::string MassDecomposition::toExpandedString() const
  {
    std::string s;
    for (map<char, Size>::const_iterator it = decomp_.begin(); it != decomp_.end(); ++it)
    {
      s += std::string(it->second, it->first);
    }
    return s;
  }

  bool MassDecomposition::containsTag(const std::string& tag) const
  {
    map<char, Size> tmp;
    for (std::string::const_iterator it = tag.begin(); it != tag.end(); ++it)
    {
      char aa = *it;
      map<char, Size>::const_iterator it2 = decomp_.find(aa);
      if (it2 == decomp_.end())
      {
        return false;
      }

      map<char, Size>::iterator it3 = tmp.find(aa);
      if (it3 != tmp.end())
      {
        ++(it3->second);
      }
      else
      {
        tmp[aa] = 1;
      }
    }

    // check if tag decomp_ is compatible with decomp_
    for (map<char, Size>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
    {
      if (decomp_.find(it->first)->second < it->second)
      {
        return false;
      }
    }

    return true;
  }

  bool MassDecomposition::compatible(const MassDecomposition& deco) const
  {
    for (map<char, Size>::const_iterator it = deco.decomp_.begin(); it != deco.decomp_.end(); ++it)
    {
      map<char, Size>::const_iterator it2 = decomp_.find(it->first);
      if (it2 == decomp_.end() || decomp_.find(it->first)->second < it->second)
      {
        cerr << it->first << " " << it->second << endl;
        return false;
      }
    }
    return true;
  }

  MassDecomposition MassDecomposition::operator+(const MassDecomposition& rhs) const
  {
    MassDecomposition d(*this);
    for (map<char, Size>::const_iterator it = rhs.decomp_.begin(); it != rhs.decomp_.end(); ++it)
    {
      map<char, Size>::iterator it2 = d.decomp_.find(it->first);
      if (it2 == d.decomp_.end())
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
