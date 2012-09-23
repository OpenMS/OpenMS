// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>

using namespace std;

namespace OpenMS
{
  MassDecomposition::MassDecomposition() :
    number_of_max_aa_(0)
  {
  }

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
  {
  }

  MassDecomposition & MassDecomposition::operator=(const MassDecomposition & rhs)
  {
    if (&rhs != this)
    {
      decomp_ = rhs.decomp_;
      number_of_max_aa_ = rhs.number_of_max_aa_;
    }
    return *this;
  }

  MassDecomposition & MassDecomposition::operator+=(const MassDecomposition & d)
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

  bool MassDecomposition::operator<(const MassDecomposition & rhs) const
  {
    return decomp_ < rhs.decomp_;
  }

  bool MassDecomposition::operator==(const String & deco) const
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

  MassDecomposition MassDecomposition::operator+(const MassDecomposition & rhs) const
  {
    MassDecomposition d(*this);
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
