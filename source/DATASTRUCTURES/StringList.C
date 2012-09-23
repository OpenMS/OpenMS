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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <QStringList>

using namespace std;

namespace OpenMS
{

  StringList::StringList()
  {
  }

  StringList::StringList(const StringList & rhs) :
    vector<String>(rhs)
  {
  }

  StringList::StringList(const vector<String> & rhs) :
    vector<String>(rhs)
  {
  }

  StringList::StringList(const vector<string> & rhs) :
    vector<String>(rhs.begin(), rhs.end())
  {
  }

  StringList::StringList(const QStringList & rhs)
  {
    for (int i = 0; i < rhs.size(); ++i)
    {
      this->push_back(rhs[i].toStdString());
    }
  }

  StringList & StringList::operator=(const StringList & rhs)
  {
    vector<String>::operator=(rhs);
    return *this;
  }

  StringList & StringList::operator=(const vector<String> & rhs)
  {
    vector<String>::operator=(rhs);
    return *this;
  }

  StringList & StringList::operator=(const vector<string> & rhs)
  {
    this->resize(rhs.size());
    for (Size i = 0; i < rhs.size(); ++i)
    {
      this->operator[](i) = rhs[i];
    }
    return *this;
  }

  StringList StringList::create(const String & list, const char splitter)
  {
    StringList out;
    list.split(splitter, out);
    return out;
  }

  StringList StringList::create(const char * const * list, UInt size)
  {
    StringList out;
    for (UInt i = 0; i < size; ++i)
    {
      out.push_back(list[i]);
    }
    return out;
  }

  bool StringList::contains(const String & s) const
  {
    for (Size i = 0; i < this->size(); ++i)
    {
      if (this->operator[](i) == s)
        return true;
    }
    return false;
  }

  void StringList::toUpper()
  {
    for (Size i = 0; i < this->size(); ++i)
    {
      this->operator[](i).toUpper();
    }
  }

  void StringList::toLower()
  {
    for (Size i = 0; i < this->size(); ++i)
    {
      this->operator[](i).toLower();
    }
  }

  String StringList::concatenate(const String & glue) const
  {

    if (size() == 0)
      return "";

    String output = *(begin());
    for (const_iterator it = begin() + 1; it != end(); ++it)
    {
      output += (glue + *it);
    }

    return output;
  }

  // ----------------- Output operator ----------------------

  ostream & operator<<(std::ostream & os, const StringList & p)
  {
    os << "[";
    if (p.size() > 0)
    {
      os << p[0];
    }

    for (Size i = 1; i < p.size(); ++i)
    {
      os << ", " << p[i];
    }
    os << "]";
    return os;
  }

  StringList::Iterator StringList::search(const String & text, bool trim)
  {
    return search(begin(), text, trim);
  }

  StringList::ConstIterator StringList::search(const String & text, bool trim) const
  {
    return search(begin(), text, trim);
  }

  StringList::Iterator StringList::search(const Iterator & start, const String & text, bool trim)
  {
    String pattern = text;
    if (trim)
    {
      pattern.trim();
    }

    String tmp;

    Iterator it = start;

    while (it != end())
    {
      tmp = *it;
      if (trim)
      {
        if (tmp.trim().hasPrefix(pattern))
          return it;
      }
      else
      {
        if (tmp.hasPrefix(pattern))
          return it;
      }
      ++it;
    }

    //nothing found
    return end();
  }

  StringList::ConstIterator StringList::search(const ConstIterator & start, const String & text, bool trim) const
  {
    String pattern = text;
    if (trim)
    {
      pattern.trim();
    }

    String tmp;

    ConstIterator it = start;

    while (it != end())
    {
      tmp = *it;
      if (trim)
      {
        if (tmp.trim().hasPrefix(pattern))
          return it;
      }
      else
      {
        if (tmp.hasPrefix(pattern))
          return it;
      }
      ++it;
    }

    //nothing found
    return end();
  }

  StringList::Iterator StringList::searchSuffix(const Iterator & start, const String & text, bool trim)
  {
    String pattern = text;
    if (trim)
    {
      pattern.trim();
    }

    String tmp;

    Iterator it = start;

    while (it != end())
    {
      tmp = *it;
      if (trim)
      {
        if (tmp.trim().hasSuffix(pattern))
          return it;
      }
      else
      {
        if (tmp.hasSuffix(pattern))
          return it;
      }
      ++it;
    }

    //nothing found
    return end();
  }

  StringList::ConstIterator StringList::searchSuffix(const ConstIterator & start, const String & text, bool trim) const
  {
    String pattern = text;
    if (trim)
    {
      pattern.trim();
    }

    String tmp;

    ConstIterator it = start;

    while (it != end())
    {
      tmp = *it;
      if (trim)
      {
        if (tmp.trim().hasSuffix(pattern))
          return it;
      }
      else
      {
        if (tmp.hasSuffix(pattern))
          return it;
      }
      ++it;
    }

    //nothing found
    return end();
  }

  StringList::Iterator StringList::searchSuffix(const String & text, bool trim)
  {
    return searchSuffix(begin(), text, trim);
  }

  StringList::ConstIterator StringList::searchSuffix(const String & text, bool trim) const
  {
    return searchSuffix(begin(), text, trim);
  }

} // namespace OpenMS
