// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <boost/mem_fn.hpp>
#include <QStringList>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace std;

namespace OpenMS
{
  
  StringList StringListUtils::fromQStringList(const QStringList & rhs)
  {
    StringList sl;
    sl.reserve(rhs.size());
    
    for (QStringList::const_iterator it = rhs.begin(); it != rhs.end(); ++it)
    {
      sl.push_back(it->toStdString());
    }
    
    return sl;
  }

  void StringListUtils::toUpper(StringList & sl)
  {
    std::for_each(sl.begin(), sl.end(), boost::mem_fn(&String::toUpper));
  }

  void StringListUtils::toLower(StringList & sl)
  {
    std::for_each(sl.begin(), sl.end(), boost::mem_fn(&String::toLower));
  }
  
  StringListUtils::Iterator StringListUtils::searchPrefix(const Iterator & start, const Iterator & end, const String & text, bool trim)
  {
    return find_if(start, end, PrefixPredicate_(text, trim));
  }
  
  StringListUtils::ConstIterator StringListUtils::searchPrefix(const ConstIterator & start, const ConstIterator & end, const String & text, bool trim)
  {
    return find_if(start, end, PrefixPredicate_(text, trim));
  }
  
  StringListUtils::ConstIterator StringListUtils::searchPrefix(const StringList & container, const String & text, bool trim)
  {
    return searchPrefix(container.begin(), container.end(), text, trim);
  }

  StringListUtils::Iterator StringListUtils::searchPrefix(StringList & container, const String & text, bool trim)
  {
    return searchPrefix(container.begin(), container.end(), text, trim);
  }
  
  StringListUtils::Iterator StringListUtils::searchSuffix(const Iterator & start, const Iterator & end, const String & text, bool trim)
  {
    return find_if(start, end, SuffixPredicate_(text, trim));
  }
  
  StringListUtils::ConstIterator StringListUtils::searchSuffix(const ConstIterator & start, const ConstIterator & end, const String & text, bool trim)
  {
    return find_if(start, end, SuffixPredicate_(text, trim));
  }

  StringListUtils::ConstIterator StringListUtils::searchSuffix(const StringList & container, const String & text, bool trim)
  {
    return searchSuffix(container.begin(), container.end(), text, trim);
  }
  
  StringListUtils::Iterator StringListUtils::searchSuffix(StringList & container, const String & text, bool trim)
  {
    return searchSuffix(container.begin(), container.end(), text, trim);
  }
  
} // namespace OpenMS
