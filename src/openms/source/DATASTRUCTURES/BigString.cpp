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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/BigString.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/config.h>

#include <iostream>

namespace OpenMS
{

  typedef std::pair<String, String> FASTAEntry;

  BigString::BigString() :
    big_string_("$"),
    separator_('$'),
    count_(1),
    len_(1)
  {
    sep_indices_.push_back(0);
    FASTA_header_.push_back("");
  }

  BigString::BigString(const BigString & bs) :
    big_string_(bs.big_string_),
    separator_(bs.separator_),
    count_(bs.count_),
    len_(bs.len_),
    sep_indices_(bs.sep_indices_),
    FASTA_header_(bs.FASTA_header_)
  {

  }

  BigString::~BigString()
  {

  }

  void BigString::add(FASTAEntry const & new_entry)
  {
    big_string_ += new_entry.second;
    big_string_ += separator_;
    ++count_;
    len_ += new_entry.second.length() + 1;
    sep_indices_.push_back(len_ - 1);
    FASTA_header_.push_back(new_entry.first);
  }

  void BigString::setSeparator(const char sep)
  {
    separator_ = sep;
  }

  char BigString::getSeparator()
  {
    return separator_;
  }

  Size BigString::size()
  {
    return count_;
  }

  Size BigString::length()
  {
    return len_;
  }

  void BigString::getPeptide(FASTAEntry & entry, Size start, Size length)
  {
    Size index_start = getIndex_(start);
    if (index_start != getIndex_(start + length))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "desired peptide is part of 2 fasta entries", "");
    }
    entry.first = FASTA_header_[index_start];
    entry.second = big_string_.substr(start, length);
    return;
  }

  const String & BigString::getBigString() const
  {
    return big_string_;
  }

  Size BigString::getIndex_(Size index, Size start, Size end)
  {
    if (end - start <= 1)
    {
      if (sep_indices_[start] >= index)
      {
        return start;
      }
      else
      {
        return start + 1;
      }
    }
    Size half = (Size) ((end - start) / 2) + start;

    if (index > sep_indices_[half])
    {
      return getIndex_(index, half, end);
    }
    else if (index < sep_indices_[half])
    {
      return getIndex_(index, start, half);
    }
    else
    {
      return half;
    }
  }

  Size BigString::getIndex_(Size index)
  {
    return getIndex_(index, 0, sep_indices_.size());
  }

}
