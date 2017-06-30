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
// $Maintainer: Timo Sachsenberg,Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/FastaIteratorIntern.h>
#include <OpenMS/FORMAT/FASTAFile.h>

namespace OpenMS
{

  typedef std::pair<String, String> FASTAEntry;

  FastaIteratorIntern::FastaIteratorIntern() :
    fasta_file_("")
  {
  }

  FastaIteratorIntern::~FastaIteratorIntern()
  {

  }

  FastaIteratorIntern::FastaIteratorIntern(const FastaIteratorIntern & source) :
    PepIterator(source),
    fasta_file_(source.fasta_file_),
    entrys_(source.entrys_),
    it_(source.it_)
  {

  }

  FASTAEntry FastaIteratorIntern::operator*()
  {
    if (fasta_file_ == "")
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    return *it_;
  }

  PepIterator & FastaIteratorIntern::operator++()
  {
    if (fasta_file_ == "")
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    ++it_;
    return *this;
  }

  PepIterator * FastaIteratorIntern::operator++(int)
  {
    if (fasta_file_ == "")
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    PepIterator * old = new FastaIteratorIntern(*this);
    ++it_;
    return old;
  }

  void FastaIteratorIntern::setFastaFile(const String & f)
  {

    FASTAFile ffile;
    std::vector<FASTAFile::FASTAEntry> entries;

    ffile.load(f, entries);
    entrys_.clear();
    entrys_.resize(entries.size(), std::make_pair("", ""));
    for (Size i = 0; i < entries.size(); ++i)
    {
      entrys_[i].first = (entries[i].identifier + " " + entries[i].description);
      entrys_[i].second = entries[i].sequence;
    }

    fasta_file_ = f;
    it_ = entrys_.begin();
  }

  String FastaIteratorIntern::getFastaFile()
  {
    return fasta_file_;
  }

  bool FastaIteratorIntern::begin()
  {
    if (fasta_file_ == "")
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    return true;
  }

  bool FastaIteratorIntern::isAtEnd()
  {
    return it_ == entrys_.end();
  }

}
