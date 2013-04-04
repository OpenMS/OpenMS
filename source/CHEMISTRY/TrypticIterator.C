// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TrypticIterator.h>
#include <OpenMS/FORMAT/FastaIterator.h>
#include <fstream>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>


namespace OpenMS
{


  typedef std::pair<String, String> FASTAEntry;

  ///Constructor to intialize algorithm
  TrypticIterator::TrypticIterator() :
    PepIterator()
  {
    b_ = 0;
    e_ = 0;
    is_at_end_ = false;

  }

  ///destructor
  TrypticIterator::~TrypticIterator()
  {

  }

  TrypticIterator::TrypticIterator(const TrypticIterator & source) :
    PepIterator(source),
    f_file_(source.f_file_),
    actual_pep_(source.actual_pep_),
    is_at_end_(source.is_at_end_),
    f_entry_(source.f_entry_),
    b_(source.b_),
    e_(source.e_)
  {
    f_iterator_ = source.f_iterator_;
  }

  FASTAEntry TrypticIterator::operator*()
  {
    if (actual_pep_ == "")
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    return FASTAEntry(f_entry_.first, actual_pep_);
  }

  PepIterator & TrypticIterator::operator++()
  {
    if (actual_pep_ == "")
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    actual_pep_ = next_();

    if (f_iterator_->isAtEnd() && !hasNext_())
    {
      is_at_end_ = true;
    }
    return *this;
  }

  PepIterator * TrypticIterator::operator++(int)
  {
    if (actual_pep_ == "")
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    PepIterator * old = new TrypticIterator(*this);
    actual_pep_ = next_();
    if (f_iterator_->isAtEnd() && !hasNext_())
    {
      is_at_end_ = true;
    }
    return old;
  }

  void TrypticIterator::setFastaFile(const String & f)
  {
    std::fstream fs;
    fs.open(f.c_str());
    if (!fs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, f);
    }
    f_file_ = f;
  }

  String TrypticIterator::getFastaFile()
  {
    return f_file_;
  }

  bool TrypticIterator::begin()
  {
    if (f_file_ == "")
    {
      throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    f_iterator_ = Factory<PepIterator>::create("FastaIterator");
    f_iterator_->setFastaFile(f_file_);
    if (!f_iterator_->begin())
    {
      return false;

    }

    f_entry_ = **f_iterator_;
    actual_pep_ = next_();

    return true;
  }

  std::string TrypticIterator::next_()
  {

    std::string seq = f_entry_.second;

    while (b_ < seq.length())
    {
      bool isInSpec = (e_ == 0 || isDigestingEnd(seq[e_ - 1], seq[e_]) || e_ >= (seq.length() - 1));
      if (isInSpec)
      {
        if (e_ < seq.length())
        {
          e_++;
          if (e_ > (b_ + 1))
          {
            if (e_ == seq.length())
            {
              return seq.substr(b_, e_ - b_);
            }
            else
            {
              return seq.substr(b_, e_ - b_ - 1);
            }
          }
        }
      }
      else
      {
        e_++;
      }
      if (e_ >= seq.length())
      {
        goToNextAA_();
      }
    }

    if (f_iterator_->isAtEnd())
    {
      return "";

    }
    ++(*f_iterator_);
    f_entry_ = **f_iterator_;
    b_ = 0;
    e_ = 0;
    return next_();
  }

  bool TrypticIterator::hasNext_()
  {
    unsigned int bold = b_;
    unsigned int eold = e_;
    std::string res = next_();
    b_ = bold;
    e_ = eold;
    if (res.length() == 0)
    {
      return false;
    }
    return true;
  }

  void TrypticIterator::goToNextAA_()
  {
    std::string seq = f_entry_.second;
    b_++;
    while ((b_ < seq.length()) && !isDigestingEnd(seq[b_ - 1], seq[b_]))
    {
      b_++;
    }
    e_ = b_;

  }

  bool TrypticIterator::isAtEnd()
  {
    return is_at_end_;
  }

  bool TrypticIterator::isDigestingEnd(char aa1, char aa2)
  {
    return (aa1 == 'K' || aa1 == 'R') && aa2 != 'P';
  }

} // namespace OpenMS
