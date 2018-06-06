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
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <iostream>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <cstdlib>

using namespace std;

namespace OpenMS
{
  Bzip2Ifstream::Bzip2Ifstream(const char * filename) :
    n_buffer_(0), stream_at_end_(false)
  {
    file_ = fopen(filename, "rb");       //read binary: always open in binary mode because windows and mac open in text mode

    //aborting, ahhh!
    if (!file_)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    bzip2file_ = BZ2_bzReadOpen(&bzerror_, file_, 0, 0, nullptr, 0);
    if (bzerror_ != BZ_OK)
    {
      close();
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "bzip2 compression failed: ");
    }
  }

  Bzip2Ifstream::Bzip2Ifstream() :
    file_(nullptr), bzip2file_(nullptr), n_buffer_(0), bzerror_(0), stream_at_end_(true)
  {
  }

  Bzip2Ifstream::~Bzip2Ifstream()
  {
    close();
  }

  size_t Bzip2Ifstream::read(char * s, size_t n)
  {
    if (bzip2file_ != nullptr)
    {
      bzerror_ = BZ_OK;
      n_buffer_ = BZ2_bzRead(&bzerror_, bzip2file_, s, (unsigned int)n /* size of buf */);
      if (bzerror_ == BZ_OK)
      {
        return n_buffer_;
      }
      else if (bzerror_ != BZ_STREAM_END)
      {
        close();
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, " ", "bzip2 compression failed: ");
      }
      else
      {
        close();
        return n_buffer_;
      }
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "no file for decompression initialized");
    }
  }

  void Bzip2Ifstream::open(const char * filename)
  {
    close();
    file_ = fopen(filename, "rb");       //read binary: always open in binary mode because windows and mac open in text mode

    //aborting, ahhh!
    if (!file_)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    bzip2file_ = BZ2_bzReadOpen(&bzerror_, file_, 0, 0, nullptr, 0);
    if (bzerror_ != BZ_OK)
    {
      close();
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "bzip2 compression failed: ");
    }
    stream_at_end_ = false;
  }

  void Bzip2Ifstream::close()
  {
    if (bzip2file_ != nullptr)
    {
      BZ2_bzReadClose(&bzerror_, bzip2file_);
    }
    if (file_ != nullptr)
    {
      fclose(file_);
    }
    file_ = nullptr;
    bzip2file_ = nullptr;
    stream_at_end_ = true;
  }

} //namespace OpenMS
