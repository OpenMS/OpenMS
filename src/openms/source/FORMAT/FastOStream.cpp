// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: $
// $Authors: Tom Lukas Lankenau, Anton Kriese $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FastOStream.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstring>

namespace OpenMS
{
  //FastOStream::FastOStream(const std::string& path):
  //    string_constructed(true)
  //{
  //  os_ = new std::fstream( path.c_str() ); // xxx Mag er nicht
  //  //os_->open()
  //}

  //FastOStream::FastOStream(const String& path):
  //    string_constructed(true)
  //{
  //  os_ = new std::ostream; // xxx Mag er auch nicht =(
  //  os_->open(path.c_str());
  //}

  //FastOStream::~FastOStream(){
  //  if(string_constructed){
  //    delete os_;
  //  }
  //}

  FastOStream& FastOStream::operator << (const OpenMS::String& s)
  {
    os_.rdbuf()->sputn(s.c_str(), s.size());
    return *this;
  }

  FastOStream& FastOStream::operator << (const std::string& s)
  {
    os_.rdbuf()->sputn(s.c_str(), s.size());
    return *this;
  }

  FastOStream& FastOStream::operator << (const char* const s)
  {
    os_.rdbuf()->sputn(s, strlen(s));
    return *this;
  }

  FastOStream& operator << (FastOStream& os, const DataValue& p)
  {
    /// for doubles or lists of doubles, you get full precision. Use DataValue::toString(false) if you only need low precision

    switch (p.value_type_)
    {
      case DataValue::STRING_VALUE: os << *(p.data_.str_); break;
      case DataValue::STRING_LIST: os << *(p.data_.str_list_); break;
      case DataValue::INT_LIST: os << *(p.data_.int_list_); break;
      case DataValue::DOUBLE_LIST: os << *(p.data_.dou_list_); break;
      case DataValue::INT_VALUE: os << p.data_.ssize_; break; // using our String conversion (faster than std::ofstream)
      case DataValue::DOUBLE_VALUE: os << p.data_.dou_; break; // using our String conversion (faster than std::ofstream)
      case DataValue::EMPTY_VALUE: break;
    }
    return os;
  }

  std::ostream& FastOStream::getStream()
  {
    return os_;
  }
}
