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

#pragma once

#include <vector>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <ostream>

namespace OpenMS
{
  class FastOStream
  {
public:
    // Rule of six
    FastOStream() = delete;
    FastOStream(const FastOStream& rhs) = delete;
    FastOStream& operator=(const FastOStream& rhs) = delete;
    FastOStream(FastOStream&& rhs) = delete;
    FastOStream& operator=(FastOStream&& rhs) = delete;
    ~FastOStream() = default;

    FastOStream(std::ostream& os) : os_(os) {};

    FastOStream& operator << (const String& s);

    FastOStream& operator << (const std::string& s);

    FastOStream& operator << (const char* const s);

    void write(const String& s);

    void write(const std::string& s);

    void write(const char* const s);

    template <typename T>
    FastOStream& operator << (const T& s)
    {
      if constexpr (std::is_arithmetic<typename std::decay<T>::type>::value)
      {
        buffer_.clear();
        buffer_ += s; // use internal buffer, instead of String(s), which would allocate
        os_.rdbuf()->sputn(buffer_.c_str(), buffer_.size());
      }
      else
      {
        os_ << s;
      }
      return *this;
    }

    template <typename T>
    inline FastOStream& operator << (const std::vector<T>& v)
    {
      *this << "[";

      if (!v.empty())
      {
        auto end = v.end() - 1;
        for (auto it = v.begin(); it < end; ++it)
        { // convert T to String using FastOStream::operator<< (very fast!)
          *this << *it << ", ";
        }
        *this << v.back();
      }

      *this << "]";
      return *this;
    }

    std::ostream& getStream();

private:
    std::ostream& os_;
    String buffer_;
  };
}
