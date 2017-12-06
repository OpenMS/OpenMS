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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <limits>
#include <string>

using namespace std;

namespace OpenMS
{

  SVOutStream::SVOutStream(const String& file_out,
                           const String& sep,
                           const String& replacement,
                           String::QuotingMethod quoting)
    :
    ostream(nullptr), ofs_(nullptr), sep_(sep), replacement_(replacement), nan_("nan"),
    inf_("inf"), quoting_(quoting), modify_strings_(true), newline_(true)
  {
    ofs_ = new std::ofstream;
    ofs_->open(file_out.c_str());

    if (!ofs_->is_open())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_out);
    }

    // bind to filestream
    this->rdbuf(ofs_->rdbuf());

    // use high decimal precision (appropriate for double):
    precision(std::numeric_limits<double>::digits10);
  }

  SVOutStream::SVOutStream(ostream& out, const String& sep,
                           const String& replacement,
                           String::QuotingMethod quoting) :
    ostream(out.rdbuf()), ofs_(nullptr), sep_(sep), replacement_(replacement), nan_("nan"),
    inf_("inf"), quoting_(quoting), modify_strings_(true), newline_(true)
  {
    // use high decimal precision (appropriate for double):
    precision(std::numeric_limits<double>::digits10);
  }


  SVOutStream::~SVOutStream()
  {
    if (ofs_)
    {
      ofs_->close();
      delete ofs_;
    }
  }

  SVOutStream& SVOutStream::operator<<(String str)
  {
    if (str.find('\n') != String::npos)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "argument must not contain newline characters");
    }

    if (!newline_)
    {
      (ostream&) *this << sep_;
    }
    else
    {
      newline_ = false;
    }

    if (!modify_strings_)
    {
      (ostream&) *this << str;
    }
    else if (quoting_ != String::NONE)
    {
      (ostream&) *this << str.quote('"', quoting_);
    }
    else
    {
      (ostream&) *this << str.substitute(sep_, replacement_);
    }
    return *this;
  }

  SVOutStream& SVOutStream::operator<<(const char* c_str)
  {
    return operator<<(String(c_str));
  }

  SVOutStream& SVOutStream::operator<<(const char c)
  {
    return operator<<(String(c));
  }

  SVOutStream& SVOutStream::operator<<(const std::string& str)
  {
    return operator<<((String&)str);
  }

  SVOutStream& SVOutStream::operator<<(ostream& (*fp)(ostream&))
  {
    // check for "std::endl":
    // this doesn't work in LLVM/clang's libc++ (used on Mac OS X 10.9):
    // ostream& (*const endlPointer)(ostream&) = &endl;
    // if (fp == endlPointer) newline_ = true;
    fp(ss_);
    if (ss_.str() == "\n")
    {
      newline_ = true;
      ss_.str("");
    }
    (ostream&) *this << fp;
    return *this;
  }

  SVOutStream& SVOutStream::operator<<(enum Newline)
  {
    newline_ = true;
    (ostream&) *this << "\n";
    return *this;
  }

  SVOutStream& SVOutStream::write(const String& str)
  {
    ostream::write(str.c_str(), str.size());
    return *this;
  }

  bool SVOutStream::modifyStrings(bool modify)
  {
    bool old = modify_strings_;
    modify_strings_ = modify;
    return old;
  }

} // namespace OpenMS
