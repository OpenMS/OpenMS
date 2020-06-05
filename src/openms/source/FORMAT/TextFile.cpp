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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>

using namespace std;

namespace OpenMS
{

  TextFile::TextFile()
  {

  }

  TextFile::~TextFile()
  {
  }

  TextFile::TextFile(const String& filename, bool trim_lines, Int first_n, bool skip_empty_lines)
  {
    load(filename, trim_lines, first_n, skip_empty_lines);
  }

  void TextFile::load(const String& filename, bool trim_lines, Int first_n, bool skip_empty_lines)
  {
    // stream in binary mode prevents interpretation and merging of \r on Windows & MacOS
    // .. so we can deal with it ourselves in a consistent way
    ifstream is(filename.c_str(), ios_base::in | ios_base::binary);
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    buffer_.clear();

    String str;
    bool had_enough = false;
    while (getLine(is, str) && !had_enough)
    {
      if (trim_lines) str.trim();
      // skip? (only after trimming!)
      if (skip_empty_lines && str.empty()) continue;

      buffer_.push_back(str);
      if (first_n > -1 && static_cast<Int>(buffer_.size()) == first_n)
      {
        had_enough = true;
        break;
      }
    } // while
  }

  void TextFile::store(const String& filename)
  {
    ofstream os;
    // stream not opened in binary mode, thus "\n" will be evaluated platform dependent (e.g. resolve to \r\n on Windows)
    os.open(filename.c_str(), ofstream::out);

    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    for (Iterator it = buffer_.begin(); it != buffer_.end(); ++it)
    {
      if (it->hasSuffix("\n"))
      {
        if (it->hasSuffix("\r\n"))
        {
          os << it->chop(2) << "\n";
        }
        else
        {
          os << *it;
        }
      }
      else
      {
        os << *it << "\n";
      }
    }
    os.close();
  }

  std::istream& TextFile::getLine(std::istream& is, std::string& t)
  {
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    std::istream::sentry se(is, true);
    if (!se)
    { // the stream has an error
      return is;
    }
    std::streambuf* sb = is.rdbuf();

    for (;;)
    {
        int c = sb->sbumpc(); // get and advance to next char
        switch (c) {
        case '\n':
            return is;
        case '\r': // consume next '\n' (if any) and return
            if (sb->sgetc() == '\n') // peek current char
            {
              sb->sbumpc(); // consume it
            }
            return is;
        case std::streambuf::traits_type::eof():
            is.setstate(std::ios::eofbit); // still allows: while(is == true)
            if (t.empty())
            { // only if we just started a new line, we set the is.fail() == true, ie. is == false
              is.setstate(std::ios::badbit);
            }
            return is;
        default:
            t += (char)c;
        }
    }
}

  TextFile::ConstIterator TextFile::begin() const
  {
    return buffer_.begin();
  }

  TextFile::ConstIterator TextFile::end() const
  {
    return buffer_.end();
  }

  TextFile::Iterator TextFile::begin()
  {
    return buffer_.begin();
  }

  TextFile::Iterator TextFile::end()
  {
    return buffer_.end();
  }
} // namespace OpenMS
