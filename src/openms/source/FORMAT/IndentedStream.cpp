// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IndentedStream.h>

#include <OpenMS/CONCEPT/Colorizer.h>

#include <algorithm>
#include <sstream>

using namespace std;

namespace OpenMS
{
  
  IndentedStream::IndentedStream(std::ostream& stream, const UInt indentation, const UInt max_lines) :
      stream_(&stream), indentation_(indentation), max_lines_(max_lines), max_line_width_(ConsoleUtils::getInstance().getConsoleWidth())
  {
  }

  /// D'tor flushes the stream

  IndentedStream::~IndentedStream()
  {
    stream_->flush();
  }

  IndentedStream& IndentedStream::operator<<(Colorizer& colorizer)
  {
    // manipulate the internal data of colorizer (if any)
    stringstream reformatted;
    // use a clone of ourselves, but dump data to a stringstream
    IndentedStream formatter(reformatted, indentation_, max_lines_);
    formatter.current_column_pos_ = current_column_pos_; // advance the formatter to the same column position that we have
    formatter << colorizer.getInternalChars_().str(); // push the data (invoking line breaks if required)
    // update colorizer with new data (with indentation)
    colorizer.setInternalChars_(reformatted.str()); // skip the initialization 'x' chars

    // apply Color to our internal stream
    // Do NOT push the data into the IndentedStream since this prevents detection
    // of stdout/stderr (and its redirection status) by Colorizer. If the underlying stream_
    // is indeed cout but redirected to a file, you would get ANSI symbols in there (not desiable)
    *stream_ << colorizer;

    // update our column position based on the new data
    current_column_pos_ = formatter.current_column_pos_; // this does not take into account ANSI codes which are added by Colorizer -- nothing we can do about it.

    return *this;
  }

  IndentedStream& IndentedStream::operator<<(IndentedStream & self)
  {
    return self;
  }
  
  IndentedStream& IndentedStream::operator<<(StreamManipulator manip)
  {
    // call the function on the internal stream
    manip(*stream_);
    return *this;
  }

  IndentedStream& IndentedStream::indent(const UInt new_indent)
  {
    indentation_ = new_indent;
    return *this;
  }
} // namespace OpenMS
