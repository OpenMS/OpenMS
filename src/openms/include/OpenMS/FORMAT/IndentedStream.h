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

#include <OpenMS/APPLICATIONS/ConsoleUtils.h>

#include <sstream>

using namespace std;

namespace OpenMS
{
  class Colorizer;

  /**
    @brief Class for writing data which spans multiple lines with an indentation for each line (all except the first).

    Internally, ConsoleUtils is used to determine the width of the current console.
    
    The stream that is written to can be any ostream (including stdout or stderr).

    If a single item passed to IndentedStream::operator<< spans multiple indented lines (e.g. a large string),
    at most @p max_lines will be retained (excess lines will be replaced by '...').

    You can manually insert extra linebreaks by pushing '\n' into the stream (they can be part of a larger string).

    The class supports coloring its output if the underlying @p stream is either std::cout or cerr by passing a Colorizer.

    Upon destruction of IndentedStream, the underlying @p stream is flushed.
  */
  class OPENMS_DLLAPI IndentedStream
  {
  public:
    /**
      @brief C'tor
      
      @param stream Underlying stream to write to (its lifetime must exceed the one of this IndentedStream)
      @param indentation Number of spaces in front of each new line written to @p stream
      @param max_lines Shorten excessive single items to at most this many number of lines (replacing excess with '...')
    */
    IndentedStream(std::ostream& stream, const UInt indentation, const UInt max_lines);

    /// D'tor flushes the stream
    ~IndentedStream();

    /// Support normal usage of Colorizer (for coloring cout/cerr). The underlying stream will receive ANSI codes unless its a redirected(!) cout/cerr.
    /// Warning: the ANSI codes are NOT considered to advance the cursor and will lead to broken formatting if
    ///          the underlying @p stream is NOT cout/cerr.
    ///          I.e. using an IndentedStream with an underlying stringstream in combination with a Colorizer will mess up the formatting.
    IndentedStream& operator<<(Colorizer& colorizer);
    
    /// Support calling our member functions within a stream
    IndentedStream& operator<<(IndentedStream& self);

    template<typename T>
    IndentedStream& operator<<(const T& data)
    {
      // convert data to string
      std::stringstream str_data;
      str_data << data;

      auto result = ConsoleUtils::breakStringList(str_data.str(), indentation_, max_lines_, current_column_pos_);
      if (result.empty())
      {
        return *this;
      }
      
      if (result.size() == 1)
      { // no new linebreak. advance our position
        current_column_pos_ += result.back().size();
      }
      else
      { // new line: this is our new position
        current_column_pos_ = result.back().size();
      }
      

      // push result into stream
      *stream_ << result[0];
      for (size_t i = 1; i < result.size(); ++i)
      {
        *stream_ << '\n';
        *stream_ << result[i];
      }      

      return *this;
    }
        
    /// Function pointer to a function that takes an ostream, and returns it, e.g. std::endl
    typedef std::ostream& (*StreamManipulator)(std::ostream&);

    /// Overload for function pointers, e.g. std::endl
    IndentedStream& operator<<(StreamManipulator manip);

    /// Support new indentation, on the fly.
    /// This will take effect when the next line break is encountered (either manual or automatic linebreak (at the right side of the console).
    IndentedStream& indent(const UInt new_indent);

  private:
    std::ostream* stream_;        ///< the underlying stream to print to
    UInt indentation_;            ///< number of spaces in prefix of each line
    UInt max_lines_;              ///< maximum number of lines a single item is split into before excess lines are replaced by '...'
    UInt max_line_width_;         ///< width of console/output
    Size current_column_pos_ = 0; ///< length of last(=current) line
  };

} // namespace OpenMS
