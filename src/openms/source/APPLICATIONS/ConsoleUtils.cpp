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
// $Maintainer: Chris Bielow $
// $Authors:  Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/ConsoleUtils.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#ifdef OPENMS_WINDOWSPLATFORM
#include <windows.h> // for GetConsoleScreenBufferInfo()
#undef min
#undef max
#endif

#include <cstdlib>
#include <cstdio>

namespace OpenMS
{

  ConsoleUtils::ConsoleUtils() :
    console_width_(80)
  {
    // initialize the console width
    readConsoleSize_();
  }

  ConsoleUtils::ConsoleUtils(ConsoleUtils const& other) :
    console_width_(other.console_width_)
  {
  }

  void ConsoleUtils::operator=(const ConsoleUtils& other)
  {
    console_width_ = other.console_width_;
  }

  int ConsoleUtils::readConsoleSize_()
  {
    // avoid calling this function more than once
    static bool been_here = false;
    if (been_here)
      return console_width_;

    been_here = true;

    // determine column width of current console
    try
    {
      console_width_ = -1;
      char* p_env;
      p_env = getenv("COLUMNS");
      if (p_env != nullptr)
      {
        console_width_ = String(p_env).toInt();
      }
      else
      {
        LOG_DEBUG << "output shaping: COLUMNS env does not exist!" << std::endl;
#ifdef OPENMS_WINDOWSPLATFORM
        HANDLE hOut;
        CONSOLE_SCREEN_BUFFER_INFO SBInfo;
        hOut = GetStdHandle(STD_OUTPUT_HANDLE);
        GetConsoleScreenBufferInfo(hOut, &SBInfo);
        console_width_ = SBInfo.dwSize.X;
#else // Linux / MacOS
      // try "stty size" command
      // don't use QProcess, since stty will not work there
        FILE* fp = popen("stty size", "r");
        if (fp != nullptr)
        {
          char buff[100];
          if (fgets(buff, sizeof(buff), fp) != nullptr)
          {
            String output(buff);
            StringList components;
            output.split(' ', components);
            if (components.size() == 2)
              console_width_ = components[1].toInt();
          }
          else
          {
            // TODO: throw ?
            LOG_DEBUG << "Could not read 100 characters from file." << std::endl;
          }
          pclose(fp); //moved pclose outside of fgets condition - cant move it out of not nullpointer condition because pclose(null pointer is undefined behaviour)
        }
        else
        {
          LOG_DEBUG << "output shaping: stty size command failed." << std::endl;
        }
#endif
      }
      --console_width_; // to add the \n at the end of each line without forcing another line break on windows
    }
    catch (...)
    {
    }
    // if console_width_ is still -1, we do not use command line reshaping
    if (console_width_ < 10)
    {
      LOG_DEBUG << "Console width could not be determined or is smaller than 10. Not using output shaping!" << std::endl;
      console_width_ = std::numeric_limits<int>::max();
    }

    return console_width_;
  }

  String ConsoleUtils::breakString(const String& input, const Size indentation, const Size max_lines)
  {
    static ConsoleUtils instance;
    return instance.breakString_(input, indentation, max_lines);
  }

  String ConsoleUtils::breakString_(const OpenMS::String& input, const Size indentation, const Size max_lines)
  {

    // get the line length
    const Int line_len = ConsoleUtils::readConsoleSize_();

    StringList result;
    Size short_line_len = line_len - indentation;
    if (short_line_len < 1)
    {
      std::cerr << "INTERNAL ERROR: cannot split lines into empty strings! see breakString_()";
      return input;
    }
    for (Size i = 0; i < input.size(); )
    {
      String line = input.substr(i, result.size() == 0 ? line_len : short_line_len); // first line has full length
      Size advance_size = line.size();
      if (line.hasSubstring("\n"))
      {
        advance_size = 0;
        while (line.hasPrefix("\n"))
        {
          line = line.substr(1);
          ++advance_size;
        } // advance by # of \n's
        if (line.hasSubstring("\n")) line = line.prefix('\n');
        advance_size += line.size(); // + actual chars
      }

      // check if we are using the full length and split a word at the same time
      // cut a little earlier in that case for nicer looks
      if (line.size() ==  (result.size() == 0 ? line_len : short_line_len) && short_line_len > 8 && line.rfind(' ') != String::npos)
      {
        String last_word = line.suffix(' ');
        if (last_word.length() < 4)
        { // shorten by last word (will move to the next line)
          line = line.prefix(line.size() - last_word.length());
          advance_size -= last_word.size(); // + actual chars
        }
      }

      i += advance_size;
      String s_intend = (result.size() == 0 ? "" : String(indentation, ' ')); // first line no indentation
      String r = s_intend + (result.size() == 0 ? line : line.trim()); // intended lines get trimmed
      result.push_back(r); //(r.fillRight(' ', (UInt) line_len));
    }
    if (result.size() > max_lines) // remove lines from end if we get too many (but leave the last one)...
    {
      String last = result.back();
      result.erase(result.begin() + max_lines - 2, result.end());
      result.push_back((String(indentation, ' ') + String("..."))); //.fillRight(' ',(UInt) line_len));
      result.push_back(last);
    }
    // remove last " " from last line to prevent automatic linebreak
    //if (result.size()>0 && result[result.size()-1].hasSuffix(" ")) result[result.size()-1] = result[result.size()-1].substr(0,result[result.size()-1].size()-1);
    return ListUtils::concatenate(result, "\n");
  }

}
