// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors:  Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/ConsoleUtils.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#ifdef OPENMS_WINDOWSPLATFORM
#include <windows.h> // for GetConsoleScreenBufferInfo()
#undef min
#undef max
#endif

#include <cstdlib>
#include <cstdio>

namespace OpenMS
{

  ConsoleUtils::ConsoleUtils()
  {
    // initialize the console width
    readConsoleSize_();
  }

  const ConsoleUtils& ConsoleUtils::getInstance()
  {
    static ConsoleUtils cu;
    return cu;
  }

  int ConsoleUtils::readConsoleSize_()
  {
    // avoid calling this function more than once
    static bool been_here = false;
    if (been_here)
    {
      return console_width_;
    }
    
    been_here = true;

    // determine column width of current console
    try
    {
      console_width_ = -1;
      char* p_env = getenv("COLUMNS");
      if (p_env)
      {
        console_width_ = String(p_env).toInt();
      }
      else
      {
        OPENMS_LOG_DEBUG << "output shaping: COLUMNS env does not exist!" << std::endl;

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
            OPENMS_LOG_DEBUG << "Could not read 100 characters from file." << std::endl;
          }
          pclose(fp); //moved pclose outside of fgets condition - cant move it out of not nullpointer condition because pclose(null pointer is undefined behaviour)
        }
        else
        {
          OPENMS_LOG_DEBUG << "output shaping: stty size command failed." << std::endl;
        }
#endif
      }
      --console_width_; // to add the \n at the end of each line without forcing another line break on windows
    }
    catch (...)
    {
    }
    // if console_width_ is still -1 or too small, we do not use command line reshaping
    if (console_width_ < 10)
    {
      OPENMS_LOG_DEBUG << "Console width could not be determined or is smaller than 10. Not using output shaping!" << std::endl;
      console_width_ = std::numeric_limits<int>::max();
    }

    return console_width_;
  }

  String ConsoleUtils::breakString(const String& input, const Size indentation, const Size max_lines, const Size first_line_prefill)
  {
    return ListUtils::concatenate(getInstance().breakString_(input, indentation, max_lines, first_line_prefill), '\n');
  }

  StringList ConsoleUtils::breakStringList(const String& input, const Size indentation, const Size max_lines, const Size first_line_prefill)
  {
    return getInstance().breakString_(input, indentation, max_lines, first_line_prefill);
  }

  StringList ConsoleUtils::breakString_(const OpenMS::String& input, const Size indentation, const Size max_lines, Size first_line_prefill) const
  {
    StringList result;
    if (input.empty())
    {
      return result;
    }
    Size short_line_len = console_width_ - indentation;
    if (short_line_len < 1)
    {
      //std::cerr << "INTERNAL ERROR: cannot split lines into empty strings! see breakString_()";
      result.push_back(input);
      return result;
    }
    if ((int)first_line_prefill > console_width_)
    { // first line is already longer than console width. Assume console did an automatic linebreak.
      first_line_prefill = first_line_prefill % console_width_;
    }

    String line; /// our current line as extracted from @p input
    for (Size i = 0; i < input.size(); /* i+=? computed below */)
    {
      // first line has full length (minus the prefilled chars)
      const Size remaining_line_chars = result.empty() ? console_width_ - first_line_prefill : short_line_len;
      // the first line does not need indentation
      const Size prefix_size_current_line = result.empty() ? 0 : indentation;
      
      line = input.substr(i, remaining_line_chars);
      
      // how many chars to advance in 'input'
      Size advance_size = line.size();

      // break by internal '\n' as well
      if (auto pos = line.find('\n', 0); pos != String::npos)
      {
        line = line.substr(0, pos); // do NOT include the '\n' (it is implicitly represented by adding a new string to 'result')
        advance_size = line.size() + 1; // skip the '\n' in the input
      }

      // check if we are using the full length and split a word at the same time
      // --> cut a little earlier in that case for nicer looks
      // e.g. "test this br" + '\n' + "eak" would become "test this " + '\n' + "break"
      if (line.size() == remaining_line_chars && short_line_len > 8 && line.rfind(' ') != String::npos)
      {
        String last_word = line.suffix(' ');
        if (last_word.length() < 4)
        { // shorten by last word (will move to the next line)
          line = line.prefix(line.size() - last_word.length());
          advance_size -= last_word.size(); // + actual chars
        }
      }

      i += advance_size;
      String s_intend = String(prefix_size_current_line, ' ');
      String r = s_intend + line;
      result.push_back(r); //(r.fillRight(' ', (UInt) line_len));
    }

    if (input.back() == '\n')
    { // last char input was a linebreak (which would put the cursor at column 0 in the next line)
      // --> but we want indentation!
      result.push_back(String(indentation, ' '));
    }

    if (result.size() > max_lines) // remove lines from end if we get too many (but leave the last one)...
    {
      String last = result.back();
      result.erase(result.begin() + max_lines - 2, result.end());
      result.push_back((String(indentation, ' ') + String("..."))); //.fillRight(' ',(UInt) line_len));
      result.push_back(last);
    }
    // remove last " " from last line to prevent automatic line break
    //if (result.size()>0 && result[result.size()-1].hasSuffix(" ")) result[result.size()-1] = result[result.size()-1].substr(0,result[result.size()-1].size()-1);
    return result;
  }

}
