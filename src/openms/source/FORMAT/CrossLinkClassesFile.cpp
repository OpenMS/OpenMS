// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/CrossLinkClassesFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <string>
#include <fstream>

namespace OpenMS
{
  /// Default constructor
  CrossLinkClassesFile::CrossLinkClassesFile()
  {
  }

  /// Destruktor
  CrossLinkClassesFile::~CrossLinkClassesFile()
  {
  }


  bool CrossLinkClassesFile::load(const String & filename)
  {
    int state = 0;
    String line;
    String current_class;
    std::vector< std::vector < String > > attributes;
    std::ifstream infile(filename.c_str());
    while (std::getline(infile, line))
    {
        // Ignore empty lines
        if (line.empty())
        {
            continue;
        }

        switch(state)
        {
          // Expect new name of cross-link class
          case 0:

            if (this->classes.find(line) != this->classes.end())
            {
                // TODO Allow overriding of standard cross-link classes
                LOG_ERROR << "ERROR: Name of cross-link class not unique: " << line << std::endl;
                return false;
            }
            current_class  = line;
            state = 1;
            break;

          // Expect attributes of the class
          case 1:
            if(line.hasPrefix("$"))  // End of cross-link class definition.
            {
                this->classes[current_class].push_back(attributes);
                attributes.clear();
                state = 0;
            }
            else if(line.hasPrefix(">"))  // New clause for class definition found
            {
              this->classes[current_class].push_back(attributes);
              attributes.clear();
            }
            else
            {
              StringList split;
              StringUtils::split(line,"|", split);
              size_t split_size = split.size();

              for (size_t i = 0; i < split_size; ++i)
              {
                 split[i] = split[i].removeWhitespaces();
              }
              String identifier = split[0];
              String predicate = split[1];

              // TODO Also check the keys that are tested
              if (identifier != "PEPID" && identifier != "ALPHA" && identifier != "BETA")
              {
                LOG_ERROR << "ERROR: Unknown identifier in attribute specification: '"<< identifier <<"'. Choose among 'PEPID', 'ALPHA', or 'BETA'." << std::endl;
                return false;
              }
              if(predicate != "IS" && predicate != "ISNOT" && predicate != "HAS" && predicate != "HASNOT")
              {
                 LOG_ERROR << "ERROR: Predicate is invalid. Choose among 'IS', 'ISNOT', 'HAS', or 'HASNOT'" << std::endl;
                 return false;
              }
              if ((predicate == "IS" || predicate == "ISNOT") && split_size != 4)
              {
                  LOG_ERROR << "ERROR: Predicates 'IS' and 'ISNOT' require 4 fields in total" << std::endl;
                  return false;
              }
              if ((predicate == "HAS" || predicate == "HASNOT") && split_size != 3)
              {
                  LOG_ERROR << "ERROR: Predicates 'HAS' and 'HASNOT' require 3 fields in total" << std::endl;
                  return false;
              }
              attributes.push_back(split);
            }
          default:
            break;
        }
    }
    if (state == 1)
    {
      LOG_ERROR << "ERROR: File ended, but no '$' symbol has been encountered." << std::endl;
      return false;
    }
    return true;
  }
} // namespace OpenMS
