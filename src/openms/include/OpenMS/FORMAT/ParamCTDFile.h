// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Ruben Grünberg $
// $Authors: Ruben Grünberg $
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS
{

  /**
   @brief A struct to pass information about the tool as one parameter
   */
  struct ToolInfo
  {
    std::string version_;
    std::string name_;
    std::string docurl_;
    std::string category_;
    std::string description_;
    std::vector<std::string> citations_;
  };

  /**
  @brief Serializes a Param class in paramCTD file format.
         Note: only storing is currently possible

*/
  class OPENMS_DLLAPI ParamCTDFile
  {
  public:
    ParamCTDFile() = default; ///Constructor
    
    ~ParamCTDFile() = default; ///Destructor

    /**
       @brief Write CTD file

       @param filename The name of the file the param data structure should be stored in.
       @param param The param data structure that should be stored.
       @param ToolInfo Additional information about the Tool for which the param data should be stored.

       @exception std::ios::failure is thrown if the file could not be created
     */
    void store(const std::string& filename, const Param& param, const ToolInfo& tool_info) const;

    /**
       @brief Write CTD to output stream.

       @param os_ptr The stream to which the param data should be written.
       @param param The param data structure that should be writte to stream.
       @param tool_info Additional information about the Tool for which the param data should be written.
     */
    void writeCTDToStream(std::ostream* os_ptr, const Param& param, const ToolInfo& tool_info) const;

  private:
    /**
      @brief Escapes certain characters in a string that are not allowed in XML
             Escaped characters are: & < > " '

      @param to_escape The string in which the characters should be escaped

      @returns The escaped string
     */
    static std::string escapeXML(const std::string& to_escape);

    /**
      @brief Replace all occurrences of a character in a string with a string

      @param replace_in The string in which the characters should be replaced.
      @param to_replace The character that should be replaced.
      @param replace_with The string the character should be replaced with.
     */
    static void replace(std::string& replace_in, char to_replace, const std::string& replace_with);

    const std::string schema_location_ = "/SCHEMAS/Param_1_7_0.xsd";
    const std::string schema_version_ = "1.7.0";
  };
}
