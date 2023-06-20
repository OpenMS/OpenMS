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
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/ParamCTDFile.h>

namespace OpenMS
{

  /**
  @brief Load from JSON (in a Common Workflow Language (CWL) compatible way) into the Param class.

        The JSON file must contain one top level mapping of param value names to actual values.
        These values can be one of the following types:
            - null
            - boolean
            - int
            - long
            - float
            - double
            - string
            - a CWL style file path ({ "class": "File", "path": "./myFolder/myFile.txt" })
            - an array of these

        param value names match the command line option without the leading '-'. Optionally the ':'
        can be replaced with a double underscore "__".
@code
{
    "in": {
        "class": "File",
        "path": "./myFolder/myFile.txt"
    },
    "out_prefix": "test_cwl_",
    "algorithm:threshold": 5,
    "algorithm:score_type": "ID"
}
@endcode

Same file with "__" instead of ':' as the section separator.
@code
{
    "in": {
        "class": "File",
        "path": "./myFolder/myFile.txt"
    },
    "out_prefix": "test_cwl_",
    "algorithm__threshold": 5,
    "algorithm__score_type": "ID"
}
@endcode
*/
  class OPENMS_DLLAPI ParamCWLFile
  {
  public:
    ParamCWLFile() = default; ///< Constructor

    ~ParamCWLFile() = default; ///< Destructor

    /**
      @brief Read JSON file that is formatted in CWL conforming style.

      @param filename The file from where to read the Param object.
      @param param A param object with pre-filled defaults, which are updated by the values in the JSON file
      @return returns true if file was successfully loaded.

      @exception Exception::FileNotFound is thrown if the file could not be found
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    bool load(const std::string& filename, Param& param);
  };
} // namespace OpenMS
