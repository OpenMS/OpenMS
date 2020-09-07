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
// $Maintainer: Timo Sachsenberg $
// $Authors:  Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

namespace OpenMS
{

  /**
    @brief Struct that captures all information of a command line parameter
  */
  struct OPENMS_DLLAPI ParameterInformation
  {
    /// Parameter types
    enum ParameterTypes
    {
      NONE = 0, ///< Undefined type
      STRING, ///< String parameter
      INPUT_FILE, ///< String parameter that denotes an input file
      OUTPUT_FILE, ///< String parameter that denotes an output file
      DOUBLE, ///< Floating point number parameter
      INT, ///< Integer parameter
      STRINGLIST, ///< More than one String Parameter
      INTLIST, ///< More than one Integer Parameter
      DOUBLELIST, ///< More than one String Parameter
      INPUT_FILE_LIST, ///< More than one String Parameter that denotes input files
      OUTPUT_FILE_LIST, ///< More than one String Parameter that denotes output files
      FLAG, ///< Parameter without argument
      TEXT, ///< Left aligned text, see addText_
      NEWLINE ///< An empty line, see addEmptyLine_
    };

    /// name of the parameter (internal and external)
    String name;
    /// type of the parameter
    ParameterTypes type;
    /// default value of the parameter stored as string
    DataValue default_value;
    /// description of the parameter
    String description;
    /// argument in the description
    String argument;
    /// flag that indicates if this parameter is required i.e. it must differ from the default value
    bool required;
    /// flag the indicates that the parameter is advanced (this is used for writing the INI file only)
    bool advanced;
    /// StringList for special tags
    StringList tags;

    ///@name Restrictions for different parameter types
    //@{
    std::vector<String> valid_strings;
    Int min_int;
    Int max_int;
    double min_float;
    double max_float;
    //@}

    /// Constructor that takes all members in declaration order
    ParameterInformation(const String& n, ParameterTypes t, const String& arg, const DataValue& def, const String& desc, bool req, bool adv, const StringList& tag_values = StringList());

    ParameterInformation();

    ParameterInformation& operator=(const ParameterInformation& rhs);

  };

} // namespace

