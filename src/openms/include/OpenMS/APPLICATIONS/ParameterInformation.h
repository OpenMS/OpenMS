// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors:  Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/ParamValue.h>

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
      STRING, ///< std::string parameter
      INPUT_FILE, ///< std::string parameter that denotes an input file
      OUTPUT_FILE, ///< std::string parameter that denotes an output file
      OUTPUT_PREFIX, ///< std::string parameter that denotes an output file prefix
      OUTPUT_DIR, ///< std::string parameter that denotes an output directory
      DOUBLE, ///< Floating point number parameter
      INT, ///< Integer parameter
      STRINGLIST, ///< More than one std::string Parameter
      INTLIST, ///< More than one Integer Parameter
      DOUBLELIST, ///< More than one std::string Parameter
      INPUT_FILE_LIST, ///< More than one std::string Parameter that denotes input files
      OUTPUT_FILE_LIST, ///< More than one std::string Parameter that denotes output files
      FLAG, ///< Parameter without argument
      TEXT, ///< Left aligned text, see addText_
      NEWLINE ///< An empty line, see addEmptyLine_
    };

    /// name of the parameter (internal and external)
    std::string name;
    /// type of the parameter
    ParameterTypes type;
    /// default value of the parameter stored as string
    ParamValue default_value;
    /// description of the parameter
    std::string description;
    /// argument in the description
    std::string argument;
    /// flag that indicates if this parameter is required i.e. it must differ from the default value
    bool required;
    /// flag the indicates that the parameter is advanced (this is used for writing the INI file only)
    bool advanced;
    /// StringList for special tags
    StringList tags;

    ///@name Restrictions for different parameter types
    //@{
    StringList valid_strings;
    Int min_int;
    Int max_int;
    double min_float;
    double max_float;
    //@}

    /// Constructor that takes all members in declaration order
    ParameterInformation(const std::string& n, ParameterTypes t, const std::string& arg, const ParamValue& def, const std::string& desc, bool req, bool adv, const StringList& tag_values = StringList());

    ParameterInformation();

    ParameterInformation(const ParameterInformation& rhs) = default;

    ParameterInformation& operator=(const ParameterInformation& rhs);

  };

} // namespace

