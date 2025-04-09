// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors:  Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
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
      STRING, ///< String parameter
      INPUT_FILE, ///< String parameter that denotes an input file
      OUTPUT_FILE, ///< String parameter that denotes an output file
      OUTPUT_PREFIX, ///< String parameter that denotes an output file prefix
      OUTPUT_DIR, ///< String parameter that denotes an output directory
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
    ParamValue default_value;
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
    StringList valid_strings;
    Int min_int;
    Int max_int;
    double min_float;
    double max_float;
    //@}

    /// Constructor that takes all members in declaration order
    ParameterInformation(const String& n, ParameterTypes t, const String& arg, const ParamValue& def, const String& desc, bool req, bool adv, const StringList& tag_values = StringList());

    ParameterInformation();

    ParameterInformation(const ParameterInformation& rhs) = default;

    ParameterInformation& operator=(const ParameterInformation& rhs);

  };

} // namespace

