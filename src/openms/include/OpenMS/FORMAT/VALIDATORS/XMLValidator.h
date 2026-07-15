// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <iostream>

namespace OpenMS
{
  /**
    @brief Validator for XML files.

    Validates an XML file against a given schema.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI XMLValidator
  {
public:
    /// Constructor
    XMLValidator();

    /**
      @brief Returns if an XML file is valid for given a schema file

      Error messages are printed to the error stream, unless redirected with the attribute @p os .

      @param[in] filename The file to validated.
      @param[in] schema The filename of the schema that should be used for validation.
      @param[in] os The stream where error messages should be send to.

      @exception Exception::FileNotFound is thrown if the file cannot be found
      @exception Exception::ParseError is thrown if the parser could not be initialized
    */
    bool isValid(const std::string& filename, const std::string& schema, std::ostream& os = std::cerr);

protected:
    /// Flag if the validated file is valid
    bool valid_;
    /// File name of validated file (for error messages)
    std::string filename_;
    //output stream reference (for error messages)
    std::ostream* os_;

    /// Record a validation problem. Called by the internal Xerces ErrorHandler
    /// bridge (defined in the .cpp); keeps Xerces out of this public header.
    void logError_(const std::string& kind, const std::string& message, unsigned long line, unsigned long column);

    /// The internal Xerces ErrorHandler bridge needs access to @c logError_.
    friend class XMLValidatorErrorHandler_;
  };

} // namespace OpenMS
