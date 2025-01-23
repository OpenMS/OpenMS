// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <xercesc/sax/ErrorHandler.hpp>

#include <iostream>

namespace OpenMS
{
  /**
    @brief Validator for XML files.

    Validates an XML file against a given schema.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI XMLValidator :
    private xercesc::ErrorHandler
  {
public:
    /// Constructor
    XMLValidator();

    /**
      @brief Returns if an XML file is valid for given a schema file

      Error messages are printed to the error stream, unless redirected with the attribute @p os .

      @param filename The file to validated.
      @param schema The filename of the schema that should be used for validation.
      @param os The stream where error messages should be send to.

      @exception Exception::FileNotFound is thrown if the file cannot be found
      @exception Exception::ParseError is thrown if the parser could not be initialized
    */
    bool isValid(const String& filename, const String& schema, std::ostream& os = std::cerr);

protected:
    /// Flag if the validated file is valid
    bool valid_;
    /// File name of validated file (for error messages)
    String filename_;
    //output stream reference (for error messages)
    std::ostream* os_;

    /// @name Implementation of Xerces ErrorHandler methods
    //@{
    void warning(const xercesc::SAXParseException& exception) override;
    void error(const xercesc::SAXParseException& exception) override;
    void fatalError(const xercesc::SAXParseException& exception) override;
    void resetErrors() override;
    //@}
  };

} // namespace OpenMS

