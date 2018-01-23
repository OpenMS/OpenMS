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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_VALIDATORS_XMLVALIDATOR_H
#define OPENMS_FORMAT_VALIDATORS_XMLVALIDATOR_H

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

#endif // OPENMS_FORMAT_VALIDATORS_XMLVALIDATOR_H
