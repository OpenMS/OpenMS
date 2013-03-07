// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XMLFILE_H
#define OPENMS_FORMAT_XMLFILE_H

// OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <xercesc/framework/XMLFormatter.hpp>

namespace OpenMS
{
  namespace Internal
  {
    class XMLHandler;

    ///Base class for loading/storing XML files that have a handler derived from XMLHandler.
    class OPENMS_DLLAPI XMLFile
    {

public:

      ///Default constructor
      XMLFile();
      /// Constructor that sets the schema location
      XMLFile(const String & schema_location, const String & version);
      ///Destructor
      virtual ~XMLFile();

      /**
          @brief Checks if a file validates against the XML schema

          Error messages are printed to the error stream, unless redirected with the attribute @p os .

      @exception Exception::FileNotFound is thrown if the file cannot be found
          @exception Exception::NotImplemented is thrown if there is no schema available for the file type
      */
      bool isValid(const String & filename, std::ostream & os = std::cerr);

      ///return the version of the schema
      const String & getVersion() const;

protected:
      /**
        @brief Parses the XML file given by @p filename using the handler given by @p handler.

        @exception Exception::FileNotFound is thrown if the file is not found
        @exception Exception::ParseError is thrown if an error occurred during the parsing
      */
      void parse_(const String & filename, XMLHandler * handler);

      /**
        @brief Stores the contents of the XML handler given by @p handler in the file given by @p filename.

        @exception Exception::UnableToCreateFile is thrown if the file cannot be created
      */
      void save_(const String & filename, XMLHandler * handler) const;

      /// XML schema file location
      String schema_location_;

      /// Version string
      String schema_version_;

      /// Encoding string that replaces the encoding (system dependend or specified in the XML). Disabled if empty. Used as a workaround for XTandem output xml.
      String enforced_encoding_;

      void enforceEncoding_(const String& encoding);
    };

    // implementation of an XMLFormatTarget
    class OPENMS_DLLAPI OpenMSXMLFormatTarget :
      public xercesc::XMLFormatTarget
    {

public:

      OpenMSXMLFormatTarget(std::string & str) :
        XMLFormatTarget(),
        str_(str)
      {
      }

      virtual void writeChars(const XMLByte * const toWrite, const XMLSize_t count, xercesc::XMLFormatter * const /*formatter*/)
      {
        str_.append((const char * const)toWrite, count);
      }

      std::string & str_;
    };

    /**
       @brief Escapes a string to be storable into an XML File

           Some characters must be escaped which are allowed in user params. E.g. > and & are not in XML and
     need to be escaped. Parsing those escaped strings is automatically done by xerces
    */
    void OPENMS_DLLAPI writeXMLEscape(const String & to_escape, std::ostream & os);

    /**
    @brief Escapes a string and returns the escaped string

        Some characters must be escaped which are allowed in user params. E.g. > and & are not in XML and
    need to be escaped. Parsing those escaped strings is automatically done by xerces
    */
    String OPENMS_DLLAPI writeXMLEscape(const String & to_escape);

  }   // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FOMAT_XMLFILE_H
