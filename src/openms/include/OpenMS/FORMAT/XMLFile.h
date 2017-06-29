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

#ifndef OPENMS_FORMAT_XMLFILE_H
#define OPENMS_FORMAT_XMLFILE_H

// OpenMS includes
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

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

        @param filename The name of the file to validate.
        @param os The ostream where error messages should be send.
       
        @exception Exception::FileNotFound is thrown if the file cannot be found
        @exception Exception::NotImplemented is thrown if there is no schema available for the file type
      */
      bool isValid(const String & filename, std::ostream & os);

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
        @brief Parses the in-memory buffer given by @p buffer using the handler given by @p handler.

        @note Currently the buffer needs to be plain text, gzip buffer is not supported.

        @exception Exception::ParseError is thrown if an error occurred during the parsing
      */
      void parseBuffer_(const std::string & buffer, XMLHandler * handler);

      /**
        @brief Stores the contents of the XML handler given by @p handler in the file given by @p filename.

        @exception Exception::UnableToCreateFile is thrown if the file cannot be created
      */
      void save_(const String & filename, XMLHandler * handler) const;

      /// XML schema file location
      String schema_location_;

      /// Version string
      String schema_version_;

      /// Encoding string that replaces the encoding (system dependent or specified in the XML). Disabled if empty. Used as a workaround for XTandem output xml.
      String enforced_encoding_;

      void enforceEncoding_(const String& encoding);
    };

    /**
      @brief Encodes tabs '\t' in the string as &amp;#x9; and returns the encoded string.
      
      @param to_encode The String to encode.
      @return The encoded string.
    */
    String OPENMS_DLLAPI encodeTab(const String & to_encode);
  }   // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_XMLFILE_H

