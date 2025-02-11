// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

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
      XMLFile(const String& schema_location, const String& version);
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
      bool isValid(const String& filename, std::ostream& os);

      ///return the version of the schema
      const String& getVersion() const;

protected:
      /**
        @brief Parses the XML file given by @p filename using the handler given by @p handler.

        @exception Exception::FileNotFound is thrown if the file is not found
        @exception Exception::ParseError is thrown if an error occurred during the parsing
      */
      void parse_(const String& filename, XMLHandler* handler);

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
      void save_(const String& filename, XMLHandler* handler) const;

      /// XML schema file location
      String schema_location_;

      /// Version string
      String schema_version_;

      /// Encoding string that replaces the encoding (system dependent or specified in the XML). Disabled if empty. Used as a workaround for XTandem output xml.
      String enforced_encoding_;

      void enforceEncoding_(const String& encoding);
    };

    /**
      @brief Encodes tabs '\\t' in the string as &amp;\#x9; and returns the encoded string.

      @param to_encode The String to encode.
      @return The encoded string.
    */
    String OPENMS_DLLAPI encodeTab(const String& to_encode);
  }   // namespace Internal
} // namespace OpenMS


