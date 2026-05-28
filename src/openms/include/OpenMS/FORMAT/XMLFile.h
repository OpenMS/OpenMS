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

    /**
      @brief Base class for loading and storing XML files via Xerces, with optional schema validation and transparent gzip / bzip2 / zip decompression on load.

      Concrete file-format classes derive from @c XMLFile and provide an
      @ref XMLHandler subclass that drives the SAX parser. The base class hides the
      Xerces setup, the compression-magic detection on the input side, the suffix-based
      compression selection on the output side, and the optional encoding override used
      for non-conforming XML producers (e.g. X!Tandem).

      Compression is detected on @ref parse_ by reading the first two bytes of the file
      and matching one of three magic numbers: @c "BZ" (bzip2), @c 0x1F8B (gzip),
      @c "PK" (zip). On @ref save_, compression is selected by filename suffix
      (@c .gz, @c .bz2) instead.

      @ingroup FileIO
    */
    class OPENMS_DLLAPI XMLFile
    {

public:

      /// Construct an @ref XMLFile without schema info; @c schema_location_ remains unset, so @ref isValid cannot be used until derived-class logic initializes @c schema_location_ before calling @ref isValid.
      XMLFile();

      /**
        @brief Construct with a schema location for later @ref isValid calls.

        @param[in] schema_location Path of the XML schema (resolved at @ref isValid time via @c File::find).
        @param[in] version         Schema version string returned by @ref getVersion.
      */
      XMLFile(const String& schema_location, const String& version);

      /// Virtual destructor — defaulted; allows safe deletion through a base-class pointer.
      virtual ~XMLFile();

      /**
        @brief Check if @p filename validates against the bound XML schema.

        Resolves the configured @c schema_location_ via @c File::find and delegates to
        @c XMLValidator::isValid. Error messages are written to @p os.

        @param[in]     filename The file to validate.
        @param[in,out] os       Error-message sink.

        @return @c true if the file validates against the schema, @c false otherwise.

        @throws OpenMS::Exception::NotImplemented if no schema is bound (default-constructed instance — @c schema_location_ is empty).
        @throws OpenMS::Exception::FileNotFound   if @p filename cannot be found.
      */
      bool isValid(const String& filename, std::ostream& os);

      /// Return the schema version string passed to the parameterised constructor; empty for default-constructed instances.
      const String& getVersion() const;

protected:
      /**
        @brief Parse the XML file at @p filename through @p handler.

        Reads the first two bytes of @p filename to detect bzip2 (@c "BZ"),
        gzip (@c 0x1F8B), or zip (@c "PK") and wraps the file in a
        @ref CompressedInputSource accordingly; any other magic falls through to a
        plain @c LocalFileInputSource. If @ref enforceEncoding_ has set
        @c enforced_encoding_, the override is applied to the source before parsing.

        On exit (success or exception), @c handler->reset() is called via an RAII guard
        so the handler can be reused on subsequent calls.

        @param[in]     filename Path to the XML file.
        @param[in,out] handler  SAX handler driven by the parser.

        @throws OpenMS::Exception::FileNotFound if @p filename does not exist.
        @throws OpenMS::Exception::ParseError   if Xerces initialisation fails or a parse error occurs.
      */
      void parse_(const String& filename, XMLHandler* handler);

      /**
        @brief Parse an in-memory XML buffer through @p handler.

        Wraps @p buffer in a @c xercesc::MemBufInputSource (with the literal id
        @c "inMemory"). If @ref enforceEncoding_ has set @c enforced_encoding_, the
        override is applied to the source before parsing. The RAII handler reset described
        for @ref parse_ also applies here.

        @note The buffer must be plain text; gzip / bzip2 / zip buffers are not supported on this path (only on @ref parse_).

        @param[in]     buffer   In-memory XML text.
        @param[in,out] handler  SAX handler driven by the parser.

        @throws OpenMS::Exception::ParseError if Xerces initialisation fails or a parse error occurs.
      */
      void parseBuffer_(const std::string & buffer, XMLHandler * handler);

      /**
        @brief Stores the contents of the XML handler given by @p handler in the file given by @p filename.

        If @p filename ends with <tt>.gz</tt>, the output is written with gzip compression.
        If @p filename ends with <tt>.bz2</tt>, the output is written with bzip2 compression.
        Otherwise, uncompressed output is written.

        @param[in] filename The output filename (extension determines compression: .gz, .bz2, or none)
        @param[in] handler The XML handler containing the content to write

        @exception Exception::UnableToCreateFile is thrown if the file cannot be created
      */
      void save_(const String& filename, XMLHandler* handler) const;

      String schema_location_;        ///< Path of the XML schema for validation; empty when the default constructor was used (@ref isValid then throws @c NotImplemented).
      String schema_version_;         ///< Schema version string returned by @ref getVersion.
      String enforced_encoding_;      ///< Optional XML encoding override applied to the @c InputSource in @ref parse_ and @ref parseBuffer_; empty disables the override. Used as a workaround for XTandem output XML which carries an encoding the parser otherwise stumbles on.

      /**
        @brief Set or clear the XML-encoding override applied to subsequent @ref parse_ / @ref parseBuffer_ calls.

        Stores @p encoding in @ref enforced_encoding_. Pass an empty string to disable the
        override and let Xerces use the encoding declared inside the XML.

        @param[in] encoding Encoding name (e.g. @c "ISO-8859-1") to apply to the @c InputSource before parsing; empty string disables the override so Xerces consults the XML declaration instead.
      */
      void enforceEncoding_(const String& encoding);
    };

    /**
      @brief Encodes tabs '\\t' in the string as &amp;\#x9; and returns the encoded string.

      @param[in] to_encode The String to encode.
      @return The encoded string.
    */
    String OPENMS_DLLAPI encodeTab(const String& to_encode);
  }   // namespace Internal
} // namespace OpenMS

