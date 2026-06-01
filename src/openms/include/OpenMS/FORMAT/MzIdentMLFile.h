// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch, Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief File adapter for MzIdentML files

      This file adapter exposes the internal MzIdentML processing capabilities to the library. The file
      adapter interface is kept the same as idXML file adapter for downward capability reasons.
      For now, read-in will be performed with DOM write-out with STREAM

      @note due to the limited capabilities of idXML/PeptideIdentification/ProteinIdentification not all
        MzIdentML features can be supported. Development for these structures will be discontinued, a new
        interface with appropriate structures will be provided.
      @note If a critical error occurs due to the missing functionality, Exception::NotImplemented is thrown.

      @note All PSM will be read into PeptideIdentification, even the passThreshold=false, threshold will be
        read into ProteinIdentification (i.e. one id run), considered at writing also will only be the
        threshold set in ProteinIdentification
      @note All PSM will be read into PeptideIdentification, even the passThreshold=false

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzIdentMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    MzIdentMLFile();
    ///Destructor
    ~MzIdentMLFile() override;

    /**
        @brief Loads the identifications from a MzIdentML file.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, std::vector<ProteinIdentification>& poid, PeptideIdentificationList& peid);

    /**
        @brief Stores the identifications in a MzIdentML file.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const std::vector<ProteinIdentification>& poid, const PeptideIdentificationList& peid) const;

    /**
        @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

        @param[in] filename File name of the file to be checked.
        @param[out] errors Errors during the validation are returned in this output parameter.
        @param[out] warnings Warnings during the validation are returned in this output parameter.

        @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings);

    /**
        @brief Checks if a file validates against the mzIdentML schema that matches its declared version.

        Unlike the base-class @ref Internal::XMLFile::isValid (which always validates against the schema
        the file adapter was constructed with), this overload first detects the mzIdentML version declared
        in @p filename and validates against the corresponding schema (1.1.0, 1.2.0 or 1.3.0). This keeps
        older, still-valid mzIdentML files from being reported as invalid just because OpenMS now writes
        1.3.0. The schema version actually used is returned via @p used_version.

        @param[in]     filename     File name of the file to be checked.
        @param[in,out] os           Error-message sink for validation errors.
        @param[out]    used_version Receives the schema version the file was validated against.

        @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isValid(const String& filename, std::ostream& os, String& used_version);

    /**
        @brief Detect the mzIdentML version declared in @p filename (e.g. "1.1.0").

        Reads the @c version attribute / target namespace from the document header. Falls back to the
        adapter's default version (1.3.0) if no version can be determined or no matching schema exists.
    */
    String detectVersion(const String& filename) const;

private:

  };

} // namespace OpenMS

