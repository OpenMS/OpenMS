// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

// OpenMS
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FileTypes.h>

namespace OpenMS
{
  /**
     @brief Manage source document information.

     This class stored information about the source document.
     Primarily this is the document id e.g. a LSID.

     For source files additional information can be stored:
     - file name
     - file type

     @ingroup Metadata
  */
  class OPENMS_DLLAPI DocumentIdentifier
  {
public:
    /** @name Constructors and Destructors
    */
    //@{

    /// Default constructor
    DocumentIdentifier();

    /// Copy constructor
    DocumentIdentifier(const DocumentIdentifier &) = default;

    /// Move constructor
    DocumentIdentifier(DocumentIdentifier&&) = default;

    /// Destructor
    virtual ~DocumentIdentifier();

    /// Assignment operator
    DocumentIdentifier & operator=(const DocumentIdentifier &) = default;

    /// Move assignment operator
    DocumentIdentifier& operator=(DocumentIdentifier&&) & = default;

    /// Equality operator
    bool operator==(const DocumentIdentifier & rhs) const;

    //@}

    /** @name Acessors
     */
    //@{

    /// set document identifier (e.g. an LSID)
    void setIdentifier(const std::string & id);

    /// retrieve document identifier (e.g. an LSID)
    const std::string & getIdentifier() const;

    /// exchange content with @p from
    void swap(DocumentIdentifier & from);


    /// Set the loaded file's absolute local path, or its remote URI without userinfo, query or fragment.
    void setLoadedFilePath(const std::string & file_name);

    /// Get the loaded file's absolute local path or redacted remote resource URI.
    const std::string & getLoadedFilePath() const;

    /// Detect local file type from content, or remote resource type from its filename (ignoring query and fragment).
    void setLoadedFileType(const std::string & file_name);

    /// get the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded from
    const FileTypes::Type & getLoadedFileType() const;

    //@}

protected:
    /// the ID (e.g. LSID)
    std::string id_;
    /// the path to the loaded file
    std::string file_path_;
    /// the type of the loaded file
    FileTypes::Type file_type_;
  };
} // namespace OpenMS
