// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

// OpenMS
#include <OpenMS/DATASTRUCTURES/String.h>
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
    void setIdentifier(const String & id);

    /// retrieve document identifier (e.g. an LSID)
    const String & getIdentifier() const;

    /// exchange content with @p from
    void swap(DocumentIdentifier & from);


    /// set the file_name_ according to absolute path of the file loaded from preferably done whilst loading
    void setLoadedFilePath(const String & file_name);

    /// get the file_name_ which is the absolute path to the file loaded from
    const String & getLoadedFilePath() const;

    /// set the file_type according to the type of the file loaded from (see FileHandler::Type) preferably done whilst loading
    void setLoadedFileType(const String & file_name);

    /// get the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded from
    const FileTypes::Type & getLoadedFileType() const;

    //@}

protected:
    /// the ID (e.g. LSID)
    String id_;
    /// the path to the loaded file
    String file_path_;
    /// the type of the loaded file
    FileTypes::Type file_type_;
  };
} // namespace OpenMS

