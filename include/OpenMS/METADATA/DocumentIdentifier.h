// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_DOCUMENTIDENTIFIER_H
#define OPENMS_METADATA_DOCUMENTIDENTIFIER_H

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
    /// default constructor
    DocumentIdentifier();

    /// Copy constructor
    DocumentIdentifier(const DocumentIdentifier & source);

    /// Assignment operator
    DocumentIdentifier & operator=(const DocumentIdentifier & source);

    /// Equality operator
    bool operator==(const DocumentIdentifier & rhs) const;

    /// destructor
    virtual ~DocumentIdentifier();
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


    /// set the file_name_ according to absolute path of the file loaded from preferrably done whilst loading
    void setLoadedFilePath(const String & file_name);

    /// get the file_name_ which is the absolute path to the file loaded from
    const String & getLoadedFilePath() const;

    /// set the file_type according to the type of the file loaded from (see FileHandler::Type) preferrably done whilst loading
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

#endif // OPENMS_METADATA_DOCUMENTIDENTIFIER_H
