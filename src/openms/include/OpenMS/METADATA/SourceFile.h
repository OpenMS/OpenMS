// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/CVTermList.h>

namespace OpenMS
{
  /**
      @brief Description of a file location, used to store the origin of (meta) data.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI SourceFile :
    public CVTermList
  {
public:
    ///Type of the checksum
    enum ChecksumType
    {
      UNKNOWN_CHECKSUM, ///< Unknown checksum type
      SHA1, ///< Secure Hash Algorithm-1
      MD5, ///< Message-Digest algorithm 5
      SIZE_OF_CHECKSUMTYPE
    };

    /// Names of checksum types
    static const std::string NamesOfChecksumType[SIZE_OF_CHECKSUMTYPE];

    /// Constructor
    SourceFile();
    /// Copy constructor
    SourceFile(const SourceFile&) = default;
    /// Move constructor
    SourceFile(SourceFile&&) = default;
    /// Destructor
    ~SourceFile() override;

    /// Assignment operator
    SourceFile& operator=(const SourceFile&) = default;
    /// Move assignment operator
    SourceFile& operator=(SourceFile&&) = default;

    /// Equality operator
    bool operator==(const SourceFile& rhs) const;
    /// Equality operator
    bool operator!=(const SourceFile& rhs) const;

    /// returns the file name
    const String& getNameOfFile() const;
    /// sets the file name
    void setNameOfFile(const String& name_of_file);

    /// returns the file path
    const String& getPathToFile() const;
    /// sets the file path
    void setPathToFile(const String& path_path_to_file);

    /// returns the file size in MB
    float getFileSize() const;
    /// sets the file size in MB
    void setFileSize(float file_size);

    /// returns the file type
    const String& getFileType() const;
    /// sets the file type
    void setFileType(const String& file_type);

    /// returns the file's checksum
    const String& getChecksum() const;
    /// sets the file's checksum
    void setChecksum(const String& checksum, ChecksumType type);
    /// returns the checksum type
    ChecksumType getChecksumType() const;

    /// Returns the native ID type of the spectra
    const String& getNativeIDType() const;
    /// Sets the native ID type of the spectra
    void setNativeIDType(const String& type);

    /// Returns the nativeID of the spectra
    const String& getNativeIDTypeAccession() const;
    /// Sets the native ID of the spectra
    void setNativeIDTypeAccession(const String& accesssion);

protected:
    String name_of_file_;
    String path_to_file_;
    double file_size_;
    String file_type_;
    String checksum_;
    ChecksumType checksum_type_ = SourceFile::ChecksumType::UNKNOWN_CHECKSUM;
    String native_id_type_;
    String native_id_type_accession_;
  };
} // namespace OpenMS

