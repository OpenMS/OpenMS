// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

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
    enum class ChecksumType
    {
      UNKNOWN_CHECKSUM, ///< Unknown checksum type
      SHA1, ///< Secure Hash Algorithm-1
      MD5, ///< Message-Digest algorithm 5
      SIZE_OF_CHECKSUMTYPE
    };

    /// Names of checksum types
    static const std::string NamesOfChecksumType[static_cast<size_t>(ChecksumType::SIZE_OF_CHECKSUMTYPE)];

    /// returns all checksum type names known to OpenMS
    static StringList getAllNamesOfChecksumType();

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
    const std::string& getNameOfFile() const;
    /// sets the file name
    void setNameOfFile(const std::string& name_of_file);

    /// returns the file path
    const std::string& getPathToFile() const;
    /// sets the file path
    void setPathToFile(const std::string& path_path_to_file);

    /// returns the file size in MB
    float getFileSize() const;
    /// sets the file size in MB
    void setFileSize(float file_size);

    /// returns the file type
    const std::string& getFileType() const;
    /// sets the file type
    void setFileType(const std::string& file_type);

    /// returns the file's checksum
    const std::string& getChecksum() const;
    /// sets the file's checksum
    void setChecksum(const std::string& checksum, ChecksumType type);
    /// returns the checksum type
    ChecksumType getChecksumType() const;

    /// Returns the native ID type of the spectra
    const std::string& getNativeIDType() const;
    /// Sets the native ID type of the spectra
    void setNativeIDType(const std::string& type);

    /// Returns the nativeID of the spectra
    const std::string& getNativeIDTypeAccession() const;
    /// Sets the native ID of the spectra
    void setNativeIDTypeAccession(const std::string& accesssion);

protected:
    std::string name_of_file_;
    std::string path_to_file_;
    double file_size_;
    std::string file_type_;
    std::string checksum_;
    ChecksumType checksum_type_ = SourceFile::ChecksumType::UNKNOWN_CHECKSUM;
    std::string native_id_type_;
    std::string native_id_type_accession_;
  };
} // namespace OpenMS

