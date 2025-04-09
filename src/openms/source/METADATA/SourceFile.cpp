// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SourceFile.h>

using namespace std;

namespace OpenMS
{
  const std::string SourceFile::NamesOfChecksumType[] = {"Unknown", "SHA-1", "MD5"};

  SourceFile::SourceFile() :
    CVTermList(),
    name_of_file_(),
    path_to_file_(),
    file_size_(),
    file_type_(),
    checksum_(),
    
    native_id_type_(""),
    native_id_type_accession_("")
  {
  }

  SourceFile::~SourceFile() = default;

  bool SourceFile::operator==(const SourceFile& rhs) const
  {
    return CVTermList::operator==(rhs) &&
           name_of_file_ == rhs.name_of_file_ &&
           path_to_file_ == rhs.path_to_file_ &&
           file_size_ == rhs.file_size_ &&
           file_type_ == rhs.file_type_ &&
           checksum_ == rhs.checksum_ &&
           checksum_type_ == rhs.checksum_type_ &&
           native_id_type_ == rhs.native_id_type_ &&
           native_id_type_accession_ == rhs.native_id_type_accession_;
  }

  bool SourceFile::operator!=(const SourceFile& rhs) const
  {
    return !(operator==(rhs));
  }

  const String& SourceFile::getNameOfFile() const
  {
    return name_of_file_;
  }

  void SourceFile::setNameOfFile(const String& name_of_file)
  {
    name_of_file_ = name_of_file;
  }

  const String& SourceFile::getPathToFile() const
  {
    return path_to_file_;
  }

  void SourceFile::setPathToFile(const String& path_to_file)
  {
    path_to_file_ = path_to_file;
  }

  float SourceFile::getFileSize() const
  {
    return file_size_;
  }

  void SourceFile::setFileSize(float file_size)
  {
    file_size_ = static_cast<double>(file_size);
  }

  const String& SourceFile::getFileType() const
  {
    return file_type_;
  }

  void SourceFile::setFileType(const String& file_type)
  {
    file_type_ = file_type;
  }

  const String& SourceFile::getChecksum() const
  {
    return checksum_;
  }

  SourceFile::ChecksumType SourceFile::getChecksumType() const
  {
    return checksum_type_;
  }

  void SourceFile::setChecksum(const String& checksum, ChecksumType type)
  {
    checksum_ = checksum;
    checksum_type_ = type;
  }

  const String& SourceFile::getNativeIDType() const
  {
    return native_id_type_;
  }

  void SourceFile::setNativeIDType(const String& type)
  {
    native_id_type_ = type;
  }

  const String& SourceFile::getNativeIDTypeAccession() const
  {
    return native_id_type_accession_;
  }

  void SourceFile::setNativeIDTypeAccession(const String& accession)
  {
    native_id_type_accession_ = accession;
  }

}

