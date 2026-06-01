// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/DocumentIdentifier.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <filesystem>

namespace OpenMS
{

  DocumentIdentifier::DocumentIdentifier() :
    id_(),
    file_path_(),
    file_type_(FileTypes::UNKNOWN)
  {
  }

  DocumentIdentifier::~DocumentIdentifier() = default;

  void DocumentIdentifier::setIdentifier(const String & id)
  {
    id_ = id;
  }

  const String & DocumentIdentifier::getIdentifier() const
  {
    return id_;
  }

  void DocumentIdentifier::setLoadedFilePath(const String & file_name)
  {
    // only change the path if we need to, otherwise low and upper case might be altered by Qt, making comparison in tests more tricky
    // i.e., a call to this will report unmatched strings
    // Treat paths starting with '/' as absolute on all platforms (matching Qt behavior).
    // On Windows, std::filesystem considers "/data/foo" relative (no drive letter),
    // but Qt treated it as absolute.
    bool is_relative = to_path(file_name).is_relative();
#ifdef OPENMS_WINDOWSPLATFORM
    if (!file_name.empty() && file_name[0] == '/') is_relative = false;
#endif
    if (is_relative)
    {
      file_path_ = File::absolutePath(file_name);
    }
    else
    {
      file_path_ = file_name;
    }
  }

  const String & DocumentIdentifier::getLoadedFilePath() const
  {
    return file_path_;
  }

  void DocumentIdentifier::setLoadedFileType(const String & file_name)
  {
    file_type_ = FileHandler::getTypeByContent(file_name);
  }

  const FileTypes::Type & DocumentIdentifier::getLoadedFileType() const
  {
    return file_type_;
  }

  void DocumentIdentifier::swap(DocumentIdentifier & from)
  {
    std::swap(id_, from.id_);
    std::swap(file_path_, from.file_path_);
    std::swap(file_type_, from.file_type_);
  }

  bool DocumentIdentifier::operator==(const DocumentIdentifier & rhs) const
  {
    return id_ == rhs.id_;
  }

}

