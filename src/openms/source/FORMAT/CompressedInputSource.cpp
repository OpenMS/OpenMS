// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/GzipInputStream.h>
#include <OpenMS/FORMAT/Bzip2InputStream.h>
#include <OpenMS/FORMAT/ZipInputStream.h>
#include <OpenMS/FORMAT/HANDLERS/StringManager.h>
#include <OpenMS/SYSTEM/File.h>
#include <memory>

namespace OpenMS
{
  CompressedInputSource::CompressedInputSource(const std::string& file_path, const std::string& header) :
    file_path_(File::absolutePath(file_path)), head_(header)
  {
  }

  CompressedInputSource::CompressedInputSource(const char16_t* file_path, const std::string& header) :
    CompressedInputSource(Internal::StringManager::convert(file_path), header)
  {
  }

  InputStream* CompressedInputSource::makeStream() const
  {
    if (head_.compare(0, 2, "BZ") == 0)
    {
      auto stream = std::make_unique<Bzip2InputStream>(file_path_);
      return stream->getIsOpen() ? stream.release() : nullptr;
    }
    if (head_.compare(0, 2, "PK") == 0)
    {
      auto stream = std::make_unique<ZipInputStream>(file_path_);
      return stream->getIsOpen() ? stream.release() : nullptr;
    }
    auto stream = std::make_unique<GzipInputStream>(file_path_);
    return stream->getIsOpen() ? stream.release() : nullptr;
  }
}
