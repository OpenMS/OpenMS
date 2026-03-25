// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <arrow/io/api.h>

namespace OpenMS
{

/**
  @brief Factory to obtain an Arrow RandomAccessFile for a named entry inside
  a zip archive (or a direct file when the archive path is a directory).

  The implementation returns a RandomAccessFile which Arrow/Parquet can read
  from directly. For ZIP archives this currently uses libzip to open the
  archive and the named entry and exposes a RandomAccessFile backed by the
  libzip entry. When the provided archive path is a directory (unpacked
  layout), a regular Arrow ReadableFile is returned for the file on disk.

  Returns an Arrow Result wrapping a shared_ptr to the RandomAccessFile.
*/
struct OPENMS_DLLAPI ZipRandomAccessFile
{
  static arrow::Result<std::shared_ptr<arrow::io::RandomAccessFile>> Open(const String& archive_path,
                                                                         const String& entry_name,
                                                                         std::unique_ptr<File::TempDir>& temp_dir);
};

} // namespace OpenMS
