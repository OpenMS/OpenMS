// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZipIfstream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>

#ifdef __has_include
#if __has_include(<zip.h>)
#include <zip.h>
#define OPENMS_HAVE_LIBZIP 1
#endif
#endif

#include <string>
#include <vector>

namespace OpenMS
{
  ZipIfstream::ZipIfstream() = default;

  ZipIfstream::ZipIfstream(const char* filename) :
    zip_archive_(nullptr), zip_entry_(nullptr), stream_at_end_(true)
  {
    open(filename);
  }

  ZipIfstream::~ZipIfstream()
  {
    close();
  }

  void ZipIfstream::open(const char* filename)
  {
#if defined(OPENMS_HAVE_LIBZIP)
    if (zip_entry_ != nullptr)
    {
      close();
    }

    int err = 0;
    zip_t* archive = zip_open(filename, ZIP_RDONLY, &err);
    if (archive == nullptr)
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    // Find the single non-directory entry
    zip_int64_t num_entries = zip_get_num_entries(archive, 0);
    std::vector<std::string> file_names;
    zip_int64_t file_index = -1;

    for (zip_int64_t i = 0; i < num_entries; ++i)
    {
      const char* name = zip_get_name(archive, (zip_uint64_t)i, 0);
      if (name == nullptr) continue;
      std::string sname(name);
      if (!sname.empty() && sname.back() == '/') continue; // skip directories
      file_names.push_back(sname);
      file_index = i;
    }

    if (file_names.empty())
    {
      zip_close(archive);
      throw Exception::FileEmpty(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "ZIP archive contains no file entries");
    }

    if (file_names.size() > 1)
    {
      zip_close(archive);
      std::string msg = "ZIP archive contains " + std::to_string(file_names.size()) + " file entries; expected exactly 1:";
      for (const auto& n : file_names) { msg += " " + n; }
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", msg);
    }

    zip_file_t* entry = zip_fopen_index(archive, (zip_uint64_t)file_index, 0);
    if (entry == nullptr)
    {
      zip_close(archive);
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    zip_archive_ = archive;
    zip_entry_ = entry;
    stream_at_end_ = false;
#else
    (void)filename;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "ZipIfstream requires libzip, which was not found at compile time.");
#endif
  }

  size_t ZipIfstream::read(char* s, size_t n)
  {
#if defined(OPENMS_HAVE_LIBZIP)
    if (n == 0)
    {
      return 0;
    }
    if (zip_entry_ != nullptr)
    {
      zip_int64_t bytes_read = zip_fread(static_cast<zip_file_t*>(zip_entry_), s, n);
      if (bytes_read < 0)
      {
        close();
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "error reading from zip entry", "zip_fread returned error");
      }
      if (bytes_read == 0)
      {
        close();
        stream_at_end_ = true;
      }
      return static_cast<size_t>(bytes_read);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "no file for decompression initialized");
    }
#else
    (void)s; (void)n;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "ZipIfstream requires libzip");
#endif
  }

  void ZipIfstream::close()
  {
#if defined(OPENMS_HAVE_LIBZIP)
    if (zip_entry_ != nullptr)
    {
      zip_fclose(static_cast<zip_file_t*>(zip_entry_));
    }
    if (zip_archive_ != nullptr)
    {
      zip_close(static_cast<zip_t*>(zip_archive_));
    }
#endif
    zip_entry_ = nullptr;
    zip_archive_ = nullptr;
    stream_at_end_ = true;
  }

} // namespace OpenMS
