// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZipArchiveFile.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#ifdef __has_include
#if __has_include(<zip.h>)
#include <zip.h>
#define OPENMS_HAVE_LIBZIP 1
#endif
#endif

#include <filesystem>
#include <fstream>
#include <vector>

namespace OpenMS
{

void ZipArchiveFile::zipDirectory(const String& directory_path, const String& output_zip)
{
#if defined(OPENMS_HAVE_LIBZIP)
  const std::filesystem::path dirpath = std::filesystem::u8path(std::string(directory_path));
  const std::filesystem::path outpath = std::filesystem::u8path(std::string(output_zip));
  const String output_zip_abs = File::absolutePath(output_zip);
  if (File::exists(output_zip_abs))
  {
    File::remove(output_zip_abs);
  }

  int err = 0;
  zip_t* za = zip_open(outpath.string().c_str(), ZIP_CREATE | ZIP_TRUNCATE, &err);
  if (!za)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to create zip archive", output_zip_abs);
  }

  for (auto it = std::filesystem::recursive_directory_iterator(dirpath); it != std::filesystem::recursive_directory_iterator(); ++it)
  {
    if (it->is_directory()) continue;
    const auto full = it->path();
    // compute relative path inside the zip
    std::string rel = std::filesystem::relative(full, dirpath).generic_string();
    zip_source_t* zs = zip_source_file(za, full.string().c_str(), 0, 0);
    if (!zs)
    {
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to add file to zip", full.string());
    }
    zip_int64_t idx = zip_file_add(za, rel.c_str(), zs, ZIP_FL_ENC_UTF_8);
    if (idx < 0)
    {
      zip_source_free(zs);
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to add file to zip", full.string());
    }
    // Create an uncompressed archive
    if (zip_set_file_compression(za, static_cast<zip_uint64_t>(idx), ZIP_CM_STORE, 0) < 0)
    {
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to set stored (no compression) for file in zip", full.string());
    }
  }

  if (zip_close(za) < 0)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to finalize zip archive", output_zip_abs);
  }
#else
  (void)directory_path;
  (void)output_zip;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}

String ZipArchiveFile::unzipDirectory(const String& input_path, std::unique_ptr<File::TempDir>& temp_dir)
{
#if defined(OPENMS_HAVE_LIBZIP)
  if (File::isDirectory(input_path))
  {
    return input_path;
  }

  if (!File::readable(input_path))
  {
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_path);
  }

  temp_dir = std::make_unique<File::TempDir>();
  const String unpack_dir = temp_dir->getPath() + "/parquet_unpacked";
  File::makeDir(unpack_dir);

  int err = 0;
  zip_t* za = zip_open(input_path.c_str(), ZIP_RDONLY, &err);
  if (!za)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to open zip archive", input_path);
  }

  zip_int64_t num = zip_get_num_entries(za, 0);
  const std::filesystem::path base_path = std::filesystem::u8path(std::string(unpack_dir)).lexically_normal();
  for (zip_uint64_t i = 0; i < static_cast<zip_uint64_t>(num); ++i)
  {
    const char* name = zip_get_name(za, i, 0);
    if (!name) continue;
    std::string entry_name(name);

    // Protect against path traversal: do not allow absolute paths or '..' to escape the base unpack dir.
    std::filesystem::path entry_path = std::filesystem::u8path(entry_name);
    if (entry_path.is_absolute())
    {
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Zip entry has absolute path", entry_name);
    }

    std::filesystem::path outpath = (base_path / entry_path).lexically_normal();
    const std::string base_str = base_path.string();
    const std::string out_str = outpath.string();
    if (out_str.size() < base_str.size() || out_str.compare(0, base_str.size(), base_str) != 0)
    {
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Zip entry would extract outside target directory", entry_name);
    }

    // create parent dirs or directory entries
    if (!entry_name.empty() && entry_name.back() == '/')
    {
      std::filesystem::create_directories(outpath);
      continue;
    }
    if (outpath.has_parent_path()) std::filesystem::create_directories(outpath.parent_path());

    zip_file_t* zf = zip_fopen_index(za, i, 0);
    if (!zf)
    {
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to read zip entry", entry_name);
    }

    std::ofstream ofs(outpath.string(), std::ios::binary);
    if (!ofs.is_open())
    {
      zip_fclose(zf);
      zip_close(za);
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, outpath.string());
    }

    constexpr size_t BUF_SIZE = 1 << 16;
    std::vector<char> buffer(BUF_SIZE);
    zip_int64_t nread = 0;
    while ((nread = zip_fread(zf, buffer.data(), buffer.size())) > 0)
    {
      ofs.write(buffer.data(), static_cast<std::streamsize>(nread));
      if (ofs.fail())
      {
        zip_fclose(zf);
        zip_close(za);
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, outpath.string());
      }
    }
    if (nread < 0)
    {
      zip_fclose(zf);
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to read zip entry", entry_name);
    }

    ofs.close();
    zip_fclose(zf);
  }

  zip_close(za);

  return unpack_dir;
#else
  (void)input_path;
  (void)temp_dir;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}



void ZipArchiveFile::addOrReplaceFromFile(const String& archive_path, const String& entry_name, const String& source_file_path)
{
#if defined(OPENMS_HAVE_LIBZIP)
  int err = 0;
  // Use ZIP_CREATE (omit ZIP_RDONLY) to open in read/write mode where supported.
  zip_t* za = zip_open(archive_path.c_str(), ZIP_CREATE, &err);
  if (!za)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to open zip archive for add/replace", archive_path);
  }

  zip_source_t* src = zip_source_file(za, source_file_path.c_str(), 0, -1);
  if (!src)
  {
    zip_close(za);
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to create zip source from file", source_file_path);
  }

  zip_int64_t idx = zip_name_locate(za, entry_name.c_str(), 0);
  if (idx >= 0)
  {
    if (zip_replace(za, idx, src) < 0)
    {
      zip_source_free(src);
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "zip_replace failed: " + String(zip_strerror(za)), "");
    }
  }
  else
  {
    if (zip_file_add(za, entry_name.c_str(), src, ZIP_FL_OVERWRITE | ZIP_FL_ENC_UTF_8) < 0)
    {
      zip_source_free(src);
      zip_close(za);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "zip_file_add failed: " + String(zip_strerror(za)), "");
    }
  }

  // Ensure stored (no compression) for Parquet files
  zip_int64_t new_idx = zip_name_locate(za, entry_name.c_str(), 0);
  if (new_idx >= 0)
  {
    zip_set_file_compression(za, static_cast<zip_uint64_t>(new_idx), ZIP_CM_STORE, 0);
  }

  if (zip_close(za) < 0)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to finalize zip archive", archive_path);
  }
#else
  (void)archive_path;
  (void)entry_name;
  (void)source_file_path;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}

std::vector<String> ZipArchiveFile::listEntries(const String& archive_path)
{
  std::vector<String> out;
#if defined(OPENMS_HAVE_LIBZIP)
  int err = 0;
  zip_t* za = zip_open(archive_path.c_str(), ZIP_RDONLY, &err);
  if (!za) return out;
  zip_int64_t num = zip_get_num_entries(za, 0);
  for (zip_uint64_t i = 0; i < static_cast<zip_uint64_t>(num); ++i)
  {
    const char* name = zip_get_name(za, i, 0);
    if (name) out.emplace_back(name);
  }
  zip_close(za);
#else
  (void)archive_path;
#endif
  return out;
}

void ZipArchiveFile::writeSidecarIndex(const String& archive_path)
{
#if defined(OPENMS_HAVE_LIBZIP)
  // If path is a directory, scan files under it
  if (File::isDirectory(archive_path))
  {
    std::vector<std::pair<std::string, uint64_t>> entries;
    const std::filesystem::path base = std::filesystem::u8path(std::string(archive_path));
    for (auto it = std::filesystem::recursive_directory_iterator(base); it != std::filesystem::recursive_directory_iterator(); ++it)
    {
      if (it->is_directory()) continue;
      const auto full = it->path();
      const std::string rel = std::filesystem::relative(full, base).generic_string();
      const uint64_t sz = static_cast<uint64_t>(std::filesystem::file_size(full));
      entries.emplace_back(rel, sz);
    }
    const std::string outpath = std::string(archive_path) + ".idx.json";
    std::ofstream ofs(outpath, std::ios::binary);
    if (!ofs.is_open()) throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, outpath);
    ofs << "{\n";
    for (size_t i = 0; i < entries.size(); ++i)
    {
      ofs << "  \"" << entries[i].first << "\": " << entries[i].second;
      if (i + 1 < entries.size()) ofs << ",\n";
      else ofs << "\n";
    }
    ofs << "}\n";
    ofs.close();
    return;
  }

  // Otherwise, open zip and stat entries
  int err = 0;
  zip_t* za = zip_open(archive_path.c_str(), ZIP_RDONLY, &err);
  if (!za)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to open zip archive", archive_path);
  }

  zip_int64_t num = zip_get_num_entries(za, 0);

  // Build JSON sidecar in memory
  std::string json = "{\n";
  bool first = true;
  for (zip_uint64_t i = 0; i < static_cast<zip_uint64_t>(num); ++i)
  {
    zip_stat_t st;
    if (zip_stat_index(za, i, 0, &st) == 0)
    {
      const char* name = st.name;
      if (!name) continue;
      if (!first) json += ",\n";
      first = false;
      json += "  \"";
      json += name;
      json += "\": ";
      json += std::to_string(static_cast<uint64_t>(st.size));
    }
  }
  json += "\n}\n";

  // Close the read-only handle and reopen archive for modification to add the sidecar
  zip_close(za);

  int err2 = 0;
  zip_t* za2 = zip_open(archive_path.c_str(), ZIP_CREATE, &err2);
  if (!za2)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to open zip archive for sidecar write", archive_path);
  }

  const std::string sidecar_name_str = std::filesystem::path(std::string(archive_path)).filename().string() + ".idx.json";
  const char* sidecar_name = sidecar_name_str.c_str();
  zip_source_t* src = zip_source_buffer(za2, json.data(), json.size(), 0);
  if (!src)
  {
    zip_close(za2);
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to create zip source for sidecar", archive_path);
  }

  zip_int64_t side_idx = zip_name_locate(za2, sidecar_name, 0);
  if (side_idx >= 0)
  {
    if (zip_replace(za2, side_idx, src) < 0)
    {
      zip_source_free(src);
      zip_close(za2);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "zip_replace failed for sidecar: " + String(zip_strerror(za2)), "");
    }
  }
  else
  {
    if (zip_file_add(za2, sidecar_name, src, ZIP_FL_OVERWRITE | ZIP_FL_ENC_UTF_8) < 0)
    {
      zip_source_free(src);
      zip_close(za2);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "zip_file_add failed for sidecar: " + String(zip_strerror(za2)), "");
    }
  }

  // Store sidecar uncompressed
  zip_int64_t new_idx = zip_name_locate(za2, sidecar_name, 0);
  if (new_idx >= 0)
  {
    zip_set_file_compression(za2, static_cast<zip_uint64_t>(new_idx), ZIP_CM_STORE, 0);
  }

  if (zip_close(za2) < 0)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to finalize zip archive while writing sidecar", archive_path);
  }
#else
  (void)archive_path;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}

String ZipArchiveFile::extractEntryToTempFile(const String& archive_path, const String& entry_name, std::unique_ptr<File::TempDir>& temp_dir)
{
#if defined(OPENMS_HAVE_LIBZIP)
  // If archive_path is actually a directory (tests create a .oswpq directory),
  // treat it as an unpacked layout and return the matching file path directly.
  if (File::isDirectory(archive_path))
  {
    const std::filesystem::path base = std::filesystem::u8path(std::string(archive_path));
    const std::filesystem::path entry_path = (base / std::filesystem::u8path(std::string(entry_name))).lexically_normal();
    if (!std::filesystem::exists(entry_path))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, entry_path.string());
    }
    return String(entry_path.string());
  }

  if (!File::readable(archive_path))
  {
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive_path);
  }

  int err = 0;
  zip_t* za = zip_open(archive_path.c_str(), ZIP_RDONLY, &err);
  if (!za)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to open zip archive", archive_path);
  }

  zip_int64_t idx = zip_name_locate(za, entry_name.c_str(), 0);
  if (idx < 0)
  {
    zip_close(za);
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, entry_name);
  }

  zip_file_t* zf = zip_fopen_index(za, idx, 0);
  if (!zf)
  {
    zip_close(za);
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to open zip entry", entry_name);
  }

  if (!temp_dir) temp_dir = std::make_unique<File::TempDir>();
  const String base = temp_dir->getPath();

  // construct output path and create parent dirs
  const std::filesystem::path entry_path = std::filesystem::u8path(std::string(entry_name));
  if (entry_path.is_absolute())
  {
    zip_fclose(zf);
    zip_close(za);
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Zip entry has absolute path", entry_name);
  }

  const std::filesystem::path outpath = (std::filesystem::u8path(std::string(base)) / entry_path).lexically_normal();
  const std::string base_str = std::filesystem::u8path(std::string(base)).string();
  const std::string out_str = outpath.string();
  if (out_str.size() < base_str.size() || out_str.compare(0, base_str.size(), base_str) != 0)
  {
    zip_fclose(zf);
    zip_close(za);
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Zip entry would extract outside target directory", entry_name);
  }

  if (!entry_name.empty() && entry_name.back() == '/')
  {
    std::filesystem::create_directories(outpath);
    zip_fclose(zf);
    zip_close(za);
    return String(outpath.string());
  }

  if (outpath.has_parent_path()) std::filesystem::create_directories(outpath.parent_path());

  std::ofstream ofs(outpath.string(), std::ios::binary);
  if (!ofs.is_open())
  {
    zip_fclose(zf);
    zip_close(za);
    throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, outpath.string());
  }

  constexpr size_t BUF_SIZE = 1 << 16;
  std::vector<char> buffer(BUF_SIZE);
  zip_int64_t nread = 0;
  while ((nread = zip_fread(zf, buffer.data(), buffer.size())) > 0)
  {
    ofs.write(buffer.data(), static_cast<std::streamsize>(nread));
    if (ofs.fail())
    {
      zip_fclose(zf);
      zip_close(za);
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, outpath.string());
    }
  }
  if (nread < 0)
  {
    zip_fclose(zf);
    zip_close(za);
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to read zip entry", entry_name);
  }

  ofs.close();
  zip_fclose(zf);
  zip_close(za);

  return String(outpath.string());
#else
  (void)archive_path;
  (void)entry_name;
  (void)temp_dir;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}

} // namespace OpenMS
