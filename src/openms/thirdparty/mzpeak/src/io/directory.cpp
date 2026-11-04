/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <fstream>
#include <memory>

#include "mzpeak/io/directory.h"

namespace MzPeak::IO {

/******************************************************************************/
// Thin wrapper around `std::fstream`.
class DirFile_ final : public MzPeak::IO::File {
public:
  DirFile_(const fs::path& path)
      : path_(path)
      // Read-only and binary: the default in|out text mode turns \r\n into
      // \n on Windows and stops at 0x1A, which corrupts every Parquet file.
      , stream_(path, std::ios::in | std::ios::binary)
  {
  }

  ~DirFile_() = default;

  std::string name() const override { return path_.string(); }

  std::size_t size() const override { return fs::file_size(path_); }

  std::optional<std::size_t> read(uint8_t* buffer, std::size_t size) override
  {
    if (stream_.good() && !stream_.eof()) {
      stream_.read(reinterpret_cast<char*>(buffer), size);
      return stream_.gcount();
    } else {
      return {};
    }
  }

  std::optional<std::size_t> tell() const override
  {
    // `tellg` should be const, but it's not marked that way.
    std::fstream& s(const_cast<std::fstream&>(stream_));
    return s.tellg();
  }

  bool seek(std::size_t pos) override
  {
    stream_.seekg(pos);
    return !stream_.fail();
  }

  void close() override { stream_.close(); }

  bool is_open() const override { return stream_.is_open(); }

private:
  fs::path path_;
  std::fstream stream_;
};

/******************************************************************************/
std::vector<fs::path> Directory::list()
{
  std::vector<fs::path> res;

  for (auto const& entry : fs::directory_iterator{path_}) {
    if (!entry.is_directory()) {
      res.push_back(entry.path().lexically_relative(path_));
    }
  }

  return res;
}

/******************************************************************************/
std::unique_ptr<MzPeak::IO::File> Directory::read_file(const fs::path& name)
{
  fs::path path(path_ / name.lexically_normal());
  return std::make_unique<DirFile_>(path.string());
}

} // namespace MzPeak::IO
