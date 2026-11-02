/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <zip.h>

#include "mzpeak/exception.h"
#include "mzpeak/io/zip.h"
#include "mzpeak/util/compat.h"

namespace MzPeak::IO {

/******************************************************************************/
/**
 * A streaming source for files in a zip archive.
 */
class ZipFile_ final : public MzPeak::IO::File {
public:
  // Each member owns its OWN archive handle.  libzip's zip_fseek on stored
  // members shares the underlying archive read position, so seeking one
  // member while another member of the same handle is open fails ("unable
  // to seek").  The reader keeps spectra_data and spectra_peaks open at the
  // same time (RDR-3), so each gets an independent handle.
  ZipFile_(zip_t* archive, zip_file_t* file, std::size_t size, fs::path path)
      : impl_(std::make_shared<Impl>(archive, file, size, path))
  {
  }

  std::string name() const override { return Util::narrow(impl_->path_); }

  std::size_t size() const override { return impl_->size_; }

  std::optional<std::size_t> read(uint8_t* buf, std::size_t size) override
  {
    if (buf == nullptr || size == 0 || !is_open()) return {};

    zip_int64_t n = zip_fread(impl_->file_, buf, size);

    if (n <= 0) {
      return {};
    } else {
      return n;
    }
  }

  std::optional<std::size_t> tell() const override
  {
    zip_int64_t n = zip_ftell(impl_->file_);

    if (n < 0) {
      return {};
    } else {
      return n;
    }
  }

  bool seek(std::size_t pos) override
  {
    zip_int8_t errnum = zip_fseek(impl_->file_, pos, SEEK_SET);
    return errnum == 0;
  }

  void close() override { impl_->close(); }
  bool is_open() const override { return impl_->file_ != nullptr; }

private:
  class Impl {
  public:
    Impl(zip_t* archive, zip_file_t* file, std::size_t size, fs::path path)
        : size_(size)
        , archive_(archive)
        , file_(file)
        , path_(path)
    {
    }
    ~Impl() { close(); }

    void close()
    {
      if (file_ != nullptr) {
        zip_fclose(file_);
        file_ = nullptr;
      }
      if (archive_ != nullptr) {
        zip_close(archive_);
        archive_ = nullptr;
      }
    }

    std::size_t size_;
    zip_t* archive_;
    zip_file_t* file_;
    fs::path path_;

  private:
    Impl& operator=(const Impl&) = delete;
    Impl(const Impl&) = delete;
  };

  std::shared_ptr<Impl> impl_;
};

/******************************************************************************/
struct Zip::Impl {

  /**************************************************************************/
  Impl(const fs::path& path)
      : archive(nullptr)
      , path_(path)
  {
    int errnum{};
    archive = zip_open(Util::narrow(path).c_str(), ZIP_RDONLY, &errnum);

    if (archive == nullptr) {
      error("failed to open zip archive ", errnum);
    }
  }

  /**************************************************************************/
  ~Impl()
  {
    if (archive != nullptr) {
      zip_close(archive);
      archive = nullptr;
    }
  }

  /**************************************************************************/
  void error(const std::string& msg, const std::optional<int>& errnum)
  {
    std::string m(msg);
    zip_error_t error;
    zip_error_t* error_ptr;

    if (errnum) {
      zip_error_init_with_code(&error, *errnum);
      error_ptr = &error;
    } else {
      error_ptr = zip_get_error(archive);
    }

    m += zip_error_strerror(error_ptr);
    throw(InvalidFormatError(m));
  }

  /**************************************************************************/
  void error_open(const fs::path& member, const std::optional<int>& errnum)
  {
    std::string msg("failed to open member '" + member.string() +
                    "' in zip archive '" + path_.string() + "': ");
    error(msg, errnum);
  }

  /**************************************************************************/
  zip_t* archive;
  fs::path path_;

private:
  Impl(const Impl&) = delete;
  Impl& operator=(const Impl&) = delete;
};

/******************************************************************************/
Zip::Zip(const fs::path& path)
    : impl_(std::make_unique<Impl>(path))
{
}

/******************************************************************************/
// NOTE: This is needed due to the pimpl pattern and the `Impl` type
// not being complete in the header file :(
Zip::~Zip() = default;

/******************************************************************************/
std::vector<fs::path> Zip::list()
{
  zip_int64_t num = zip_get_num_entries(impl_->archive, ZIP_FL_UNCHANGED);
  zip_stat_t stat;
  int errnum;

  std::vector<fs::path> files;
  files.reserve(num);

  for (zip_int64_t index = 0; index < num; ++index) {
    errnum = zip_stat_index(impl_->archive, index, ZIP_FL_UNCHANGED, &stat);

    if (errnum == 0 && stat.valid & ZIP_STAT_NAME) {
      files.push_back(stat.name);
    }
  }

  return files;
}

/******************************************************************************/
std::unique_ptr<MzPeak::IO::File> Zip::read_file(const fs::path& path)
{
  // Open a dedicated archive handle for this member so its seeks are
  // independent of any other member the reader keeps open (RDR-26).
  int errnum{};
  zip_t* archive = zip_open(Util::narrow(impl_->path_).c_str(), ZIP_RDONLY, &errnum);
  if (archive == nullptr) {
    impl_->error_open(path, errnum);
  }

  zip_stat_t stat;
  const std::string member = Util::narrow(path);
  if (zip_stat(archive, member.c_str(), ZIP_FL_UNCHANGED, &stat) != 0 ||
      !(stat.valid & ZIP_STAT_SIZE)) {
    zip_close(archive);
    impl_->error_open(path, {});
  }

  zip_file_t* file = zip_fopen(archive, member.c_str(), 0);
  if (file == nullptr) {
    zip_close(archive);
    impl_->error_open(path, {});
  }

  return std::make_unique<ZipFile_>(archive, file, stat.size, path);
}

} // namespace MzPeak::IO
