/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <zip.h>

#include "mzpeak/zip_buffer.h"

namespace MzPeak {

/******************************************************************************/
/**
 * RDR-22 — a streaming source for a file in an in-memory zip archive.
 *
 * Each member owns its OWN archive handle (a `zip_t*` opened from a
 * `zip_source_buffer` over the owned bytes) and its `zip_file_t*`, closing
 * both on destruction.  The shared independent-handle rationale from
 * `Zip`/`ZipFile_` (RDR-26) applies: seeking one member must not disturb
 * another member the reader keeps open at the same time, so each member
 * gets its own handle.  Closing the handle frees the libzip source but NOT
 * the underlying buffer (created with `freep=0`); the buffer is owned by the
 * `ZipBuffer` archive, which outlives every `File`.
 */
class ZipBufferFile_ final : public MzPeak::File {
public:
  ZipBufferFile_(zip_t* archive, zip_file_t* file, std::size_t size, fs::path path)
      : impl_(std::make_shared<Impl>(archive, file, size, path))
  {
  }

  std::string name() const { return impl_->path_; }

  std::size_t size() const { return impl_->size_; }

  std::optional<std::size_t> read(uint8_t* buf, std::size_t size)
  {
    if (buf == nullptr || size == 0 || !is_open()) return {};

    zip_int64_t n = zip_fread(impl_->file_, buf, size);

    if (n <= 0) {
      return {};
    } else {
      return n;
    }
  }

  std::optional<std::size_t> tell() const
  {
    zip_int64_t n = zip_ftell(impl_->file_);

    if (n < 0) {
      return {};
    } else {
      return n;
    }
  }

  bool seek(std::size_t pos)
  {
    zip_int8_t errnum = zip_fseek(impl_->file_, pos, SEEK_SET);
    return errnum == 0;
  }

  void close() { impl_->close(); }
  bool is_open() const { return impl_->file_ != nullptr; }

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
    Impl(const Impl&) = default;
  };

  std::shared_ptr<Impl> impl_;
};

/******************************************************************************/
struct ZipBuffer::Impl {

  /**************************************************************************/
  // RDR-22 — the archive owns the bytes.  This vector must outlive every
  // `zip_t*` opened over it (each `zip_source_buffer` is created with
  // `freep=0`, so libzip never frees these bytes and reads them in place).
  // Since the archive owns the buffer and outlives the Index/Files, this
  // lifetime contract holds.
  Impl(std::vector<std::byte> bytes)
      : data_(std::move(bytes))
  {
  }

  /**************************************************************************/
  // RDR-22 — open a fresh archive handle over a libzip source backed by the
  // owned buffer.  On `zip_open_from_source` failure libzip does NOT consume
  // the source, so we must free it ourselves; on success the zip takes
  // ownership of the source.
  zip_t* open()
  {
    zip_error_t error;
    zip_error_init(&error);

    zip_source_t* src =
        zip_source_buffer_create(data_.data(), data_.size(), /*freep=*/0, &error);
    if (src == nullptr) {
      std::string msg("failed to create zip source from buffer: ");
      msg += zip_error_strerror(&error);
      zip_error_fini(&error);
      throw(std::invalid_argument(msg));
    }

    zip_t* archive = zip_open_from_source(src, ZIP_RDONLY, &error);
    if (archive == nullptr) {
      zip_source_free(src);
      std::string msg("failed to open zip archive from buffer: ");
      msg += zip_error_strerror(&error);
      zip_error_fini(&error);
      throw(std::invalid_argument(msg));
    }

    zip_error_fini(&error);
    return archive;
  }

  /**************************************************************************/
  // Closes @p archive (taking its error string first) before throwing, so a
  // failed open does not leak the handle.
  [[noreturn]] void error_open(zip_t* archive, const fs::path& path)
  {
    std::string msg("failed to open file in zip archive ");
    msg += path.string() + ": ";
    msg += zip_error_strerror(zip_get_error(archive));
    zip_close(archive);
    throw(std::invalid_argument(msg));
  }

  /**************************************************************************/
  std::vector<std::byte> data_;
};

/******************************************************************************/
ZipBuffer::ZipBuffer(std::vector<std::byte> bytes)
    : impl_(std::make_unique<Impl>(std::move(bytes)))
{
}

/******************************************************************************/
// NOTE: This is needed due to the pimpl pattern and the `Impl` type
// not being complete in the header file :(
ZipBuffer::~ZipBuffer() = default;

/******************************************************************************/
std::vector<fs::path> ZipBuffer::list()
{
  zip_t* archive = impl_->open();

  zip_int64_t num = zip_get_num_entries(archive, ZIP_FL_UNCHANGED);
  zip_stat_t stat;
  int errnum;

  std::vector<fs::path> files;
  files.reserve(num);

  for (zip_int64_t index = 0; index < num; ++index) {
    errnum = zip_stat_index(archive, index, ZIP_FL_UNCHANGED, &stat);

    if (errnum == 0 && stat.valid & ZIP_STAT_NAME) {
      files.push_back(stat.name);
    }
  }

  zip_close(archive);
  return files;
}

/******************************************************************************/
std::unique_ptr<MzPeak::File> ZipBuffer::read_file(const fs::path& path)
{
  // RDR-22 — open a dedicated archive handle over the owned buffer for this
  // member so its seeks are independent of any other member the reader keeps
  // open (mirrors RDR-26).
  zip_t* archive = impl_->open();

  zip_stat_t stat;
  if (zip_stat(archive, path.c_str(), ZIP_FL_UNCHANGED, &stat) != 0 ||
      !(stat.valid & ZIP_STAT_SIZE)) {
    impl_->error_open(archive, path); // throws
  }

  zip_file_t* file = zip_fopen(archive, path.c_str(), 0);
  if (file == nullptr) {
    impl_->error_open(archive, path); // throws
  }

  return std::make_unique<ZipBufferFile_>(archive, file, stat.size, path);
}

} // namespace MzPeak
