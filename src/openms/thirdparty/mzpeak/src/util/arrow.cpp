/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <arrow/buffer.h>
#include <arrow/io/api.h>
#include <arrow/io/buffered.h>
#include <memory>
#include <parquet/properties.h>

#include "mzpeak/exception.h"
#include "mzpeak/util/arrow.h"
#include "mzpeak/util/decode_pool.h"

namespace MzPeak::Util {

/******************************************************************************/
/**
 * Arrow file access.
 */
class ArrowFile_ final : public arrow::io::RandomAccessFile {
public:
  /// Constructor.
  ArrowFile_(std::shared_ptr<IO::File> file)
      : file_(std::move(file))
  {
  }

  /// Destructor.
  ~ArrowFile_() = default;

  /// Return the total file size in bytes.
  arrow::Result<int64_t> GetSize() override { return file_->size(); }

  /// Seek in file/stream.
  arrow::Status Seek(int64_t position) override
  {
    if (!file_->seek(position)) {
      std::string msg("unable to seek");
      return arrow::Status(arrow::StatusCode::IOError, msg);
    }

    return arrow::Status::OK();
  }

  /// Report the current position.
  arrow::Result<int64_t> Tell() const override
  {
    std::optional<std::size_t> n = file_->tell();

    if (n.has_value()) {
      return arrow::Result<int64_t>(n.value());
    } else {
      // Arrow error result:
      return arrow::Result<int64_t>();
    }
  }

  /// Read data from current file position.
  arrow::Result<int64_t> Read(int64_t nbytes, void* out) override
  {
    std::optional<std::size_t> n = file_->read(static_cast<uint8_t*>(out), nbytes);

    if (n.has_value()) {
      return arrow::Result<int64_t>(n.value());
    } else {
      return arrow::Result<int64_t>();
    }
  }

  /// Read into a buffer.
  arrow::Result<std::shared_ptr<arrow::Buffer>> Read(int64_t nbytes) override
  {
    using arrow_buffer_t = std::shared_ptr<arrow::Buffer>;

    if (nbytes < 0) {
      return arrow::Status(arrow::StatusCode::Invalid, "negative read size");
    }

    // Each buffer-returning read MUST own its storage.  Returning a shared
    // scratch buffer that a later read overwrites would corrupt data Arrow
    // still holds (e.g. a footer buffer kept while metadata is read).
    //
    // From the DecodePool: these are the column chunks Parquet pre-buffers,
    // read on Arrow's I/O threads, freed on a reader's, and kept by each
    // reader until its next decode (~4 MB each).  In Arrow's default pool
    // (jemalloc) they cost ~1 GB of RSS beyond the live ones at 128 threads.
    arrow::Result<std::unique_ptr<arrow::ResizableBuffer>> alloc(
        arrow::AllocateResizableBuffer(nbytes, &decode_pool()));
    if (!alloc.ok()) throw ParquetError(alloc.status().ToString());
    std::shared_ptr<arrow::ResizableBuffer> buffer(std::move(alloc.ValueOrDie()));

    std::optional<std::size_t> n = file_->read(buffer->mutable_data(), nbytes);

    if (!n.has_value()) {
      return arrow::Result<arrow_buffer_t>();
    }

    // Report the number of bytes actually read, not the requested size.
    // Otherwise a short read (any file smaller than the footer read size)
    // hands Arrow an over-long buffer and it looks for the Parquet footer
    // magic at the wrong offset ("magic bytes not found").
    arrow::Status status = buffer->Resize(static_cast<int64_t>(*n), false);
    if (!status.ok()) throw ParquetError(status.ToString());

    return arrow::Result<arrow_buffer_t>(buffer);
  }

  /// Close the file/stream.
  arrow::Status Close() override
  {
    file_->close();
    return arrow::Status::OK();
  }

  /// Return `true` if the file/stream is closed.
  bool closed() const override { return !file_->is_open(); }

private:
  std::shared_ptr<IO::File> file_;
};

/******************************************************************************/
struct Arrow::Impl {
  Impl(std::unique_ptr<IO::File> file)
      : file_(std::move(file))
      , reader_(std::make_shared<ArrowFile_>(file_))
  {
  }

  std::shared_ptr<IO::File> file_;
  std::shared_ptr<ArrowFile_> reader_;
};

/******************************************************************************/
Arrow::Arrow(std::unique_ptr<IO::File> file)
    : impl_(std::make_unique<Impl>(std::move(file)))
{
}

/******************************************************************************/
Arrow::~Arrow() = default;

/******************************************************************************/
std::shared_ptr<Arrow::random_access_t> Arrow::reader() const
{
  return impl_->reader_;
}

} // namespace MzPeak::Util
