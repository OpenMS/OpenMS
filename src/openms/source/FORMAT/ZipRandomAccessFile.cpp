// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZipRandomAccessFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <filesystem>
#include <vector>
#include <arrow/buffer.h>

#if __has_include(<zip.h>)
#include <zip.h>
#define OPENMS_HAVE_LIBZIP 1
#endif

namespace OpenMS
{

// Implements an Arrow RandomAccessFile backed by a libzip entry. The
// implementation mirrors the working test helper in Libzip_test.cpp and
// relies on zip_fseek/zip_ftell/zip_fread for seekable access on stored
// (uncompressed) entries.
class ZipEntryRandomAccessFile : public arrow::io::RandomAccessFile
{
public:
	ZipEntryRandomAccessFile(zip_t* archive, zip_file_t* entry, int64_t entry_size)
		: archive_(archive), entry_(entry), entry_size_(entry_size), closed_(false) {}

	~ZipEntryRandomAccessFile() override
	{
		if (!closed_)
		{
			if (entry_) zip_fclose(entry_);
			if (archive_) zip_close(archive_);
		}
	}

	arrow::Status Close() override
	{
		if (!closed_)
		{
			if (entry_) { zip_fclose(entry_); entry_ = nullptr; }
			if (archive_) { zip_close(archive_); archive_ = nullptr; }
			closed_ = true;
		}
		return arrow::Status::OK();
	}

	bool closed() const override { return closed_; }

	arrow::Result<int64_t> Tell() const override
	{
		if (closed_) return arrow::Status::Invalid("File is closed");
		zip_int64_t pos = zip_ftell(entry_);
		if (pos < 0) return arrow::Status::IOError("zip_ftell failed");
		return static_cast<int64_t>(pos);
	}

	arrow::Result<int64_t> GetSize() override { return entry_size_; }

	arrow::Status Seek(int64_t position) override
	{
		if (closed_) return arrow::Status::Invalid("File is closed");
		if (position < 0 || position > entry_size_) return arrow::Status::IOError("Seek out of bounds");
		if (zip_fseek(entry_, position, SEEK_SET) < 0) return arrow::Status::IOError("zip_fseek failed");
		return arrow::Status::OK();
	}

	arrow::Result<int64_t> Read(int64_t nbytes, void* out) override
	{
		if (closed_) return arrow::Status::Invalid("File is closed");
		if (nbytes < 0) return arrow::Status::Invalid("nbytes must be non-negative");
		if (nbytes == 0) return int64_t{0};
		zip_int64_t bytes_read = zip_fread(entry_, out, static_cast<zip_uint64_t>(nbytes));
		if (bytes_read < 0) return arrow::Status::IOError("zip_fread failed");
		return static_cast<int64_t>(bytes_read);
	}

	arrow::Result<std::shared_ptr<arrow::Buffer>> Read(int64_t nbytes) override
	{
		if (closed_) return arrow::Status::Invalid("File is closed");
		if (nbytes < 0) return arrow::Status::Invalid("nbytes must be non-negative");
		if (nbytes == 0) return arrow::Buffer::Copy(nullptr, 0);

		// Read into a temporary vector then copy into an Arrow buffer. This avoids
		// relying on AllocateResizableBuffer which may not be available in all
		// Arrow versions present in build environments.
		std::vector<uint8_t> tmp(static_cast<size_t>(nbytes));
		ARROW_ASSIGN_OR_RAISE(int64_t bytes_read, Read(nbytes, tmp.data()));
		if (bytes_read < 0) return arrow::Status::IOError("Read failed");
		// Resize to actual bytes read to avoid returning uninitialized data when
		// a short read (e.g., EOF) occurred.
		tmp.resize(static_cast<size_t>(bytes_read));
		auto buf = arrow::Buffer::FromVector(std::move(tmp));
		return buf;
	}

	arrow::Result<int64_t> ReadAt(int64_t position, int64_t nbytes, void* out) override
	{
		ARROW_RETURN_NOT_OK(Seek(position));
		return Read(nbytes, out);
	}

	arrow::Result<std::shared_ptr<arrow::Buffer>> ReadAt(int64_t position, int64_t nbytes) override
	{
		ARROW_RETURN_NOT_OK(Seek(position));
		return Read(nbytes);
	}

private:
	zip_t* archive_;
	zip_file_t* entry_;
	int64_t entry_size_;
	bool closed_;
};


arrow::Result<std::shared_ptr<arrow::io::RandomAccessFile>> ZipRandomAccessFile::Open(const String& archive_path,
																																											const String& entry_name,
											std::unique_ptr<File::TempDir>& temp_dir)
{
    // temp_dir is only used on the extraction fallback path; silence unused-parameter
    // warnings when building with libzip support enabled.
    (void)temp_dir;
	// If the archive_path is a directory, return a plain ReadableFile for the path on disk
	if (File::isDirectory(archive_path))
	{
		// Build the entry path relative to the provided directory
		const std::filesystem::path base{std::string(archive_path)};
		const std::filesystem::path entry_path = (base / std::filesystem::path(std::string(entry_name))).lexically_normal();
		if (!std::filesystem::exists(entry_path))
		{
			return arrow::Status::IOError("Entry not found in directory: " + entry_path.string());
		}
		return arrow::io::ReadableFile::Open(entry_path.string());
	}

#if !defined(OPENMS_HAVE_LIBZIP)
	return arrow::Status::NotImplemented("libzip not available, cannot open ZIP entry directly");
#else
	int err = 0;
	zip_t* za = zip_open(archive_path.c_str(), ZIP_RDONLY, &err);
	if (!za)
	{
		return arrow::Status::IOError("Failed to open zip archive: " + std::string(archive_path));
	}

	zip_stat_t st;
	zip_stat_init(&st);
	if (zip_stat(za, entry_name.c_str(), 0, &st) != 0)
	{
		zip_close(za);
		return arrow::Status::IOError("Entry not found in ZIP: " + std::string(entry_name));
	}

	zip_file_t* entry = zip_fopen(za, entry_name.c_str(), 0);
	if (!entry)
	{
		zip_close(za);
		return arrow::Status::IOError("Failed to open ZIP entry: " + std::string(entry_name));
	}

	auto wrapper = std::make_shared<ZipEntryRandomAccessFile>(za, entry, static_cast<int64_t>(st.size));
	return std::static_pointer_cast<arrow::io::RandomAccessFile>(wrapper);
#endif
}

} // namespace OpenMS
