// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DecompressingInputStream.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <bzlib.h>
#include <zlib.h>

#include <algorithm>
#include <climits>
#include <cstring>
#include <vector>

namespace OpenMS
{
  struct DecompressingInputStream::Impl
  {
    std::unique_ptr<InputStream> source;
    Format format;
    std::string name;
    z_stream z{};
    bz_stream bz{};
    bool active = false;        ///< a decompressor is initialized
    std::vector<unsigned char> buffer = std::vector<unsigned char>(64 * 1024); ///< compressed bytes from the source
    Size next = 0;              ///< the first unconsumed byte in buffer
    Size available = 0;         ///< the number of unconsumed bytes from next on
    bool source_end = false;    ///< the source has no more bytes
    bool finished = false;      ///< no more decompressed bytes follow
    UInt64 position = 0;

    const char* formatName() const { return format == Format::GZIP ? "gzip" : "bzip2"; }

    [[noreturn]] void fail(const std::string& message) const
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name, message);
    }

    void begin()
    {
      // windowBits 15 + 16: gzip format
      const bool ok = format == Format::GZIP ? inflateInit2(&z, 15 + 16) == Z_OK : BZ2_bzDecompressInit(&bz, 0, 0) == BZ_OK;
      if (!ok) { fail(std::string("Cannot initialize the ") + formatName() + " decompression"); }
      active = true;
    }

    void end()
    {
      if (!active) { return; }
      if (format == Format::GZIP) { inflateEnd(&z); }
      else { BZ2_bzDecompressEnd(&bz); }
      active = false;
    }

    // Makes at least @p count unconsumed bytes available, unless the source ends first; whether it did
    bool ensure(Size count)
    {
      while (available < count && !source_end)
      {
        if (available == 0)
        {
          next = 0; // all consumed: fill the whole buffer
        }
        else if (next + available == buffer.size())
        {
          std::memmove(buffer.data(), buffer.data() + next, available);
          next = 0;
        }
        const Size read = source->readBytes(buffer.data() + next + available, buffer.size() - next - available);
        if (read == 0) { source_end = true; }
        available += read;
      }
      return available >= count;
    }

    // Decompresses into [to, to + space); the number of bytes written. Sets finished at the end of the data.
    Size decompress(unsigned char* to, Size space)
    {
      if (available == 0 && !ensure(1))
      {
        fail(std::string("Truncated ") + formatName() + " data");
      }
      const Size in = std::min<Size>(available, UINT_MAX);
      const Size out = std::min<Size>(space, UINT_MAX);
      Size consumed = 0;
      Size produced = 0;
      bool member_end = false;
      if (format == Format::GZIP)
      {
        z.next_in = buffer.data() + next;
        z.avail_in = static_cast<uInt>(in);
        z.next_out = to;
        z.avail_out = static_cast<uInt>(out);
        const int ret = inflate(&z, Z_NO_FLUSH);
        consumed = in - z.avail_in;
        produced = out - z.avail_out;
        member_end = ret == Z_STREAM_END;
        if (ret != Z_OK && !member_end) { fail("Corrupt gzip data" + std::string(z.msg ? std::string(": ") + z.msg : "")); }
      }
      else
      {
        bz.next_in = reinterpret_cast<char*>(buffer.data() + next);
        bz.avail_in = static_cast<unsigned int>(in);
        bz.next_out = reinterpret_cast<char*>(to);
        bz.avail_out = static_cast<unsigned int>(out);
        const int ret = BZ2_bzDecompress(&bz);
        consumed = in - bz.avail_in;
        produced = out - bz.avail_out;
        member_end = ret == BZ_STREAM_END;
        if (ret != BZ_OK && !member_end) { fail("Corrupt bzip2 data"); }
      }
      next += consumed;
      available -= consumed;
      if (member_end)
      {
        // Concatenated members are read in sequence; other trailing bytes are ignored (as gzread() does)
        if (ensure(3) || available > 0)
        {
          if (DecompressingInputStream::detect(buffer.data() + next, available) == format)
          {
            if (format == Format::GZIP) { inflateReset(&z); }
            else { end(); begin(); }
            return produced;
          }
        }
        end();
        finished = true;
      }
      else if (consumed == 0 && produced == 0)
      {
        fail(std::string("Corrupt ") + formatName() + " data: the decompression made no progress");
      }
      return produced;
    }
  };

  std::optional<DecompressingInputStream::Format> DecompressingInputStream::detect(const unsigned char* head, Size size)
  {
    if (size >= 2 && head[0] == 0x1F && head[1] == 0x8B) { return Format::GZIP; }
    if (size >= 3 && head[0] == 'B' && head[1] == 'Z' && head[2] == 'h') { return Format::BZIP2; }
    return std::nullopt;
  }

  DecompressingInputStream::DecompressingInputStream(std::unique_ptr<InputStream> source, Format format, const std::string& name) :
    impl_(std::make_unique<Impl>())
  {
    impl_->source = std::move(source);
    impl_->format = format;
    impl_->name = name;
    impl_->begin();
  }

  DecompressingInputStream::~DecompressingInputStream()
  {
    impl_->end();
  }

  UInt64 DecompressingInputStream::curPos() const
  {
    return impl_->position;
  }

  void DecompressingInputStream::cancel()
  {
    impl_->source->cancel();
  }

  Size DecompressingInputStream::readBytes(unsigned char* to_fill, Size max_to_read)
  {
    Size written = 0;
    while (written < max_to_read && !impl_->finished)
    {
      written += impl_->decompress(to_fill + written, max_to_read - written);
    }
    impl_->position += written;
    return written;
  }

} // namespace OpenMS
