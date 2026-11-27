/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <boost/iostreams/concepts.hpp>
#include <boost/iostreams/stream.hpp>
#include <iostream>
#include <parquet/properties.h>

#include "mzpeak/io/file.h"

namespace MzPeak::IO {

/******************************************************************************/
class Source final : public boost::iostreams::source {
public:
  Source(std::shared_ptr<File> file)
      : file_(file)
  {
  }

  std::streamsize read(char* buf, std::streamsize size)
  {
    if (size == 0 || buf == nullptr || !file_->is_open()) return -1;
    std::optional<std::size_t> n =
        file_->read(reinterpret_cast<uint8_t*>(buf), size);

    if (n.has_value()) {
      return n.value();
    } else {
      return -1;
    }
  }

  void close() { file_->close(); }
  void close(std::ios_base::openmode&) { close(); }

private:
  std::shared_ptr<File> file_;
};

/******************************************************************************/
class Istream : public std::istream {
public:
  Istream(const Source& s)
      : std::istream(&buf_)
      , buf_(s)
  {
  }

private:
  boost::iostreams::stream_buffer<Source> buf_;
};

/******************************************************************************/
std::unique_ptr<std::istream> to_istream(std::unique_ptr<File> r)
{
  Source source(std::move(r));
  return std::make_unique<Istream>(source);
}

} // namespace MzPeak::IO
