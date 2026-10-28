/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <memory>
#include <optional>

namespace MzPeak::IO {

/******************************************************************************/
/**
 * An interface for files that can be read from.
 */
class File {
public:
  /// Destructor.
  virtual ~File() {}

  /**
   * The name of this file.
   */
  virtual std::string name() const = 0;

  /**
   * Return the size of the file, in bytes.
   */
  virtual std::size_t size() const = 0;

  /**
   * Return the position of the read pointer.
   */
  virtual std::optional<std::size_t> tell() const = 0;

  /**
   * Read from the file.
   *
   * If successful, return the number of bytes actually read.
   */
  virtual std::optional<std::size_t> read(uint8_t* buffer, std::size_t size) = 0;

  /**
   * Move the read pointer to the given file position.
   *
   * Returns `true` if the seek was successful.
   */
  virtual bool seek(std::size_t pos) = 0;

  /**
   * Close the file.
   */
  virtual void close() = 0;

  /**
   * Returns `true` if the file is open.
   */
  virtual bool is_open() const = 0;

protected:
  // Prevent construction.
  File() = default;

private:
  // Prevent copying.
  File(const File&) = default;
};

/******************************************************************************/
/**
 * Turn a `File` object into an `std::istream`.
 */
std::unique_ptr<std::istream> to_istream(std::unique_ptr<File>);

} // namespace MzPeak::IO
