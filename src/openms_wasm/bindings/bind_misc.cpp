// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Version, shared data, logging, threads, memory and the file system.

#include "binding_utils.h"

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/SYSTEM/File.h>

#include <emscripten/heap.h>
#include <emscripten/threading.h>
#include <emscripten/wasmfs.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <new>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace emscripten;

namespace OpenMS::WASM
{
  namespace
  {
    std::string version() { return VersionInfo::getVersion(); }
    std::string revision() { return VersionInfo::getRevision(); }
    std::string branch() { return VersionInfo::getBranch(); }
    std::string dataPath() { return File::getOpenMSDataPath(); }

    /// The global log stream of @p channel and the console stream it writes to by default.
    std::pair<Logger::LogStream*, std::ostream*> logChannel(const std::string& channel)
    {
      if (channel == "fatal") return {&getGlobalLogFatal(), &std::cerr};
      if (channel == "error") return {&getGlobalLogError(), &std::cerr};
      if (channel == "warn") return {&getGlobalLogWarn(), &std::cerr};
      if (channel == "info") return {&getGlobalLogInfo(), &std::cout};
      if (channel == "debug") return {&getGlobalLogDebug(), &std::cerr};
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Unknown log channel (fatal, error, warn, info, debug)", channel);
    }

    /// Enables or disables the console output (Module.print / Module.printErr) of a log channel.
    void setLogEnabled(const std::string& channel, bool enabled)
    {
      auto [log, console] = logChannel(channel);
      if (enabled && !log->hasStream(*console))
      {
        log->insert(*console);
      }
      else if (!enabled)
      {
        log->remove(*console);
      }
    }

    bool isLogEnabled(const std::string& channel)
    {
      auto [log, console] = logChannel(channel);
      return log->hasStream(*console);
    }

    /// Whether the module is built with threads, whether it can start them (browsers: on a
    /// cross-origin isolated page) and the number of threads OpenMS computes with.
    JS::ThreadingInfo threadingInfo()
    {
      val info = val::object();
#ifdef __EMSCRIPTEN_PTHREADS__
      info.set("threaded", true);
      info.set("available", emscripten_has_threading_support());
#else
      info.set("threaded", false);
      info.set("available", false);
#endif
#ifdef _OPENMP
      info.set("numThreads", omp_get_max_threads());
#else
      info.set("numThreads", 1);
#endif
      return info;
    }

    /// Whether the module includes ONNX Runtime (built with WITH_ONNX): machine-learning models
    /// such as PeptDeep's, e.g. for ProSE's 'peptdeep:enable'.
    bool hasONNXRuntime()
    {
#ifdef WITH_ONNX
      return true;
#else
      return false;
#endif
    }

    /// The size of the module's memory and the most it can grow to, in bytes. WebAssembly memory
    /// grows but never shrinks: after a computation, heapSize is its high-water mark.
    JS::MemoryInfo memoryInfo()
    {
      val info = val::object();
      info.set("heapSize", static_cast<double>(emscripten_get_heap_size()));
      info.set("heapMax", static_cast<double>(emscripten_get_heap_max()));
      return info;
    }

    //------------------------------------------------------------------------------
    // File system (WasmFS, src/openms_wasm/CMakeLists.txt)

    bool runsInNode()
    {
      const val process = val::global("process");
      return !process.isUndefined() && !process["versions"].isUndefined() && !process["versions"]["node"].isUndefined();
    }

    /// Mounts directory @p host_directory of the host at @p mount_point of the module's file system
    /// (created; its parent has to exist). Node.js only: OpenMS then reads and writes the files
    /// of the host directly, without a copy.
    void mountHostDirectory(const std::string& mount_point, const std::string& host_directory)
    {
      if (!runsInNode())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "mountHostDirectory() is available in Node.js only; in browsers, write the files with FS.writeFile()", mount_point);
      }
      backend_t backend = wasmfs_create_node_backend(host_directory.c_str());
      if (backend == nullptr)
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, host_directory);
      }
      const int error = wasmfs_create_directory(mount_point.c_str(), 0777, backend);
      if (error < 0)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mount_point, std::strerror(-error));
      }
    }

    //------------------------------------------------------------------------------
    // Files in parts: writeFileChunk(), readFileChunk() and allocateFile(), for files that are written or read in
    // chunks of a bounded size, e.g. a File of the browser staged with FileReaderSync. One copy each way (a typed
    // array set), unlike FS.write(), which copies byte by byte.
    //
    // WasmFS' memory backend keeps a file in one std::vector, and WasmFS is compiled without exceptions: an
    // allocation that fails inside it aborts (and would unwind through its locks). So the size of a file is checked
    // first (InvalidValue at 2 GiB), then the memory a write or a resize needs is allocated once here, and an error
    // (OutOfMemory) is thrown if the module's memory cannot hold it.

    /// A file descriptor, closed when it goes out of scope
    struct FileDescriptor
    {
      int fd = -1;
      ~FileDescriptor() { if (fd >= 0) ::close(fd); }
    };

    /// A byte offset or size from JavaScript: a non-negative integer (a safe one)
    off_t byteCount(double value, const char* what, const std::string& path)
    {
      if (!(value >= 0.0 && value <= 9007199254740991.0 && std::floor(value) == value))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::string(what) + " has to be a non-negative integer (" + path + ")",
                                      val(value).call<std::string>("toString")); // NaN, Infinity, 1e+21 as JavaScript writes them
      }
      return static_cast<off_t>(value);
    }

    /// The largest file of WasmFS' memory backend: its std::vector holds less than 2 GiB in wasm32 (max_size() is
    /// 2^31 - 1), and an allocation close to that aborts inside WasmFS (in a fresh module 1.99 GiB worked, 1.999 GiB
    /// aborted although the memory had room). Up to 2 GiB - 64 MiB, an allocation succeeds where the memory has room.
    constexpr off_t MAX_FILE_SIZE = (off_t(1) << 31) - (off_t(64) << 20);

    /// Throws InvalidValue if a file of @p size bytes does not fit into the memory file system (MAX_FILE_SIZE)
    void checkFileSize(off_t size, const std::string& path)
    {
      if (size <= MAX_FILE_SIZE) return;
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "a file of the module's file system has to be smaller than 2 GiB (at most " + std::to_string(MAX_FILE_SIZE)
        + " bytes): " + path, std::to_string(size));
    }

    /// The block that growing a file of @p size bytes to @p end bytes allocates: its std::vector grows to twice its
    /// capacity (taken as its size; allocateFile() makes them equal), or to @p end if that is larger, and at most to
    /// max_size(). Throws InvalidValue if that block would not fit (a file of more than 1 GiB that grows).
    off_t growthAllocation(off_t size, off_t end, const std::string& path)
    {
      const off_t block = std::min<off_t>(std::max<off_t>(end, 2 * size), (off_t(1) << 31) - 1);
      if (block > MAX_FILE_SIZE)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "growing the file " + path + " of " + std::to_string(size) + " bytes takes a block of " + std::to_string(block)
          + " bytes (it grows by doubling), and a block of the module's file system has to be smaller than 2 GiB: "
          "allocate the file at its final size first (allocateFile())", std::to_string(end));
      }
      return block;
    }

    /// Throws OutOfMemory unless the module's memory can hold @p bytes more (it grows if it can)
    void checkMemoryFor(size_t bytes)
    {
      if (bytes == 0) return;
      // volatile: the compiler may otherwise drop an allocation that is freed unused, and take it as successful
      void* volatile probe = std::malloc(bytes);
      if (probe == nullptr) throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, bytes);
      std::free(probe);
    }

    /// The bytes copied between JavaScript and a file at a time: a buffer of this size, not one of the whole chunk, so
    /// that a chunk takes no memory of the module (mimalloc does not always reuse the memory of large blocks)
    constexpr size_t COPY_BUFFER_SIZE = size_t(1) << 20;

    /// A buffer for the copies of a chunk of @p n bytes (at most COPY_BUFFER_SIZE)
    std::unique_ptr<uint8_t[]> copyBuffer(size_t n)
    {
      const size_t size = std::max<size_t>(1, std::min(n, COPY_BUFFER_SIZE));
      std::unique_ptr<uint8_t[]> buffer(new (std::nothrow) uint8_t[size]);
      if (!buffer) throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, size);
      return buffer;
    }

    /// Writes @p data (a Uint8Array) into the file @p path at byte @p offset; the file is created if missing. Writing
    /// past its end grows it: allocate the file first (allocateFile()) when its size is known, else it grows by
    /// doubling, which can take up to three times its size while it grows (and fails beyond 1 GiB).
    void writeFileChunk(const std::string& path, const JS::Bytes& data, double offset)
    {
      if (!data.instanceof(val::global("Uint8Array")))
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "writeFileChunk() takes the data as Uint8Array (" + path + ")");
      }
      const off_t at = byteCount(offset, "The offset", path);
      // (the length of a typed array can exceed 4 GiB, more than a size_t of wasm32 holds)
      const off_t length = byteCount(data["length"].as<double>(), "The length of the data", path);
      const off_t end = at + length;
      checkFileSize(end, path);
      const size_t n = static_cast<size_t>(length);
      FileDescriptor file;
      struct stat st;
      file.fd = ::open(path.c_str(), O_WRONLY | O_CREAT, 0666);
      if (file.fd < 0 || ::fstat(file.fd, &st) != 0)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::strerror(errno));
      }
      if (n == 0) return;
      if (end > st.st_size) checkMemoryFor(static_cast<size_t>(growthAllocation(st.st_size, end, path)));
      const std::unique_ptr<uint8_t[]> buffer = copyBuffer(n);
      for (size_t done = 0; done < n;)
      {
        const size_t part = std::min(n - done, COPY_BUFFER_SIZE);
        val(typed_memory_view(part, buffer.get())).call<void>("set", data.call<val>("subarray", static_cast<double>(done), static_cast<double>(done + part)));
        for (size_t written = 0; written < part;)
        {
          const ssize_t w = ::pwrite(file.fd, buffer.get() + written, part - written, at + static_cast<off_t>(done + written));
          if (w <= 0)
          {
            throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::strerror(w < 0 ? errno : EIO));
          }
          written += static_cast<size_t>(w);
        }
        done += part;
      }
    }

    /// Up to @p length bytes of the file @p path from byte @p offset, as a new Uint8Array (fewer at its end, none past it)
    JS::Bytes readFileChunk(const std::string& path, double offset, double length)
    {
      const off_t at = byteCount(offset, "The offset", path);
      const off_t wanted = byteCount(length, "The length", path);
      FileDescriptor file;
      struct stat st;
      file.fd = ::open(path.c_str(), O_RDONLY);
      if (file.fd < 0 || ::fstat(file.fd, &st) != 0)
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path);
      }
      const size_t n = at >= st.st_size ? 0 : static_cast<size_t>(std::min<off_t>(wanted, st.st_size - at));
      val result = val::global("Uint8Array").new_(n);
      const std::unique_ptr<uint8_t[]> buffer = copyBuffer(n);
      size_t read = 0;
      while (read < n)
      {
        const ssize_t r = ::pread(file.fd, buffer.get(), std::min(n - read, COPY_BUFFER_SIZE), at + static_cast<off_t>(read));
        if (r < 0) throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path + ": " + std::strerror(errno));
        if (r == 0) break; // shorter than its size said
        result.call<void>("set", val(typed_memory_view(static_cast<size_t>(r), buffer.get())), static_cast<double>(read));
        read += static_cast<size_t>(r);
      }
      return read == n ? result : result.call<val>("slice", 0, static_cast<double>(read));
    }

    /// Creates the file @p path with @p size bytes (zeros), or resizes it (keeping its contents up to that size), with
    /// one allocation: then writeFileChunk() writes into it without growing it.
    void allocateFile(const std::string& path, double size)
    {
      const off_t bytes = byteCount(size, "The size", path);
      checkFileSize(bytes, path);
      FileDescriptor file;
      struct stat st;
      file.fd = ::open(path.c_str(), O_WRONLY | O_CREAT, 0666);
      if (file.fd < 0 || ::fstat(file.fd, &st) != 0)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::strerror(errno));
      }
      if (bytes == st.st_size) return;
      if (bytes > st.st_size) checkMemoryFor(static_cast<size_t>(growthAllocation(st.st_size, bytes, path)));
      if (::ftruncate(file.fd, bytes) != 0)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::strerror(errno));
      }
    }

    /// Creates the file @p path anew with @p size bytes (zeros) in a block of @p capacity bytes (at least @p size): an
    /// existing file at @p path is replaced (its block is freed first; if it is open elsewhere, when that is closed).
    /// mimalloc does not give a freed block to a larger one: files that come one after the other, each removed before
    /// the next, keep the heap flat only if their blocks have the same size, e.g. the size of the largest of them.
    /// The file grows within its capacity (up to @p capacity bytes) without a new block, but writeFileChunk() checks a
    /// write past the size as growth of a block of that size (WasmFS does not tell the capacity of a file).
    void allocateFileWithCapacity(const std::string& path, double size, double capacity)
    {
      const off_t bytes = byteCount(size, "The size", path);
      const off_t reserved = byteCount(capacity, "The capacity", path);
      checkFileSize(reserved, path);
      if (reserved < bytes)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "The capacity has to be at least the size (" + std::to_string(bytes) + " bytes): " + path, std::to_string(reserved));
      }
      if (::unlink(path.c_str()) != 0 && errno != ENOENT)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::strerror(errno));
      }
      FileDescriptor file;
      file.fd = ::open(path.c_str(), O_WRONLY | O_CREAT, 0666);
      if (file.fd < 0)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::strerror(errno));
      }
      // the file is empty: its std::vector allocates exactly the capacity, and keeps it when it shrinks to the size
      checkMemoryFor(static_cast<size_t>(reserved));
      if (::ftruncate(file.fd, reserved) != 0 || ::ftruncate(file.fd, bytes) != 0)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::strerror(errno));
      }
    }

    /// Removes the mount point @p mount_point (a directory created by mountHostDirectory()).
    void unmountDirectory(const std::string& mount_point)
    {
      const int error = wasmfs_unmount(mount_point.c_str());
      if (error < 0)
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mount_point, std::strerror(-error));
      }
    }
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_misc)
{
  function("getVersion", OMS_FN(&version));
  function("getRevision", OMS_FN(&revision));
  function("getBranch", OMS_FN(&branch));
  function("getOpenMSDataPath", OMS_FN(&dataPath));
  function("setLogEnabled", OMS_FN(&setLogEnabled));
  function("isLogEnabled", OMS_FN(&isLogEnabled));
  function("getThreadingInfo", OMS_FN(&threadingInfo));
  function("getMemoryInfo", OMS_FN(&memoryInfo));
  function("hasONNXRuntime", OMS_FN(&hasONNXRuntime));
  function("mountHostDirectory", OMS_FN(&mountHostDirectory));
  function("unmount", OMS_FN(&unmountDirectory));
  function("writeFileChunk", OMS_FN(&writeFileChunk));
  function("readFileChunk", OMS_FN(&readFileChunk));
  function("allocateFile", OMS_FN(&allocateFile));
  function("allocateFile", OMS_FN(&allocateFileWithCapacity));
}
