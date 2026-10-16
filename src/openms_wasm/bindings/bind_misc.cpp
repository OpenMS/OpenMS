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

#include <cerrno>
#include <cstring>
#include <iostream>

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
}
