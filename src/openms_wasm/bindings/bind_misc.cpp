// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Version, shared data and logging.

#include "binding_utils.h"

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>

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
}
