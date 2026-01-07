// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ThreadLogContext.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Colorizer.h>

#include <iostream>
#include <memory>

namespace OpenMS
{

namespace
{
  /// Internal structure holding thread-local log streams (implementation detail)
  struct ThreadLogStreamsImpl
  {
    /// Each thread gets its own set of LogStream instances
    std::unique_ptr<Logger::LogStream> fatal;
    std::unique_ptr<Logger::LogStream> error;
    std::unique_ptr<Logger::LogStream> warn;
    std::unique_ptr<Logger::LogStream> info;
    std::unique_ptr<Logger::LogStream> debug;
    bool initialized = false;

    ThreadLogStreamsImpl()
    {
      initialize();
    }

    ~ThreadLogStreamsImpl()
    {
      // Flush all streams before destruction
      flushAll();
    }

    void initialize()
    {
      // Create LogStreamBuf with same configuration as the global streams
      // Each thread gets its own buffer, cache, and incomplete_line_ state
      fatal = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf("FATAL_ERROR", &red), true, &std::cerr);
      error = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf("ERROR", &red), true, &std::cerr);
      warn = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf("WARNING", &yellow), true, &std::cout);
      info = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf("INFO", nullptr), true, &std::cout);
      // Debug is disabled by default (no stream attached), same as global
      debug = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf("DEBUG", &magenta), false);
      initialized = true;
    }

    void flushAll()
    {
      if (fatal) fatal->flush();
      if (error) error->flush();
      if (warn) warn->flush();
      if (info) info->flush();
      if (debug) debug->flush();
    }

    void reset()
    {
      flushAll();
      initialized = false;
      fatal.reset();
      error.reset();
      warn.reset();
      info.reset();
      debug.reset();
      initialize();
    }
  };

  // Thread-local storage for log streams
  // Each thread gets its own instance, automatically constructed on first access
  // and destroyed when the thread exits
  thread_local std::unique_ptr<ThreadLogStreamsImpl> tls_streams;

  ThreadLogStreamsImpl& getThreadLogStreams()
  {
    if (!tls_streams)
    {
      tls_streams = std::make_unique<ThreadLogStreamsImpl>();
    }
    return *tls_streams;
  }

} // anonymous namespace

  Logger::LogStream& ThreadLogContext::fatal()
  {
    return *getThreadLogStreams().fatal;
  }

  Logger::LogStream& ThreadLogContext::error()
  {
    return *getThreadLogStreams().error;
  }

  Logger::LogStream& ThreadLogContext::warn()
  {
    return *getThreadLogStreams().warn;
  }

  Logger::LogStream& ThreadLogContext::info()
  {
    return *getThreadLogStreams().info;
  }

  Logger::LogStream& ThreadLogContext::debug()
  {
    return *getThreadLogStreams().debug;
  }

  bool ThreadLogContext::isInitialized()
  {
    return tls_streams && tls_streams->initialized;
  }

  void ThreadLogContext::ensureInitialized()
  {
    (void)getThreadLogStreams(); // Force initialization
  }

  void ThreadLogContext::reset()
  {
    if (tls_streams)
    {
      tls_streams->reset();
    }
    else
    {
      ensureInitialized();
    }
  }

  void ThreadLogContext::flushAll()
  {
    if (tls_streams)
    {
      tls_streams->flushAll();
    }
  }

} // namespace OpenMS
