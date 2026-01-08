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

// Forward declarations to access global LogStream instances
extern Logger::LogStream OpenMS_Log_fatal;
extern Logger::LogStream OpenMS_Log_error;
extern Logger::LogStream OpenMS_Log_warn;
extern Logger::LogStream OpenMS_Log_info;
extern Logger::LogStream OpenMS_Log_debug;

namespace
{
  /// Internal structure holding thread-local log streams (implementation detail)
  struct ThreadLogStreamsImpl
  {
    /// Each thread gets its own set of LogStream instances with thread-local buffers
    /// but shared stream_list_ with the global LogStream instances
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
      // Create LogStreamBuf with thread-local buffers but shared stream_list_
      // This way each thread has its own buffer/cache (for thread safety)
      // but shares the output stream configuration with the global LogStreams
      fatal = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf(OpenMS_Log_fatal.rdbuf(), &red), true);
      error = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf(OpenMS_Log_error.rdbuf(), &red), true);
      warn = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf(OpenMS_Log_warn.rdbuf(), &yellow), true);
      info = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf(OpenMS_Log_info.rdbuf(), nullptr), true);
      debug = std::make_unique<Logger::LogStream>(
          new Logger::LogStreamBuf(OpenMS_Log_debug.rdbuf(), &magenta), true);
      initialized = true;
    }

    void flushAll()
    {
      // Use flushIncomplete() to ensure incomplete lines (text not terminated by newline)
      // are also flushed to the streams
      if (fatal) fatal->flushIncomplete();
      if (error) error->flushIncomplete();
      if (warn) warn->flushIncomplete();
      if (info) info->flushIncomplete();
      if (debug) debug->flushIncomplete();
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
