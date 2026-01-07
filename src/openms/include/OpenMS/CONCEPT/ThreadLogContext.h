// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <memory>

namespace OpenMS
{
  // Forward declarations
  namespace Logger
  {
    class LogStream;
  }

  /**
    @brief Thread-local logging context for thread-safe logging.

    This class provides thread-local LogStream instances to eliminate data races
    in the logging system. Each thread gets its own set of log streams (fatal, error,
    warn, info, debug), ensuring that concurrent logging from multiple threads
    does not cause race conditions on shared buffers, caches, or incomplete line state.

    The implementation uses C++11 thread_local storage, ensuring automatic initialization
    on first access and automatic cleanup when threads exit.

    @see GitHub Issue #8596 - Thread Safety Bug in LogStream

    Usage:
    @code
    // Access the thread-local info logger
    ThreadLogContext::info() << "This is thread-safe" << std::endl;

    // Or via the thread-safe macros (when OPENMS_THREADLOCAL_LOGGING is enabled)
    OPENMS_LOG_INFO << "Also thread-safe" << std::endl;
    @endcode

    @ingroup Concept
  */
  class OPENMS_DLLAPI ThreadLogContext
  {
  public:
    /// Access the thread-local fatal error log stream
    static Logger::LogStream& fatal();

    /// Access the thread-local error log stream
    static Logger::LogStream& error();

    /// Access the thread-local warning log stream
    static Logger::LogStream& warn();

    /// Access the thread-local info log stream
    static Logger::LogStream& info();

    /// Access the thread-local debug log stream
    static Logger::LogStream& debug();

    /// Check if the current thread's log context has been initialized
    static bool isInitialized();

    /// Force initialization of the current thread's log context.
    /// This is normally called automatically on first access.
    static void ensureInitialized();

    /// Reset the current thread's log context to default configuration.
    /// This flushes all streams and reinitializes them.
    static void reset();

    /// Flush all log streams for the current thread
    static void flushAll();
  };

} // namespace OpenMS
