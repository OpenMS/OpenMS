// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <cstdlib>  // for getenv in terminate()
#include <exception>
//#include <sys/types.h>
#include <csignal> // for SIGSEGV and kill
#include <iostream>
#include <new>

#ifndef OPENMS_WINDOWSPLATFORM
  #ifdef OPENMS_HAS_UNISTD_H
  #include <unistd.h> // for getpid
  #endif
#endif

#define OPENMS_CORE_DUMP_ENVNAME "OPENMS_DUMP_CORE"

namespace OpenMS::Exception
{

    GlobalExceptionHandler::GlobalExceptionHandler() throw()
    {
      std::set_terminate(terminate);
      //std::set_unexpected(terminate); // removed in c++17
      std::set_new_handler(newHandler);
    }

    void GlobalExceptionHandler::newHandler()
    {
      throw OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    void GlobalExceptionHandler::terminate() throw()
    {
      // add cerr to the log stream
      // and write all available information on
      // the exception to the log stream (potentially with an assigned file!)
      // and cerr

      std::cout << std::endl;
      std::cout << "---------------------------------------------------" << std::endl;
      std::cout << "FATAL: uncaught exception!" << std::endl;
      std::cout << "---------------------------------------------------" << std::endl;
      // The exception that was not caught reports itself: it may have been constructed on another thread (and carried
      // here by std::exception_ptr, e.g. from a parallel loop), so the last entry of this thread need not be about it.
      bool reported = false;
      if (const std::exception_ptr uncaught = std::current_exception())
      {
        try
        {
          std::rethrow_exception(uncaught);
        }
        catch (const BaseException& e)
        {
          std::cout << "exception of type " << e.getName() << " occurred in line " << e.getLine() << ", function "
                    << e.getFunction() << " of " << e.getFile() << std::endl;
          std::cout << "error message: " << e.what() << std::endl;
          reported = true;
        }
        catch (...)
        {
        }
      }
      if (!reported && (line_() != -1) && (name_() != "unknown"))
      {
        std::cout << "last entry in the exception handler: " << std::endl;
        std::cout << "exception of type " << name_().c_str() << " occurred in line "
                  << line_() << ", function " << function_() << " of " << file_().c_str() << std::endl;
        std::cout << "error message: " << what_().c_str() << std::endl;
      }
      std::cout << "---------------------------------------------------" << std::endl;

#ifndef OPENMS_WINDOWSPLATFORM
      // if the environment variable declared in OPENMS_CORE_DUMP_ENVNAME
      // is set, provoke a core dump (this is helpful to get a stack traceback)
      if (getenv(OPENMS_CORE_DUMP_ENVNAME) != nullptr)
      {
#ifdef OPENMS_HAS_KILL
        std::cout << "dumping core file.... (to avoid this, unset " << OPENMS_CORE_DUMP_ENVNAME
                  << " in your environment)" << std::endl;
        // provoke a core dump
        kill(getpid(), SIGSEGV);
#endif
      }
#endif

      // otherwise exit as default terminate() would:
      abort();
    }

    void GlobalExceptionHandler::set(const std::string & file, int line, const std::string & function, const std::string & name, const std::string & message) throw()
    {
      GlobalExceptionHandler::name_() = name;
      GlobalExceptionHandler::line_() = line;
      GlobalExceptionHandler::what_() = message;
      GlobalExceptionHandler::file_() = file;
      GlobalExceptionHandler::function_() = function;
    }

    void GlobalExceptionHandler::setName(const std::string & name) throw()
    {
      GlobalExceptionHandler::name_() = name;
    }

    void GlobalExceptionHandler::setMessage(const std::string & message) throw()
    {
      GlobalExceptionHandler::what_() = message;
    }

    void GlobalExceptionHandler::setFile(const std::string & file) throw()
    {
      GlobalExceptionHandler::file_() = file;
    }

    void GlobalExceptionHandler::setFunction(const std::string & function) throw()
    {
      GlobalExceptionHandler::function_() = function;
    }

    void GlobalExceptionHandler::setLine(int line) throw()
    {
      GlobalExceptionHandler::line_() = line;
    }

    GlobalExceptionHandler & GlobalExceptionHandler::getInstance()
    {
      // A function-local static is initialised once, also when several threads construct their first
      // exception at the same time (the former check-then-new on a static pointer was a data race).
      static GlobalExceptionHandler globalExceptionHandler_;
      return globalExceptionHandler_;
    }

    // The last exception is recorded per thread: every exception constructor writes these fields,
    // and threads that throw at the same time (e.g. the chunks of a parallel file read that fail to
    // parse) would otherwise write the same strings at once. terminate() runs on the thread whose
    // exception was not caught, so it reads that thread's entry.
    std::string & GlobalExceptionHandler::file_()
    {
      static thread_local std::string file = "unknown";
      return file;
    }

    int & GlobalExceptionHandler::line_()
    {
      static thread_local int line = -1;
      return line;
    }

    std::string & GlobalExceptionHandler::function_()
    {
      static thread_local std::string function = "unknown";
      return function;
    }

    std::string & GlobalExceptionHandler::name_()
    {
      static thread_local std::string name = "unknown exception";
      return name;
    }

    std::string & GlobalExceptionHandler::what_()
    {
      static thread_local std::string what = " - ";
      return what;
    }


} // namespace OpenMS::Exception
