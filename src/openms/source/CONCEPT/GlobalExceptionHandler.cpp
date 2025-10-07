// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <cstdlib>  // for getenv in terminate()
//#include <sys/types.h>
#include <csignal> // for SIGSEGV and kill
#include <iostream>

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
      if ((line_() != -1) && (name_() != "unknown"))
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
      static GlobalExceptionHandler * globalExceptionHandler_;

      if (globalExceptionHandler_ == nullptr)
      {
        globalExceptionHandler_ = new GlobalExceptionHandler;
      }
      return *globalExceptionHandler_;
    }

    std::string & GlobalExceptionHandler::file_()
    {
      static std::string * file_ = nullptr;
      if (file_ == nullptr)
      {
        file_  = new std::string;
        *file_ = "unknown";
      }
      return *file_;
    }

    int & GlobalExceptionHandler::line_()
    {
      static int * line_ = nullptr;
      if (line_ == nullptr)
      {
        line_  = new int;
        *line_ = -1;
      }
      return *line_;
    }

    std::string & GlobalExceptionHandler::function_()
    {
      static std::string * function_ = nullptr;
      if (function_ == nullptr)
      {
        function_  = new std::string;
        *function_ = "unknown";
      }
      return *function_;
    }

    std::string & GlobalExceptionHandler::name_()
    {
      static std::string * name_ = nullptr;
      if (name_ == nullptr)
      {
        name_  = new std::string;
        *name_ = "unknown exception";
      }
      return *name_;
    }

    std::string & GlobalExceptionHandler::what_()
    {
      static std::string * what_ = nullptr;
      if (what_ == nullptr)
      {
        what_  = new std::string;
        *what_ = " - ";
      }
      return *what_;
    }


} // namespace OpenMS::Exception
