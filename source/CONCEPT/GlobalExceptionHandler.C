// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <cstdlib>	// for getenv in terminate()
//#include <sys/types.h>
#include <csignal> // for SIGSEGV and kill

#ifndef OPENMS_WINDOWSPLATFORM
#ifdef OPENMS_HAS_UNISTD_H
#include <unistd.h> // for getpid
#endif
#ifdef OPENMS_HAS_PROCESS_H
#include <process.h>
#endif
#endif

#define OPENMS_CORE_DUMP_ENVNAME "OPENMS_DUMP_CORE"

namespace OpenMS
{

namespace Exception
{

GlobalExceptionHandler::GlobalExceptionHandler() throw()
{
  std::set_terminate(terminate);
  std::set_unexpected(terminate);
  std::set_new_handler(newHandler);
}

void GlobalExceptionHandler::newHandler()
#ifndef OPENMS_COMPILER_MSVC
  throw(OutOfMemory)
#endif
{
  throw OutOfMemory(__FILE__, __LINE__, __PRETTY_FUNCTION__);
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
    std::cout << "exception of type " << name_().c_str() << " occured in line "
              << line_() <<", function " << function_() << " of " << file_().c_str() << std::endl;
    std::cout << "error message: " << what_().c_str() << std::endl;
  }
  std::cout << "---------------------------------------------------" << std::endl;

#ifndef OPENMS_WINDOWSPLATFORM
  // if the environment variable declared in OPENMS_CORE_DUMP_ENVNAME
  // is set, provoke a core dump (this is helpful to get a stack traceback)
  if (getenv(OPENMS_CORE_DUMP_ENVNAME) != 0)
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

void GlobalExceptionHandler::set(const std::string& file, int line, const std::string& function, const std::string& name, const std::string& message) throw()
{
  GlobalExceptionHandler::name_() = name;
  GlobalExceptionHandler::line_() = line;
  GlobalExceptionHandler::what_() = message;
  GlobalExceptionHandler::file_() = file;
  GlobalExceptionHandler::function_() = function;
}

void GlobalExceptionHandler::setName(const std::string& name) throw()
{
  GlobalExceptionHandler::name_() = name;
}

void GlobalExceptionHandler::setMessage(const std::string& message) throw()
{
  GlobalExceptionHandler::what_() = message;
}

void GlobalExceptionHandler::setFile(const std::string& file) throw()
{
  GlobalExceptionHandler::file_() = file;
}

void GlobalExceptionHandler::setFunction(const std::string& function) throw()
{
  GlobalExceptionHandler::function_() = function;
}

void GlobalExceptionHandler::setLine(int line) throw()
{
  GlobalExceptionHandler::line_() = line;
}

} // namespace Exception

} // namespace OpenMS
