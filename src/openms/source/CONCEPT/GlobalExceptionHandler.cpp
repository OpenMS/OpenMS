// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <cstdlib>  // for getenv in terminate()
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
        std::cout << "exception of type " << name_().c_str() << " occured in line "
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

  } // namespace Exception

} // namespace OpenMS
