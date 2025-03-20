// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>

namespace OpenMS
{

  namespace Exception
  {

/**
  @brief OpenMS global exception handler

  @ingroup Exceptions
*/
    class OPENMS_DLLAPI GlobalExceptionHandler
    {
private:
      /**	@name	Constructors
      */
      //@{

      /**	@brief Default constructor.

          This constructor installs the OPENMS specific handlers for
          <tt>terminate</tt>, <tt>unexpected</tt>, and <tt>new_handler</tt>.
          <tt>terminate</tt> or <tt>unexpected</tt> are called to abort a
          program if an exception was not caught or a function exits via an
          exception that is not allowed by its exception specification. Both
          functions are replaced by a function of GlobalExceptionHandler that
          tries to determine the last exception thrown. This mechanism only
          works, if all exceptions are derived from Base.

          The default <tt>new_handler</tt> is replaced by #newHandler and
          throws an exception of type OutOfMemory instead of
          <tt>bad_alloc</tt> (the default behaviour defined in the ANSI C++
          standard).
      */
      GlobalExceptionHandler()
      throw();

      //@}

      /// private version of c'tor to avoid initialization
      GlobalExceptionHandler(const GlobalExceptionHandler &) {}

      ~GlobalExceptionHandler() {}
public:

      /// The accessor for the singleton. It also serves as a replacement for the constructor.
      static GlobalExceptionHandler & getInstance()
      {
        static GlobalExceptionHandler * globalExceptionHandler_;

        if (globalExceptionHandler_ == nullptr)
        {
          globalExceptionHandler_ = new GlobalExceptionHandler;
        }
        return *globalExceptionHandler_;
      }

      /**	@name	Accessors
      */
      //@{

      /**
      */
      static void setName(const std::string & name)
      throw();

      /**
      */
      static void setMessage(const std::string & message)
      throw();

      /**
      */
      static void setLine(int line)
      throw();

      /**
      */
      static void setFile(const std::string & file)
      throw();

      /**
      */
      static void setFunction(const std::string & function)
      throw();

      /**
      */
      static void set(const std::string & file, int line, const std::string & function,
                      const std::string & name, const std::string & message)
      throw();
      //@}
protected:

      /// The OPENMS replacement for terminate
      static void terminate()
      throw();

      /// The OPENMS new handler
      static void newHandler();

      /**
        @name	Wrapper for static member.

        @note To avoid problems when accessing uninitialised static
        members we replaced them with static functions returning
        references to the members.
       */
      //@{

      /// wrapper for static member file_
      static std::string & file_()
      {
        static std::string * file_ = nullptr;
        if (file_ == nullptr)
        {
          file_  = new std::string;
          *file_ = "unknown";
        }
        return *file_;
      }

      /// wrapper for static member line_
      static int & line_()
      {
        static int * line_ = nullptr;
        if (line_ == nullptr)
        {
          line_  = new int;
          *line_ = -1;
        }
        return *line_;
      }

      /// wrapper for static member function_
      static std::string & function_()
      {
        static std::string * function_ = nullptr;
        if (function_ == nullptr)
        {
          function_  = new std::string;
          *function_ = "unknown";
        }
        return *function_;
      }

      /// wrapper for static member name_
      static std::string & name_()
      {
        static std::string * name_ = nullptr;
        if (name_ == nullptr)
        {
          name_  = new std::string;
          *name_ = "unknown exception";
        }
        return *name_;
      }

      /// wrapper for static member what_
      static std::string & what_()
      {
        static std::string * what_ = nullptr;
        if (what_ == nullptr)
        {
          what_  = new std::string;
          *what_ = " - ";
        }
        return *what_;
      }

      //@}
    };

  } // namespace Exception

} // namespace OpenMS

