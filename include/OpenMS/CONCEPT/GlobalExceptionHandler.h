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

#ifndef OPENMS_CONCEPT_GLOBALEXCEPTIONHANDLER_H
#define OPENMS_CONCEPT_GLOBALEXCEPTIONHANDLER_H

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
    GlobalExceptionHandler(const GlobalExceptionHandler&) {}

    ~GlobalExceptionHandler() {}
  public:

    /// The accessor for the singleton. It also serves as a replacement for the constructor.
    static GlobalExceptionHandler& getInstance()
    {
      static GlobalExceptionHandler* globalExceptionHandler_;

      if (globalExceptionHandler_ == 0)
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
    static void setName(const std::string& name)
      throw();

    /**
    */
    static void setMessage(const std::string& message)
      throw();

    /**
    */
    static void setLine(int line)
      throw();

    /**
    */
    static void setFile(const std::string& file)
      throw();

    /**
    */
    static void setFunction(const std::string& function)
      throw();

    /**
    */
    static void set
      (const std::string& file, int line, const std::string& function,
       const std::string& name, const std::string& message)
      throw();
    //@}
  protected:

    /// The OPENMS replacement for terminate
    static void terminate()
      throw();

    /// The OPENMS new handler
#ifdef OPENMS_COMPILER_MSVC
    static void newHandler();
#else
    static void newHandler() throw(OutOfMemory);
#endif

    /**
      @name	Wrapper for static member.

      @note To avoid problems when accessing uninitialised static
      members we replaced them with static functions returning
      references to the members.
     */
    //@{

    /// wrapper for static member file_
    static std::string& file_()
    {
      static std::string * file_ = 0;
      if(file_ == 0)
      {
        file_  = new std::string;
        *file_ = "unknown";
      }
      return *file_;
    }

    /// wrapper for static member line_
    static int& line_()
    {
      static int * line_ = 0;
      if(line_ == 0)
      {
        line_  = new int;
        *line_ = -1;
      }
      return *line_;
    }

    /// wrapper for static member function_
    static std::string& function_()
    {
      static std::string * function_ = 0;
      if(function_ == 0)
      {
        function_  = new std::string;
        *function_ = "unknown";
      }
      return *function_;
    }

    /// wrapper for static member name_
    static std::string& name_()
    {
      static std::string * name_ = 0;
      if(name_ == 0)
      {
        name_  = new std::string;
        *name_ = "unknown exception";
      }
      return *name_;
    }

    /// wrapper for static member what_
    static std::string& what_()
    {
      static std::string * what_ = 0;
      if(what_ == 0)
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

#endif // OPENMS_CONCEPT_GLOBALEXCEPTIONHANDLER_H
