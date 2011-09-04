// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_EXCEPTION_H
#define OPENMS_CONCEPT_EXCEPTION_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>

#include <new>
#include <string>
#include <sstream>

namespace OpenMS 
{

	/**
		@defgroup Exceptions Exceptions
		
		@brief Exceptions
		
		@ingroup Concept
	*/


	/** 
		@brief %Exception namespace
		
		@ingroup Concept
	 */
	namespace Exception 
	{
		
		/**	
			@brief Exception base class.

		  This class is intended as a base class for all other exceptions.
		
			Each exception class should define a constructor taking a filename (string), line (int) and function name (string)
			as first agruments. This information is usually printed in case of an uncaught exception.
		  
			To support this feature, each @em throw directive should look as follows:
			@code throw Exception::Exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,...); @endcode
			
			@em __FILE__ and @em __LINE__ are built-in preprocessor macros that hold the desired information.
			@n @em __PRETTY_FUNCTION__ is replaced by the GNU G++ compiler with the demangled name of the current function.
			(For other compilers it is defined as "<unknown>" in config.h.)

			%OpenMS provides its own Exception::GlobalExceptionHandler::terminate() handler. This handler extracts as much 
		  information as possible from the exception, prints it to @em cerr , and finally calls exits the program
			cleanly (with exit code 1). This can be rather inconvenient for debugging, since you are told where the
			exception was thrown, but in general you do not know anything about the context.  Therefore
			terminate() can also create a core dump. Using a debugger (e.g. @em dbx or @em gdb) you can then create a 
			stack traceback. To create a core dump, you should set the environment variable @em OPENMS_DUMP_CORE to any
			(non empty) value.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI BaseException
			:	public std::exception
		{
			public:

				/**	@name	Constructors and Destructors
				*/
				//@{

				/// Default constructor
				BaseException() throw();
				
				/// Constructor
				BaseException(const char* file, int line, const char* function) throw();

				/// Constructor
				BaseException(const char* file, int line, const char* function, const std::string& name , const std::string& message) throw();

				/// Copy constructor
				BaseException(const BaseException& exception) throw();

				/// Destructor
				virtual ~BaseException() throw();
				//@}

				/**	@name	Accessors
				*/
				//@{
		
				///	Returns the name of the exception 
				const char* getName() const throw();

				///	Returns the error message of the exception
				virtual const char* what() const throw();

				/// Returns the line number where it occured
				int getLine() const throw();
		
				/// Returns the file where it occured
				const char* getFile() const throw();

				/// Returns the function where it occured
				const char* getFunction() const throw();

				/// Returns the message
				const char* getMessage() const throw();

				/// Modify the exception's error message
				void setMessage(const std::string& message) throw();

				//@}

			protected:
				
				/// The source file the exception was thrown in
				const char*	file_;

				/// The line number the exception was thrown in
				int					line_;

				/// The source file the exception was thrown in
				const char*	function_;
					
				/// The name of the exception.
				std::string name_;
				
				/// A more detailed description of the exception's cause.
				std::string what_;
		};		

		/**	
			@brief Precondition failed exception.
			
			A precondition (as defined by @ref OPENMS_PRECONDITION ) has failed.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI Precondition
			: public BaseException
		{
			public:
			  Precondition(const char* file, int line, const char* function, const std::string& condition)	throw();
		};

		/**	
			@brief Postcondition failed exception.
			
			A postcondition (as defined by @ref OPENMS_POSTCONDITION ) has failed.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI Postcondition
			: public BaseException
		{
			public:
			  Postcondition(const char* file, int line, const char* function, const std::string& condition) throw();
		};

		/**	
			@brief Not all required information provided.
			
			Information that are required are not provided.
			Especially usefull for missing MetaInfo values.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI MissingInformation
			: public BaseException
		{
			public:
				MissingInformation(const char* file, int line, const char* function, const std::string& error_message) throw();
		};


		/**	
			@brief Int underflow exception.
				
			Throw this exception to indicate an index that was smaller than
			allowed.  The constructor has two additional arguments, the values
			of which should be set to the index that caused the failure and the
			smallest allowed value to simplify debugging.

			@param	index the value of the index causing the problem
			@param	size	smallest value allowed for index

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI IndexUnderflow 
			: public BaseException
		{
			public:
				IndexUnderflow(const char* file, int line, const char* function, SignedSize index = 0, Size size = 0) throw();
		};

		/**	
			@brief UInt underflow exception.

			Throw this exception to indicate a size was smaller than allowed.
			The constructor has an additional argument: the value of of the
			requested size.  This exception is thrown, if buffer sizes are
			insufficient.

			@param	size the size causing the problem

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI SizeUnderflow 
			: public BaseException
		{
			public:
			  SizeUnderflow(const char* file, int line, const char* function, Size size = 0) throw();
		};

		/**	
			@brief Int overflow exception.
			
			Throw this exception to indicate an index that was larger than
			allowed.  The constructor has two additional arguments, the values
			of which should be set to the index that caused the failure and the
			largest allowed value to simplify debugging.
			@param	index the value of the index causing the problem
			@param	size	largest value allowed for index

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI IndexOverflow 
			: public BaseException
		{
			public:
				IndexOverflow(const char* file, int line, const char* function, SignedSize index = 0, Size size = 0) throw();
		};

		/**	
			@brief A call to an external library (other than OpenMS) went wrong.

			Throw this exception to indicate that an external library call came
			back unsuccessfull.

			@param	size the size causing the problem

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI FailedAPICall 
			: public BaseException
		{
			public:
			  FailedAPICall(const char* file, int line, const char* function, const std::string& message) throw();
		};
		
		/**	
			@brief Invalid range exception.

			Use this exception to indicate a general range problems.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI InvalidRange 
			: public BaseException
		{
			public:
				InvalidRange(const char* file, int line, const char* function) throw();
		};


		/**	
			@brief Invalid UInt exception.
				
			Throw this exception to indicate that a size was unexpected.
			The constructor has an additional argument: the value of of the
			requested size. 
			@param	size the size causing the problem

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI InvalidSize 
			: public BaseException
		{
			public:
				InvalidSize(const char* file, int line, const char* function, Size size = 0) throw();
		};


		/**	
			@brief Out of range exception.
				
			Use this exception to indicate that a given value is out of a
			defined range, i. e. not within the domain of a function.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI OutOfRange 
			: public BaseException
		{
			public:
				OutOfRange(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief Invalid value exception.
				
			Use this exception to indicate that a given value is not valid,
			when the value is only allowed to be out of a certain set of values.
			If the value has to be inside a given range, you should rather use OutOfRange.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI InvalidValue 
			: public BaseException
		{
			public:
				InvalidValue(const char* file, int line, const char* function, const std::string& message ,const std::string& value) throw();
		};

		/**	
			@brief Exception indicating that an invalid parameter was handed over to an algorithm.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI InvalidParameter 
			: public BaseException
		{
			public:
				InvalidParameter(const char* file, int line, const char* function, const std::string& message) throw();
		};

		/**	
			@brief Invalid conversion exception.
		
			This exception indicates a conversion problem when converting from
			one type to another.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI ConversionError 
			: public BaseException
		{
			public:
				ConversionError(const char* file, int line, const char* function, const std::string& error) throw();
		};

		/**	
			@brief Illegal self operation exception.	
				
			Throw this exception to indicate an invalid operation on the object
			itself. In general these operations are self assignments or related
			methods.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI IllegalSelfOperation
			: public BaseException
		{
			public:
			  IllegalSelfOperation(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief Null pointer argument is invalid exception.
				
			Use this exception to indicate a failure due to an argument not
			containing a pointer to a valid object, but a null pointer.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI NullPointer 
			: public BaseException
		{
			public:
			  NullPointer(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief Invalid iterator exception.
				
			The iterator on which an operation should be performed was invalid.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI InvalidIterator
			: public BaseException
		{
			public:
			  InvalidIterator(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief Incompatible iterator exception.
				
			The iterators could not be assigned because they are bound to
			different containers.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI IncompatibleIterators
			: public BaseException
		{
			public:
			  IncompatibleIterators(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief Not implemented exception.
				
			This exception should be thrown to indicate not yet implemented methods.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI NotImplemented
			: public BaseException
		{
			public:
			  NotImplemented(const char* file, int line, const char* function) throw();
		};

		/** 
			@brief Illegal tree operation exception.

			This exception is thrown to indicate that an illegal tree operation was requested.
			
			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI IllegalTreeOperation
			: public BaseException
		{
			public:
			  IllegalTreeOperation(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief Out of memory exception.
				
			Throw this exception to indicate that an allocation failed.
			This exception is thrown in the OPENMS new handler.
			@param	size	the number of bytes that should have been allocated
			@see GlobalException::newHandler

			@ingroup Exceptions
		*/
#ifdef _MSC_VER // disable some seqan warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4275 )
#endif
		class OPENMS_DLLAPI OutOfMemory
			: public BaseException, public std::bad_alloc
		{
			public:
				OutOfMemory(const char* file, int line, const char* function, Size size = 0) throw();
		};
#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#endif
		/**	
			@brief Buffer overflow exception.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI BufferOverflow 
			: public BaseException
		{
			public:
			  BufferOverflow(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief Division by zero error exception.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI DivisionByZero 
			: public BaseException
		{
			public:
			  DivisionByZero(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief Out of grid exception.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI OutOfGrid 
			: public BaseException
		{
			public:
			  OutOfGrid(const char* file, int line, const char* function) throw();
		};

		/**	
			@brief File not found exception.

			A given file could not be found.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI FileNotFound 
			: public BaseException
		{
			public:
				FileNotFound(const char* file, int line, const char* function, const std::string& filename) throw();
		};

		/**	
			@brief File not readable exception.

			A given file is not readable for the current user.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI FileNotReadable 
			: public BaseException
		{
			public:
				FileNotReadable(const char* file, int line, const char* function, const std::string& filename) throw();
		};

		/** 
		  @brief File not writable exception.

			A given file is not writable for the current user.

			@ingroup Exceptions
		*/
	 	class OPENMS_DLLAPI FileNotWritable
			: public BaseException
		{
			public:
				FileNotWritable(const char* file, int line, const char* function, const std::string& filename) throw();
		};

    /**
      @brief General IOException.

      General error for IO operations, that can not be associated to the more specific exceptions (e.g. FileNotWritable)

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI IOException
      : public BaseException
    {
      public:
        IOException(const char* file, int line, const char* function, const std::string& filename) throw();
    };



		/**	
			@brief File is empty.

			A given file is empty.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI FileEmpty 
			: public BaseException
		{
			public:
				FileEmpty(const char* file, int line, const char* function, const std::string& filename) throw();
		};

		/**	
			@brief Invalid 3-dimensional position exception.

			A given position in three dimensional is invalid.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI IllegalPosition 
			: public BaseException
		{
			public:
			  IllegalPosition(const char* file, int line, const char* function, float x, float y, float z) throw();
		};

		/**	
			@brief Parse Error exception.

			A given expression could not be parsed.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI ParseError
			: public BaseException
		{
			public:
			  ParseError(const char* file, int line, const char* function, const std::string& expression, const std::string& message) throw();
		};

		/**	
			@brief Unable to create file exception.

			The given file could not be created.

			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI UnableToCreateFile
			: public BaseException
		{
			public:
				UnableToCreateFile(const char* file, int line, const char* function, const std::string& filename) throw();
		};
		
    /**
			@brief A method or algorithm argument contains illegal values
			
			@ingroup Exceptions
		*/
    class OPENMS_DLLAPI IllegalArgument 
			: public BaseException
    {
			public:
				IllegalArgument(const char* file, int line, const char* function, const std::string& error_message) throw();
    };

		/**	
			@brief Element could not be found exception.

			The given element could not be found. 

			@ingroup Exceptions
		*/			
		class OPENMS_DLLAPI ElementNotFound
			: public BaseException
		{
			public:
				ElementNotFound(const char* file, int line, const char* function, const std::string& element)	throw();
		};

		/**	
			@brief Exception used if an error occurred while fitting a model to a given dataset
			
			The given element could not be found. 

			@ingroup Exceptions
		*/			
		class OPENMS_DLLAPI UnableToFit
			: public BaseException
		{
			public:
				UnableToFit(const char* file, int line, const char* function, const std::string& name , const std::string& message) throw();
		};

		/**	
			@brief Exception used if an error occurred while calibrating a dataset.
				
			The calibration can not be performed because not enough reference masses
					were detected.
				
			@ingroup Exceptions
		*/			
		class OPENMS_DLLAPI UnableToCalibrate
			: public BaseException
		{
		public:
			UnableToCalibrate(const char* file, int line, const char* function, const std::string& name , const std::string& message) throw();
		};

		/**	
			@brief Exception used if no more unique document ID's can be drawn from ID pool.
				
			The ID pool of OpenMS is either depleted or not existant.
				
			@ingroup Exceptions
		*/			
		class OPENMS_DLLAPI DepletedIDPool
			: public BaseException
		{
		public:
			DepletedIDPool(const char* file, int line, const char* function, const std::string& name , const std::string& message) throw();
		};


		/**
			@brief OpenMS global exception handler
		
			@ingroup Exceptions
		*/
		class OPENMS_DLLAPI GlobalExceptionHandler
		{
			public:
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
				static std::string file_;
				static int				 line_;
				static std::string function_;
				static std::string name_;
				static std::string what_;
		};

		///Global static instance of GlobalExceptionHandler
		extern GlobalExceptionHandler globalHandler;

		}
		
		
		/**	
			@brief Output operator for exceptions.
		
			All %OpenMS exceptions can be printed to an arbitrary output stream.
			Information written contains the exception class, the error message,
			and the location (file, line number). The following code block
			can thus be used to catch any %OpenMS exceptions and convert them to
			human readable information:
			\code
			try
			{
				.... // some code which potentially throws an exception
			}
			catch (Exception::Exception e)
			{
				Log.error() << "caught exception: " << e << std::endl;
			}
			\endcode
			
			@ingroup Exceptions
		*/
		OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Exception::BaseException& e);
	
} // namespace OPENMS

#endif // OPENMS_CONCEPT_EXCEPTION_H
