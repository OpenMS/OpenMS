// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <iosfwd>
#include <stdexcept>
#include <string>

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
      as first arguments. This information is usually printed in case of an uncaught exception.

      To support this feature, each @em throw directive should look as follows:
      @code throw Exception::Exception(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,...); @endcode

      @em __FILE__ and @em __LINE__ are built-in preprocessor macros that hold the desired information.
      @n @em OPENMS_PRETTY_FUNCTION is replaced by the GNU G++ compiler with the demangled name of the current function.
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
    class OPENMS_DLLAPI BaseException : public std::runtime_error
    {
    public:
      /**	@name	Constructors and Destructors
       */
      //@{

      /// Default constructor
      BaseException() noexcept;

      /// Constructor
      BaseException(const char* file, int line, const char* function) noexcept;

      /// Constructor
      BaseException(const char* file, int line, const char* function, const std::string& name, const std::string& message) noexcept;

      /// Copy constructor
      BaseException(const BaseException& exception) noexcept;

      /// Destructor
      ~BaseException() noexcept override;
      //@}

      /**	@name	Accessors
       */
      //@{

      /// Returns the name of the exception
      const char* getName() const noexcept;

      /// Returns the line number where it occurred
      int getLine() const noexcept;

      /// Returns the file where it occurred
      const char* getFile() const noexcept;

      /// Returns the function where it occurred
      const char* getFunction() const noexcept;

      /// Returns the message
      const char* getMessage() const noexcept;

      //@}

    protected:
      /// The source file the exception was thrown in
      const char* file_;

      /// The line number the exception was thrown in
      int line_;

      /// The source file the exception was thrown in
      const char* function_;

      /// The name of the exception.
      std::string name_;
    };

    /**
        @brief Precondition failed exception.

        A precondition (as defined by @ref OPENMS_PRECONDITION ) has failed.

        @ingroup Exceptions
    */
    class OPENMS_DLLAPI Precondition : public BaseException
    {
    public:
      Precondition(const char* file, int line, const char* function, const std::string& condition) noexcept;
    };

    /**
      @brief Postcondition failed exception.

      A postcondition (as defined by @ref OPENMS_POSTCONDITION ) has failed.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI Postcondition : public BaseException
    {
    public:
      Postcondition(const char* file, int line, const char* function, const std::string& condition) noexcept;
    };

    /**
      @brief Not all required information provided.

      Information that are required are not provided.
      Especially useful for missing MetaInfo values.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI MissingInformation : public BaseException
    {
    public:
      MissingInformation(const char* file, int line, const char* function, const std::string& error_message) noexcept;
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
    class OPENMS_DLLAPI IndexUnderflow : public BaseException
    {
    public:
      IndexUnderflow(const char* file, int line, const char* function, SignedSize index = 0, Size size = 0) noexcept;
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
    class OPENMS_DLLAPI SizeUnderflow : public BaseException
    {
    public:
      SizeUnderflow(const char* file, int line, const char* function, Size size = 0) noexcept;
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
    class OPENMS_DLLAPI IndexOverflow : public BaseException
    {
    public:
      IndexOverflow(const char* file, int line, const char* function, SignedSize index = 0, Size size = 0) noexcept;
    };


    /**
      @brief Array not sorted exception

      Throw this exception to indicate that an array/vector of elements
      was expected to be sorted, but was found to be unsorted.

      @param	message What was unsorted?

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI NotSorted : public BaseException
    {
    public:
      NotSorted(const char* file, int line, const char* function, const std::string& message) noexcept;
    };

    /**
      @brief A call to an external library (other than OpenMS) went wrong.

      Throw this exception to indicate that an external library call came
      back unsuccessful.

      @param	size the size causing the problem

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI FailedAPICall : public BaseException
    {
    public:
      FailedAPICall(const char* file, int line, const char* function, const std::string& message) noexcept;
    };

    /**
      @brief Invalid range exception.

      Use this exception to indicate a general range problems.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI InvalidRange : public BaseException
    {
    public:
      InvalidRange(const char* file, int line, const char* function) noexcept;
      InvalidRange(const char* file, int line, const char* function, const std::string& message) noexcept;
    };


    /**
      @brief Invalid UInt exception.

      Throw this exception to indicate that a size was unexpected.
      The constructor has an additional argument: the value of of the
      requested size.
      @param	size the size causing the problem

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI InvalidSize : public BaseException
    {
    public:
      InvalidSize(const char* file, int line, const char* function, Size size = 0) noexcept;
    };


    /**
      @brief Out of range exception.

      Use this exception to indicate that a given value is out of a
      defined range, i. e. not within the domain of a function.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI OutOfRange : public BaseException
    {
    public:
      OutOfRange(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Invalid value exception.

      Use this exception to indicate that a given value is not valid,
      when the value is only allowed to be out of a certain set of values.
      If the value has to be inside a given range, you should rather use OutOfRange.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI InvalidValue : public BaseException
    {
    public:
      InvalidValue(const char* file, int line, const char* function, const std::string& message, const std::string& value) noexcept;
    };

    /**
      @brief Exception indicating that an invalid parameter was handed over to an algorithm.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI InvalidParameter : public BaseException
    {
    public:
      InvalidParameter(const char* file, int line, const char* function, const std::string& message) noexcept;
    };

    /**
      @brief Invalid conversion exception.

      This exception indicates a conversion problem when converting from
      one type to another.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI ConversionError : public BaseException
    {
    public:
      ConversionError(const char* file, int line, const char* function, const std::string& error) noexcept;
    };

    /**
      @brief Illegal self operation exception.

      Throw this exception to indicate an invalid operation on the object
      itself. In general these operations are self assignments or related
      methods.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI IllegalSelfOperation : public BaseException
    {
    public:
      IllegalSelfOperation(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Null pointer argument is invalid exception.

      Use this exception to indicate a failure due to an argument not
      containing a pointer to a valid object, but a null pointer.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI NullPointer : public BaseException
    {
    public:
      NullPointer(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Invalid iterator exception.

      The iterator on which an operation should be performed was invalid.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI InvalidIterator : public BaseException
    {
    public:
      InvalidIterator(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Incompatible iterator exception.

      The iterators could not be assigned because they are bound to
      different containers.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI IncompatibleIterators : public BaseException
    {
    public:
      IncompatibleIterators(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Not implemented exception.

      This exception should be thrown to indicate not yet implemented methods.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI NotImplemented : public BaseException
    {
    public:
      NotImplemented(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Illegal tree operation exception.

      This exception is thrown to indicate that an illegal tree operation was requested.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI IllegalTreeOperation : public BaseException
    {
    public:
      IllegalTreeOperation(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Out of memory exception.

      Throw this exception to indicate that an allocation failed.
      This exception is thrown in the OPENMS new handler.
      @param	size	the number of bytes that should have been allocated
      @see GlobalException::newHandler

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI OutOfMemory : public BaseException, public std::bad_alloc
    {
    public:
      OutOfMemory(const char* file, int line, const char* function, Size size = 0) noexcept;
    };

    /**
      @brief Buffer overflow exception.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI BufferOverflow : public BaseException
    {
    public:
      BufferOverflow(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Division by zero error exception.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI DivisionByZero : public BaseException
    {
    public:
      DivisionByZero(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief Out of grid exception.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI OutOfGrid : public BaseException
    {
    public:
      OutOfGrid(const char* file, int line, const char* function) noexcept;
    };

    /**
      @brief File not found exception.

      A given file could not be found.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI FileNotFound : public BaseException
    {
    public:
      FileNotFound(const char* file, int line, const char* function, const std::string& filename) noexcept;
    };

    /**
      @brief External executable (e.g. comet.exe) not found exception.

      A given file could not be found. Usually used in Adapters for external tools.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI ExternalExecutableNotFound : public BaseException
    {
    public:
      ExternalExecutableNotFound(const char* file, int line, const char* function, const std::string& filename) noexcept;
    };

    /**
      @brief File not readable exception.

      A given file is not readable for the current user.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI FileNotReadable : public BaseException
    {
    public:
      FileNotReadable(const char* file, int line, const char* function, const std::string& filename) noexcept;
    };

    /**
      @brief File not writable exception.

      A given file is not writable for the current user.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI FileNotWritable : public BaseException
    {
    public:
      FileNotWritable(const char* file, int line, const char* function, const std::string& filename) noexcept;
    };

    /**
    @brief Filename is too long to be writable/readable by the filesystem

    This exception usually occurs when output filenames are automatically generated
    and are found to be too long. Usually 255 characters is the limit (NTFS, Ext2/3/4,...).

    @ingroup Exceptions
    */
    class OPENMS_DLLAPI FileNameTooLong : public BaseException
    {
    public:
      FileNameTooLong(const char* file, int line, const char* function, const std::string& filename, int max_length) noexcept;
    };

    /**
      @brief General IOException.

      General error for IO operations, that can not be associated to the more specific exceptions (e.g. FileNotWritable)

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI IOException : public BaseException
    {
    public:
      IOException(const char* file, int line, const char* function, const std::string& filename) noexcept;
    };

    /**
      @brief SqlOperation failed exception.

      E.g. when retrieving data from a table using the wrong column name or index.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI SqlOperationFailed : public BaseException
    {
    public:
      SqlOperationFailed(const char* file, int line, const char* function, const std::string& description) noexcept;
    };

    /**
      @brief File is empty.

      A given file is empty.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI FileEmpty : public BaseException
    {
    public:
      FileEmpty(const char* file, int line, const char* function, const std::string& filename) noexcept;
    };

    /**
      @brief Invalid 3-dimensional position exception.

      A given position in three dimensional is invalid.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI IllegalPosition : public BaseException
    {
    public:
      IllegalPosition(const char* file, int line, const char* function, float x, float y, float z) noexcept;
    };

    /**
      @brief Parse Error exception.

      A given expression could not be parsed.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI ParseError : public BaseException
    {
    public:
      ParseError(const char* file, int line, const char* function, const std::string& expression, const std::string& message) noexcept;
    };

    /**
      @brief Unable to create file exception.

      The given file could not be created.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI UnableToCreateFile : public BaseException
    {
    public:
      UnableToCreateFile(const char* file, int line, const char* function, const std::string& filename, const std::string& message = "") noexcept;
    };

    /**
      @brief Invalid file type exception.

      The file type specification is not valid.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI InvalidFileType : public BaseException
    {
    public:
      InvalidFileType(const char* file, int line, const char* function, const std::string& filename, const std::string& message = "") noexcept;
    };

    /**
      @brief A method or algorithm argument contains illegal values

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI IllegalArgument : public BaseException
    {
    public:
      IllegalArgument(const char* file, int line, const char* function, const std::string& error_message) noexcept;
    };

    /**
      @brief A tool or algorithm which was called internally raised an exception

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI InternalToolError : public BaseException
    {
    public:
      InternalToolError(const char* file, int line, const char* function, const std::string& error_message) noexcept;
    };

    /**
      @brief Element could not be found exception.

      The given element could not be found.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI ElementNotFound : public BaseException
    {
    public:
      ElementNotFound(const char* file, int line, const char* function, const std::string& element) noexcept;
    };

    /**
      @brief Exception used if an error occurred while fitting a model to a given dataset

      The given element could not be found.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI UnableToFit : public BaseException
    {
    public:
      UnableToFit(const char* file, int line, const char* function, const std::string& name, const std::string& message) noexcept;
    };

    /**
      @brief Exception used if an error occurred while calibrating a dataset.

      The calibration can not be performed because not enough reference masses
      were detected.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI UnableToCalibrate : public BaseException
    {
    public:
      UnableToCalibrate(const char* file, int line, const char* function, const std::string& name, const std::string& message) noexcept;
    };

    /**
      @brief Exception used if no more unique document ID's can be drawn from ID pool.

      The ID pool of OpenMS is either depleted or not existent.

      @ingroup Exceptions
    */
    class OPENMS_DLLAPI DepletedIDPool : public BaseException
    {
    public:
      DepletedIDPool(const char* file, int line, const char* function, const std::string& name, const std::string& message) noexcept;
    };

  } // namespace Exception

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
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Exception::BaseException& e);

} // namespace OpenMS
