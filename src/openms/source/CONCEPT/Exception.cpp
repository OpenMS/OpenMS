// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>

#include <cstdio>
#include <iostream>
#include <sstream>
#include <typeinfo>

#define DEF_EXCEPTION(a, b) \
  a::a(const char* file, int line, const char* function) noexcept :                                                           \
    BaseException(file, line, function, # a, b) \
  { \
  } \

using namespace std;

namespace OpenMS
{

  namespace Exception
  {

    BaseException::BaseException() noexcept :
      std::runtime_error("unknown error"),
      file_("?"),
      line_(-1),
      function_("?"),
      name_("Exception")
    {
      GlobalExceptionHandler::getInstance().set(file_, line_, function_, name_, what());
    }

    BaseException::BaseException(const char* file, int line, const char* function, const std::string& name, const std::string& message) noexcept :
      std::runtime_error(message),
      file_(file),
      line_(line),
      function_(function),
      name_(name)
    {
      GlobalExceptionHandler::getInstance().set(file_, line_, function_, name_, what());
    }

    BaseException::BaseException(const char* file, int line, const char* function) noexcept :
      std::runtime_error("unknown error"),
      file_(file),
      line_(line),
      function_(function),
      name_("Exception")
    {
      GlobalExceptionHandler::getInstance().set(file_, line_, function_, name_, what());
    }

    BaseException::BaseException(const BaseException& exception) noexcept :
      std::runtime_error(exception),
      file_(exception.file_),
      line_(exception.line_),
      function_(exception.function_),
      name_(exception.name_)
    {
    }

    BaseException::~BaseException() noexcept
    {
    }

    const char* BaseException::getName() const noexcept
    {
      return name_.c_str();
    }

    const char* BaseException::getFile() const noexcept
    {
      return file_;
    }

    const char* BaseException::getFunction() const noexcept
    {
      return function_;
    }

    const char* BaseException::getMessage() const noexcept
    {
      return what();
    }

    int BaseException::getLine() const noexcept
    {
      return line_;
    }

    Precondition::Precondition(const char* file, int line, const char* function, const string& condition) noexcept :
      BaseException(file, line, function, "Precondition failed", std::string(condition))
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    Postcondition::Postcondition(const char* file, int line, const char* function, const string& condition) noexcept :
      BaseException(file, line, function, "Postcondition failed", std::string(condition))
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    IndexUnderflow::IndexUnderflow(const char* file, int line, const char* function, SignedSize index, Size size) noexcept :
      BaseException(file, line, function, "IndexUnderflow", "the given index was too small: " + String(index) + " (size = " + String(size) + ")")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    IndexOverflow::IndexOverflow(const char* file, int line, const char* function, SignedSize index, Size size) noexcept :
      BaseException(file, line, function, "IndexOverflow", "the given index was too large: " + String(index) + " (size = " + String(size) + ")")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    NotSorted::NotSorted(const char* file, int line, const char* function, const std::string& message) noexcept:
      BaseException(file, line, function, "NotSorted", message)
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    FailedAPICall::FailedAPICall(const char* file, int line, const char* function, const std::string& message) noexcept :
      BaseException(file, line, function, "FailedAPICall", message)
    {
    }

    OutOfMemory::OutOfMemory(const char* file, int line, const char* function, Size size) noexcept :
      BaseException(file, line, function, "OutOfMemory", "unable to allocate enough memory (size = " + String(size) + " bytes) ")
    {
      GlobalExceptionHandler::getInstance().setMessage(static_cast<std::runtime_error>(*this).what());
    }

    SizeUnderflow::SizeUnderflow(const char* file, int line, const char* function, Size size) noexcept :
      BaseException(file, line, function, "SizeUnderflow", "the given size was too small: " + String(size))
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    InvalidSize::InvalidSize(const char* file, int line, const char* function, Size size) noexcept :
      BaseException(file, line, function, "InvalidSize", "the given size was not expected: " + String(size))
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    IllegalPosition::IllegalPosition(const char* file, int line, const char* function, float x, float y, float z) noexcept :
      BaseException(file, line, function, "IllegalPosition:", "(" + String(x) + "," + String(y) + "," + String(z) + ")")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    ParseError::ParseError(const char* file, int line, const char* function, const std::string& expression, const std::string& message) noexcept :
      BaseException(file, line, function, "Parse Error", message + " in: " + expression)
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    FileNotFound::FileNotFound(const char* file, int line, const char* function, const std::string& filename) noexcept :
      BaseException(file, line, function, "FileNotFound", "the file '" + filename + "' could not be found")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    ExternalExecutableNotFound::ExternalExecutableNotFound(const char* file, int line, const char* function, const std::string& filename) noexcept:
        BaseException(file, line, function, "ExternalExecutableNotFound", "the executable '" + filename + "' could not be found")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }
    

    FileNotReadable::FileNotReadable(const char* file, int line, const char* function, const std::string& filename) noexcept :
      BaseException(file, line, function, "FileNotReadable", "the file '" + filename + "' is not readable for the current user")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    FileNotWritable::FileNotWritable(const char* file, int line, const char* function, const std::string& filename) noexcept :
      BaseException(file, line, function, "FileNotWritable", "the file '" + filename + "' is not writable for the current user")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    FileNameTooLong::FileNameTooLong(const char* file, int line, const char* function, const std::string& filename, int max_length) noexcept :
      BaseException(file, line, function, "FileNameTooLong", 
        "the file '" + filename + "' is too long (" + String(filename.size()) + " chars) "
         + "and exceeds the allowed limit of " + String(max_length) + "; "
         + "use shorter filenames and/or fewer subdirectories.")      
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    IOException::IOException(const char* file, int line, const char* function, const std::string& filename) noexcept :
      BaseException(file, line, function, "IOException", "IO error for file '" + filename + "'")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    SqlOperationFailed::SqlOperationFailed(const char* file, int line, const char* function, const std::string& description) noexcept :
      BaseException(file, line, function, "SqlOperationFailed", "an sql operation failed ('" + description + "')")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    FileEmpty::FileEmpty(const char* file, int line, const char* function, const std::string& filename) noexcept :
      BaseException(file, line, function, "FileEmpty", "the file '" + filename + "' is empty")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    ConversionError::ConversionError(const char* file, int line, const char* function, const std::string& error) noexcept :
      BaseException(file, line, function, "ConversionError", error)
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    InvalidValue::InvalidValue(const char* file, int line, const char* function, const std::string& message, const std::string& value) noexcept :
      BaseException(file, line, function, "InvalidValue", "the value '" + value + "' was used but is not valid; " + message)
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    InvalidParameter::InvalidParameter(const char* file, int line, const char* function, const std::string& message) noexcept :
      BaseException(file, line, function, "InvalidParameter", message)
    {
    }

    UnableToCreateFile::UnableToCreateFile(const char* file, int line, const char* function, const std::string& filename, const std::string& message) noexcept :
      BaseException(file, line, function, "UnableToCreateFile", "the file '" + filename + "' could not be created. " + message)
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    InvalidFileType::InvalidFileType(const char* file, int line, const char* function, const std::string& filename, const std::string& message) noexcept :
      BaseException(file, line, function, "InvalidFileType", "the file '" + filename + "' could not be created because the type specified was not valid. " + message)
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    IllegalArgument::IllegalArgument(const char* file, int line, const char* function, const string& error_message) noexcept :
      BaseException(file, line, function, "IllegalArgument", error_message)
    {
    }

    InternalToolError::InternalToolError(const char* file, int line, const char* function, const std::string& error_message) noexcept:
      BaseException(file, line, function, "InternalToolError", error_message)
    {
    }

    MissingInformation::MissingInformation(const char* file, int line, const char* function, const string& error_message) noexcept :
      BaseException(file, line, function, "MissingInformation", error_message)
    {
    }

    ElementNotFound::ElementNotFound(const char* file, int line, const char* function, const string& element)   noexcept :
      BaseException(file, line, function, "ElementNotFound", "the element '" + element + "' could not be found")
    {
      GlobalExceptionHandler::getInstance().setMessage(what());
    }

    UnableToFit::UnableToFit(const char* file, int line, const char* function, const string& name, const string& message) noexcept :
      BaseException(file, line, function, name, message)
    {
    }

    UnableToCalibrate::UnableToCalibrate(const char* file, int line, const char* function, const string& name, const string& message) noexcept :
      BaseException(file, line, function, name, message)
    {
    }

    DepletedIDPool::DepletedIDPool(const char* file, int line, const char* function, const string& name, const string& message) noexcept :
      BaseException(file, line, function, name, message)
    {
    }

    InvalidRange::InvalidRange(const char* file, int line, const char* function) noexcept :
      BaseException(file, line, function, "InvalidRange", "the range of the operation was invalid")
    {
    }

    InvalidRange::InvalidRange(const char* file, int line, const char* function, const std::string& message) noexcept :
      BaseException(file, line, function, "InvalidRange", message)
    {
    }

    DEF_EXCEPTION(DivisionByZero, "a division by zero was requested")

    DEF_EXCEPTION(OutOfRange, "the argument was not in range")

    DEF_EXCEPTION(NullPointer, "a null pointer was specified")

    DEF_EXCEPTION(InvalidIterator, "the iterator is invalid - probably it is not bound to a container")

    DEF_EXCEPTION(IncompatibleIterators, "the iterator could not be assigned because it is bound to a different container")

    DEF_EXCEPTION(NotImplemented, "this method has not been implemented yet. Feel free to complain about it!")

    DEF_EXCEPTION(IllegalSelfOperation, "cannot perform operation on the same object")

    DEF_EXCEPTION(IllegalTreeOperation, "an illegal tree operation was requested")

    DEF_EXCEPTION(BufferOverflow, "the maximum buffersize has been reached")

    DEF_EXCEPTION(OutOfGrid, "a point was outside a grid")

  } // namespace Exception

  std::ostream& operator<<(std::ostream& os, const Exception::BaseException& e)
  {
    os << e.getName() << " @ " << e.getFile() << ":" << e.getFunction() << " (Line " << e.getLine() << "): " << e.what();

    return os;
  }

} // namespace OPENMS
