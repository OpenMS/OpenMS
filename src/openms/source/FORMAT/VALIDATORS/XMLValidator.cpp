// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
  XMLValidator::XMLValidator() :
    valid_(true),
    os_(nullptr)
  {
  }

  bool XMLValidator::isValid(const String & filename, const String & schema, std::ostream & os)
  {
    filename_ = filename;
    os_ = &os;

    //try to open file
    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    // initialize parser
    try
    {
      XMLPlatformUtils::Initialize();
    }
    catch (const XMLException & toCatch)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "", 
        String("Error during initialization: ") + Internal::StringManager().convert(toCatch.getMessage()));
    }

    SAX2XMLReader * parser = XMLReaderFactory::createXMLReader();
    parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, true);
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgXercesDynamic, false);
    parser->setFeature(XMLUni::fgXercesSchema, true);
    parser->setFeature(XMLUni::fgXercesSchemaFullChecking, true);

    //set this class as error handler
    parser->setErrorHandler(this);
    parser->setContentHandler(nullptr);
    parser->setEntityResolver(nullptr);

    //load schema
    LocalFileInputSource schema_file(Internal::StringManager().convert(schema).c_str());
    parser->loadGrammar(schema_file, Grammar::SchemaGrammarType, true);
    parser->setFeature(XMLUni::fgXercesUseCachedGrammarInParse, true);

    // try to parse file
    LocalFileInputSource source(Internal::StringManager().convert(filename.c_str()).c_str());

    try
    {
      parser->parse(source);
      delete(parser);
    }
    catch (...)
    {
      /// nothing to do here
    }

    return valid_;
  }

  void XMLValidator::warning(const SAXParseException & exception)
  {
    char * message = XMLString::transcode(exception.getMessage());
    String error_message = String("Validation warning in file '") + filename_ + "' line " + (UInt) exception.getLineNumber() + " column " + (UInt) exception.getColumnNumber() + ": " + message;
    (*os_) << error_message << endl;
    valid_ = false;
    XMLString::release(&message);
  }

  void XMLValidator::error(const SAXParseException & exception)
  {
    char * message = XMLString::transcode(exception.getMessage());
    String error_message = String("Validation error in file '") + filename_ + "' line " + (UInt) exception.getLineNumber() + " column " + (UInt) exception.getColumnNumber() + ": " + message;
    (*os_) << error_message << endl;
    valid_ = false;
    XMLString::release(&message);
  }

  void XMLValidator::fatalError(const SAXParseException & exception)
  {
    char * message = XMLString::transcode(exception.getMessage());
    String error_message = String("Validation error in file '") + filename_ + "' line " + (UInt) exception.getLineNumber() + " column " + (UInt) exception.getColumnNumber() + ": " + message;
    (*os_) << error_message << endl;
    valid_ = false;
    XMLString::release(&message);
  }

  void XMLValidator::resetErrors()
  {
    valid_ = true;
  }

} // namespace OpenMS
