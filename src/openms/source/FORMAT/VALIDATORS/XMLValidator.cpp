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
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/util/XMLString.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
  // Internal bridge: adapts the Xerces ErrorHandler interface to XMLValidator so
  // that Xerces types stay out of the public XMLValidator header.
  class XMLValidatorErrorHandler_ : public xercesc::ErrorHandler
  {
  public:
    explicit XMLValidatorErrorHandler_(XMLValidator& v) : v_(v) {}
    void warning(const SAXParseException& e) override { report_("warning", e); }
    void error(const SAXParseException& e) override { report_("error", e); }
    void fatalError(const SAXParseException& e) override { report_("error", e); }
    void resetErrors() override { v_.valid_ = true; }
  private:
    void report_(const std::string& kind, const SAXParseException& e)
    {
      v_.logError_(kind, Internal::StringManager::convert(reinterpret_cast<const char16_t*>(e.getMessage())),
                   (unsigned long)e.getLineNumber(), (unsigned long)e.getColumnNumber());
    }
    XMLValidator& v_;
  };

  XMLValidator::XMLValidator() :
    valid_(true),
    os_(nullptr)
  {
  }

  bool XMLValidator::isValid(const std::string & filename, const std::string & schema, std::ostream & os)
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
        std::string("Error during initialization: ") + Internal::StringManager().convert(toCatch.getMessage()));
    }

    SAX2XMLReader * parser = XMLReaderFactory::createXMLReader();
    parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, true);
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgXercesDynamic, false);
    parser->setFeature(XMLUni::fgXercesSchema, true);
    parser->setFeature(XMLUni::fgXercesSchemaFullChecking, true);

    //set the internal Xerces error-handler bridge
    XMLValidatorErrorHandler_ error_handler(*this);
    parser->setErrorHandler(&error_handler);
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

  void XMLValidator::logError_(const std::string& kind, const std::string& message, unsigned long line, unsigned long column)
  {
    std::string error_message = std::string("Validation ") + kind + " in file '" + filename_ + "' line " + (UInt)line + " column " + (UInt)column + ": " + message;
    (*os_) << error_message << endl;
    valid_ = false;
  }

} // namespace OpenMS
