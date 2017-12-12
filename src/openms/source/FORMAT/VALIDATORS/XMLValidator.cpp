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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/validators/common/Grammar.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
  XMLValidator::XMLValidator() :
    valid_(true),
    os_(0)
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
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Error during initialization: ") + Internal::StringManager().convert(toCatch.getMessage()));
    }

    SAX2XMLReader * parser = XMLReaderFactory::createXMLReader();
    parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, true);
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgXercesDynamic, false);
    parser->setFeature(XMLUni::fgXercesSchema, true);
    parser->setFeature(XMLUni::fgXercesSchemaFullChecking, true);

    //set this class as error handler
    parser->setErrorHandler(this);
    parser->setContentHandler(NULL);
    parser->setEntityResolver(NULL);

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
