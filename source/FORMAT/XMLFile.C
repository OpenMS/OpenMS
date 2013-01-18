// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

#include <OpenMS/FORMAT/CompressedInputSource.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <fstream>
#include <iomanip> // setprecision etc.

using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    /// This class ensures that the reset() method of the XMLHandler is called when it goes out of scope.
    /// useful when used in exeption handling
    class XMLCleaner_
    {
public:
      explicit XMLCleaner_(XMLHandler * handler) :
        p_(handler)
      {

      }

      ~XMLCleaner_()
      {
        p_->reset();
      }

private:
      XMLHandler * p_;
    };

    XMLFile::XMLFile()
    {
    }

    XMLFile::XMLFile(const String & schema_location, const String & version) :
      schema_location_(schema_location),
      schema_version_(version)
    {
    }

    XMLFile::~XMLFile()
    {
    }

    void XMLFile::enforceEncoding_(const String& encoding)
    {
      enforced_encoding_ = encoding;
    }

    void XMLFile::parse_(const String & filename, XMLHandler * handler)
    {
      // ensure handler->reset() is called to save memory (in case the XMLFile reader, e.g. FatureXMLFile, is used again)
      XMLCleaner_ clean(handler);

      //try to open file
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      }

      // initialize parser
      try
      {
        xercesc::XMLPlatformUtils::Initialize();
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + StringManager().convert(toCatch.getMessage()));
      }

      xercesc::SAX2XMLReader * parser = xercesc::XMLReaderFactory::createXMLReader();
      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces, false);
      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes, false);

      parser->setContentHandler(handler);
      parser->setErrorHandler(handler);

      //is it bzip2 or gzip compressed?
      std::ifstream file(filename.c_str());
      char bz[2];
      file.read(bz, 2);
      xercesc::InputSource * source;

      char g1 = 0x1f;
      char g2 = 0;
      g2 |= 1 << 7;
      g2 |= 1 << 3;
      g2 |= 1 << 1;
      g2 |= 1 << 0;
      //g2 = static_cast<char>(0x8b); // can make troubles if it is casted to 0x7F which is the biggest number signed char can save
      if ((bz[0] == 'B' && bz[1] == 'Z') ||    (bz[0] == g1 && bz[1] == g2))
      {
        source = new CompressedInputSource(StringManager().convert(filename.c_str()), bz);
      }
      else
      {
        source = new xercesc::LocalFileInputSource(StringManager().convert(filename.c_str()));
      }
      // what if no encoding given http://xerces.apache.org/xerces-c/apiDocs-3/classInputSource.html
      if (!enforced_encoding_.empty())
      {
        static const XMLCh* s_enc = xercesc::XMLString::transcode(enforced_encoding_.c_str());
        source->setEncoding(s_enc);
      }
      // try to parse file
      try
      {
        parser->parse(*source);
        delete(parser);
        delete source;
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const xercesc::SAXException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const XMLHandler::EndParsingSoftly & /*toCatch*/)
      {
        //nothing to do here, as this exception is used to softly abort the parsing for whatever reason.
      }
    }

    void XMLFile::save_(const String & filename, XMLHandler * handler) const
    {
      std::ofstream os(filename.c_str());

      //set high precision for writing of floating point numbers
      os.precision(writtenDigits(DoubleReal()));

      if (!os)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      }

      // write data and close stream
      handler->writeTo(os);
      os.close();
    }

    void writeXMLEscape(const String & to_escape, ostream & os)
    {
      XMLCh * xmlch = xercesc::XMLString::transcode(to_escape.c_str());

      std::string out = "";
      OpenMSXMLFormatTarget ft(out);
      xercesc::XMLFormatter f("UTF-8", /* XMLUni::fgVersion1_1 */ "1.1", &ft);
      f << xercesc::XMLFormatter::StdEscapes << xmlch;
      os << out;
      xercesc::XMLString::release(&xmlch);
      return;

    }

    String writeXMLEscape(const String & to_escape)
    {
      stringstream ss;
      writeXMLEscape(to_escape, ss);
      return String(ss.str());
    }

    bool XMLFile::isValid(const String & filename, std::ostream & os)
    {
      if (schema_location_.empty())
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      String current_location = File::find(schema_location_);
      return XMLValidator().isValid(filename, current_location, os);
    }

    const String & XMLFile::getVersion() const
    {
      return schema_version_;
    }

  }   // namespace Internal
} // namespace OpenMS
