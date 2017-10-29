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

#include <OpenMS/FORMAT/XMLFile.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

#include <OpenMS/FORMAT/CompressedInputSource.h>

#include <xercesc/framework/XMLFormatter.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <fstream>
#include <iomanip> // setprecision etc.
#include <iosfwd>

#include <boost/shared_ptr.hpp>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    /// This class ensures that the reset() method of the XMLHandler is called when it goes out of scope.
    /// useful when used in exception handling
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
      // ensure handler->reset() is called to save memory (in case the XMLFile
      // reader, e.g. FeatureXMLFile, is used again)
      XMLCleaner_ clean(handler);

      StringManager sm;
      //try to open file
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      // initialize parser
      try
      {
        xercesc::XMLPlatformUtils::Initialize();
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "", String("Error during initialization: ") + StringManager().convert(toCatch.getMessage()));
      }

      boost::shared_ptr< xercesc::SAX2XMLReader > parser(xercesc::XMLReaderFactory::createXMLReader());
      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces, false);
      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes, false);

      parser->setContentHandler(handler);
      parser->setErrorHandler(handler);

      // peak ahead into the file: is it bzip2 or gzip compressed?
      String bz;
      {
        std::ifstream file(filename.c_str());
        char tmp_bz[3];
        file.read(tmp_bz, 2);
        tmp_bz[2] = '\0';
        bz = String(tmp_bz);
      }

      boost::shared_ptr< xercesc::InputSource > source;

      char g1 = 0x1f;
      char g2 = 0;
      g2 |= 1 << 7;
      g2 |= 1 << 3;
      g2 |= 1 << 1;
      g2 |= 1 << 0;
      //g2 = static_cast<char>(0x8b); // can make troubles if it is casted to 0x7F which is the biggest number signed char can save
      if ((bz[0] == 'B' && bz[1] == 'Z') || (bz[0] == g1 && bz[1] == g2))
      {
        source.reset(new CompressedInputSource(sm.convert(filename).c_str(), bz));
      }
      else
      {
        source.reset(new xercesc::LocalFileInputSource(sm.convert(filename).c_str()));
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
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", 
            String("XMLException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const xercesc::SAXException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
            String("SAXException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const XMLHandler::EndParsingSoftly & /*toCatch*/)
      {
        // nothing to do here, as this exception is used to softly abort the
        // parsing for whatever reason.
      }
      catch (...)
      {
        // re-throw
        throw;
      }
    }

    void XMLFile::parseBuffer_(const std::string & buffer, XMLHandler * handler)
    {
      // ensure handler->reset() is called to save memory (in case the XMLFile
      // reader, e.g. FeatureXMLFile, is used again)
      XMLCleaner_ clean(handler);

      StringManager sm;

      // initialize parser
      try
      {
        xercesc::XMLPlatformUtils::Initialize();
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "", String("Error during initialization: ") + StringManager().convert(toCatch.getMessage()));
      }

      boost::shared_ptr< xercesc::SAX2XMLReader > parser(xercesc::XMLReaderFactory::createXMLReader());
      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces, false);
      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes, false);

      parser->setContentHandler(handler);
      parser->setErrorHandler(handler);

      // TODO: handle non-plaintext
      // peak ahead into the file: is it bzip2 or gzip compressed?
      // String bz = buffer.substr(0, 2);

      boost::shared_ptr< xercesc::InputSource > source;
      {
        auto fake_id = sm.convert("inMemory");
        source.reset(new xercesc::MemBufInputSource(reinterpret_cast<const unsigned char *>(buffer.c_str()), buffer.size(), fake_id.c_str()));
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
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", 
            String("XMLException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const xercesc::SAXException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
            String("SAXException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const XMLHandler::EndParsingSoftly & /*toCatch*/)
      {
        // nothing to do here, as this exception is used to softly abort the
        // parsing for whatever reason.
      }
      catch (...)
      {
        // re-throw
        throw;
      }
    }

    void XMLFile::save_(const String & filename, XMLHandler * handler) const
    {
      // open file in binary mode to avoid any line ending conversions
      std::ofstream os(filename.c_str(), std::ios::out | std::ios::binary);

      //set high precision for writing of floating point numbers
      os.precision(writtenDigits(double()));

      if (!os)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      // write data and close stream
      handler->writeTo(os);
      os.close();
    }

    String encodeTab(const String& to_encode)
    {
      if (!to_encode.has('\t')) return to_encode;
      else return String(to_encode).substitute("\t", "&#x9;");
    }

    bool XMLFile::isValid(const String & filename, std::ostream & os)
    {
      if (schema_location_.empty())
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
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
