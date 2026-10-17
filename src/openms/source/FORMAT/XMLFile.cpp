// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XMLFile.h>

#include <OpenMS/CONCEPT/Init.h>

#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/SAX2HandlerAdapter.h>
#ifdef WITH_S3
#include "DATAACCESS/S3InputSource.h"
#endif
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/FileNameUtils.h>

#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/framework/URLInputSource.hpp>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <zlib.h>
#include <bzlib.h>

#include <fstream>
#include <iomanip> // setprecision etc.
#include <streambuf>

#include <memory>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLURL.hpp>

using namespace std;

namespace OpenMS::Internal
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
    = default;

    XMLFile::XMLFile(const std::string & schema_location, const std::string & version) :
      schema_location_(schema_location),
      schema_version_(version)
    {
    }

    XMLFile::~XMLFile()
    = default;

    void XMLFile::enforceEncoding_(const std::string& encoding)
    {
      enforced_encoding_ = encoding;
    }

    void parse(xercesc::InputSource* const source, XMLHandler* handler)
    {
      unique_ptr<xercesc::SAX2XMLReader> parser(xercesc::XMLReaderFactory::createXMLReader());

      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces, false);
      parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes, false);

      // XMLHandler exposes only Xerces-free callbacks; bridge it to the SAX parser.
      SAX2HandlerAdapter adapter(*handler);
      parser->setContentHandler(&adapter);
      parser->setErrorHandler(&adapter);

      const bool remote = source->getSystemId()
                       && FileNameUtils::isRemoteResource(StringManager().convert(source->getSystemId()));

      // try to parse file
      try
      {
        parser->parse(*source);
      }
      catch (const xercesc::XMLException& toCatch)
      {
        if (remote)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                      "XMLException while reading remote XML; check network transport and XML encoding");
        }
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",std::string("XMLException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const xercesc::SAXException& toCatch)
      {
        if (remote)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "SAXException while parsing remote XML");
        }
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",std::string("SAXException: ") + StringManager().convert(toCatch.getMessage()));
      }
      catch (const XMLHandler::EndParsingSoftly& /*toCatch*/)
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

    void XMLFile::parse_(const std::string & filename, XMLHandler * handler)
    {
      // ensure handler->reset() is called to save memory (in case the XMLFile
      // reader, e.g. FeatureXMLFile, is used again)
      XMLCleaner_ clean(handler);

      StringManager sm;
      //try to open file
      if (!FileNameUtils::isRemoteResource(filename) && !File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      // initialize parser
      try
      {
        Internal::xercesInitialize(); // serialised: Xerces requires it, and parseBuffer_() runs on several threads at once
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "",std::string("Error during initialization: ") + StringManager().convert(toCatch.getMessage()));
      }


      unique_ptr<xercesc::InputSource> source;
      if (filename.compare(0, 5, "s3://") == 0)
      {
#ifdef WITH_S3
        source = std::make_unique<S3InputSource>(filename);
#else
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                    "S3 XML input requires an OpenMS build with WITH_S3 enabled");
#endif
      }
      else if (FileNameUtils::isRemoteResource(filename))
      {
        // Availability of HTTPS/FTP depends on Xerces' configured network accessor.
        try
        {
          source = std::make_unique<xercesc::URLInputSource>(xercesc::XMLURL(sm.convert(filename).c_str()));
          source->setSystemId(sm.convert(FileNameUtils::redactRemoteResource(filename)).c_str());
        }
        catch (const xercesc::XMLException&)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                      "Invalid remote XML URL");
        }
      }
      else
      {
        char header[2] = {};
        std::ifstream file(filename, std::ios::binary);
        file.read(header, sizeof(header));
        const std::string magic(header, static_cast<size_t>(file.gcount()));
        if (magic == "BZ" || magic == std::string("\x1f\x8b", 2) || magic == "PK")
        {
          source = std::make_unique<CompressedInputSource>(sm.convert(filename).c_str(), magic);
        }
        else
        {
          source = std::make_unique<xercesc::LocalFileInputSource>(sm.convert(filename).c_str());
        }
      }
      // what if no encoding given http://xerces.apache.org/xerces-c/apiDocs-3/classInputSource.html
      if (!enforced_encoding_.empty())
      {
        static const XMLCh* s_enc = xercesc::XMLString::transcode(enforced_encoding_.c_str());
        source->setEncoding(s_enc);
      }
      
      parse(source.get(), handler);
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
        Internal::xercesInitialize(); // serialised: Xerces requires it, and parseBuffer_() runs on several threads at once
      }
      catch (const xercesc::XMLException & toCatch)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "",std::string("Error during initialization: ") + StringManager().convert(toCatch.getMessage()));
      }

      // TODO: handle non-plain text
      // peak ahead into the file: is it bzip2 or gzip compressed?
      // std::string bz = StringUtils::substr(buffer, 0, 2);

      unique_ptr<xercesc::InputSource> source;
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
      
      parse(source.get(), handler);
    }

    void XMLFile::save_(const std::string & filename, XMLHandler * handler) const
    {
      // Detect compression from filename extension, in any letter case like on input (and in TOPPBase's check)
      const FileTypes::Type compression = FileNameUtils::compressionType(filename);
      const bool use_gzip = compression == FileTypes::GZ;
      const bool use_bzip2 = compression == FileTypes::BZ2;

      if (use_gzip)
      {
        // Minimal std::streambuf that writes gzip-compressed data via zlib
        class GzipOStreambuf : public std::streambuf
        {
        public:
          enum { BUF_SIZE = 65536 };

          explicit GzipOStreambuf(const char* fn) : gz_(gzopen(fn, "wb"))
          {
            if (gz_) setp(buf_, buf_ + BUF_SIZE);
          }

          ~GzipOStreambuf()
          {
            flush_buffer();
            if (gz_) gzclose(gz_);
          }

          bool is_open() const { return gz_ != nullptr; }

        protected:
          int_type overflow(int_type c) override
          {
            if (!gz_) return traits_type::eof();
            if (pptr() == epptr() && !flush_buffer()) return traits_type::eof();
            if (!traits_type::eq_int_type(c, traits_type::eof()))
            {
              *pptr() = traits_type::to_char_type(c);
              pbump(1);
            }
            return traits_type::not_eof(c);
          }

          int sync() override { return flush_buffer() ? 0 : -1; }

        private:
          bool flush_buffer()
          {
            int n = static_cast<int>(pptr() - pbase());
            if (n > 0)
            {
              if (gzwrite(gz_, pbase(), static_cast<unsigned>(n)) != n) return false;
              setp(buf_, buf_ + BUF_SIZE);
            }
            return true;
          }

          gzFile gz_;
          char buf_[BUF_SIZE];
        };

        GzipOStreambuf sbuf(filename.c_str());
        if (!sbuf.is_open())
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
        }
        std::ostream os(&sbuf);
        os.precision(writtenDigits(double()));
        handler->writeTo(os);
      }
      else if (use_bzip2)
      {
        // Minimal std::streambuf that writes bzip2-compressed data via bzlib
        class Bzip2OStreambuf : public std::streambuf
        {
        public:
          enum { BUF_SIZE = 65536 };

          explicit Bzip2OStreambuf(const char* fn) : file_(nullptr), bz_(nullptr)
          {
            file_ = fopen(fn, "wb");
            if (!file_) return;
            int err = BZ_OK;
            bz_ = BZ2_bzWriteOpen(&err, file_, 9, 0, 0);
            if (err != BZ_OK)
            {
              BZ2_bzWriteClose(&err, bz_, 1, nullptr, nullptr);
              bz_ = nullptr;
              fclose(file_);
              file_ = nullptr;
              return;
            }
            setp(buf_, buf_ + BUF_SIZE);
          }

          ~Bzip2OStreambuf()
          {
            flush_buffer();
            if (bz_)
            {
              int err = BZ_OK;
              BZ2_bzWriteClose(&err, bz_, 0, nullptr, nullptr);
            }
            if (file_) fclose(file_);
          }

          bool is_open() const { return bz_ != nullptr; }

        protected:
          int_type overflow(int_type c) override
          {
            if (!bz_) return traits_type::eof();
            if (pptr() == epptr() && !flush_buffer()) return traits_type::eof();
            if (!traits_type::eq_int_type(c, traits_type::eof()))
            {
              *pptr() = traits_type::to_char_type(c);
              pbump(1);
            }
            return traits_type::not_eof(c);
          }

          int sync() override { return flush_buffer() ? 0 : -1; }

        private:
          bool flush_buffer()
          {
            if (!bz_) return false;
            int n = static_cast<int>(pptr() - pbase());
            if (n > 0)
            {
              int err = BZ_OK;
              BZ2_bzWrite(&err, bz_, pbase(), n);
              if (err != BZ_OK && err != BZ_RUN_OK) return false;
              setp(buf_, buf_ + BUF_SIZE);
            }
            return true;
          }

          FILE* file_;
          BZFILE* bz_;
          char buf_[BUF_SIZE];
        };

        Bzip2OStreambuf sbuf(filename.c_str());
        if (!sbuf.is_open())
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
        }
        std::ostream os(&sbuf);
        os.precision(writtenDigits(double()));
        handler->writeTo(os);
      }
      else
      {
        // Uncompressed: open in binary mode to avoid any line ending conversions
        std::ofstream os(filename.c_str(), std::ios::out | std::ios::binary);
        os.precision(writtenDigits(double()));
        if (!os)
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
        }
        handler->writeTo(os);
        os.close();
      }
    }

    std::string encodeTab(const std::string& to_encode)
    {
      if (!StringUtils::has(to_encode, '\t'))
      {
        return to_encode;
      }
      else
      {
        std::string result = to_encode;
        StringUtils::substitute(result, "\t", "&#x9;");
        return result;
      }
    }

    bool XMLFile::isValid(const std::string & filename, std::ostream & os)
    {
      if (schema_location_.empty())
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
      std::string current_location = File::find(schema_location_);
      return XMLValidator().isValid(filename, current_location, os);
    }

    const std::string & XMLFile::getVersion() const
    {
      return schema_version_;
    }

} // namespace OpenMS  // namespace Internal
