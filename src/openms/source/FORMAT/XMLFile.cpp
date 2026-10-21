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
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/FileNameUtils.h>

#include <libxml/parser.h>
#include <libxml/parserInternals.h>
#include <exception>
#include <algorithm>
#include <cstring>


#include <zlib.h>
#include <bzlib.h>

#include <fstream>
#include <filesystem>
#include <iomanip> // setprecision etc.
#include <streambuf>

#include <memory>

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

    // C callbacks capture every exception: no C++ exception crosses libxml2 frames.
    struct LibxmlState
    {
      XMLHandler* handler;
      std::exception_ptr exception;
      static LibxmlState& state(void* ctx)
      {
        return *static_cast<LibxmlState*>(static_cast<xmlParserCtxt*>(ctx)->_private);
      }
      static std::u16string utf16(const xmlChar* text, int length = -1)
      {
        std::u16string result;
        if (!text) return result;
        const unsigned char* cursor = text;
        const unsigned char* end = text + (length < 0 ? xmlStrlen(text) : length);
        while (cursor < end)
        {
          unsigned int cp = *cursor++;
          if (cp >= 0xc0)
          {
            const int remaining = cp < 0xe0 ? 1 : cp < 0xf0 ? 2 : 3;
            cp &= (1u << (6 - remaining)) - 1;
            for (int i = 0; i < remaining && cursor < end; ++i) { cp = (cp << 6) | (*cursor++ & 0x3f); }
          }
          if (cp > 0xffff)
          {
            cp -= 0x10000;
            result.push_back(static_cast<char16_t>(0xd800 + (cp >> 10)));
            result.push_back(static_cast<char16_t>(0xdc00 + (cp & 0x3ff)));
          }
          else { result.push_back(static_cast<char16_t>(cp)); }
        }
        return result;
      }
      template <typename Callback> static void callback(void* ctx, Callback&& fn)
      {
        // Stopping the parser can re-enter an error callback. Keep the original
        // handler exception rather than replacing it with a parser-stop error.
        if (state(ctx).exception) { return; }
        try { fn(state(ctx)); }
        catch (...)
        {
          state(ctx).exception = std::current_exception();
          xmlStopParser(static_cast<xmlParserCtxt*>(ctx));
        }
      }
      static void startLegacy(void* ctx, const xmlChar* name, const xmlChar** attrs)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          XMLAttributes::Values values;
          if (attrs)
          {
            for (Size i = 0; attrs[i]; i += 2) values.emplace_back(utf16(attrs[i]), utf16(attrs[i + 1]));
          }
          s.handler->onStartElement(utf16(name).c_str(), XMLAttributes(values));
        });
      }
      static void endLegacy(void* ctx, const xmlChar* name)
      {
        callback(ctx, [&](LibxmlState& s) { s.handler->onEndElement(utf16(name).c_str()); });
      }
      static std::u16string qualifiedName(const xmlChar* local, const xmlChar* prefix)
      {
        if (!prefix) return utf16(local);
        auto name = utf16(prefix);
        name += u':';
        name += utf16(local);
        return name;
      }
      static void start(void* ctx, const xmlChar* local, const xmlChar* prefix,
                        const xmlChar*, int namespace_count, const xmlChar** namespaces,
                        int attribute_count, int, const xmlChar** attrs)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          XMLAttributes::Values values;
          values.reserve(namespace_count + attribute_count);
          // Existing handlers consume qualified names and can inspect xmlns declarations.
          for (int i = 0; i < namespace_count; ++i)
          {
            std::u16string name = u"xmlns";
            if (namespaces[2 * i]) { name += u':'; name += utf16(namespaces[2 * i]); }
            values.emplace_back(std::move(name), utf16(namespaces[2 * i + 1]));
          }
          for (int i = 0; i < attribute_count; ++i)
          {
            const xmlChar** attr = attrs + 5 * i;
            // SAX2 values are bounded slices, not null-terminated strings.
            values.emplace_back(qualifiedName(attr[0], attr[1]),
                                utf16(attr[3], static_cast<int>(attr[4] - attr[3])));
          }
          const auto name = qualifiedName(local, prefix);
          s.handler->onStartElement(name.c_str(), XMLAttributes(values));
        });
      }
      static void end(void* ctx, const xmlChar* local, const xmlChar* prefix, const xmlChar*)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          const auto name = qualifiedName(local, prefix);
          s.handler->onEndElement(name.c_str());
        });
      }
      static void characters(void* ctx, const xmlChar* text, int length)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          if (s.handler->tryOnCharactersUTF8(reinterpret_cast<const char*>(text), static_cast<Size>(length))) { return; }
          const auto value = utf16(text, length);
          s.handler->onCharacters(value.c_str(), value.size());
        });
      }
      static void error(void* ctx, const xmlError* error)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          const std::string message = error->message ? error->message : "libxml2 parse error";
          if (error->level == XML_ERR_WARNING) { s.handler->warning(XMLHandler::LOAD, message, error->line, error->int2); }
          else { s.handler->fatalError(XMLHandler::LOAD, message, error->line, error->int2); }
        });
      }
    };
    template <typename Read>
    void parse(Read&& read, XMLHandler* handler, const std::string& encoding, const std::string& filename, bool namespace_aware)
    {
      libxmlInitialize();
      LibxmlState state{handler, {}};
      xmlSAXHandler sax{};
      if (xmlSAXVersion(&sax, namespace_aware ? 2 : 1) != 0)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                    "libxml2 could not initialize the requested SAX callbacks");
      }
      if (namespace_aware)
      {
        sax.startElementNs = LibxmlState::start;
        sax.endElementNs = LibxmlState::end;
      }
      else
      {
        // Preserve generic OpenMS XML readers' historical namespace-disabled mode.
        sax.startElement = LibxmlState::startLegacy;
        sax.endElement = LibxmlState::endLegacy;
      }
      // OpenMS has no comment/PI callbacks. The libxml2 defaults build DOM nodes
      // even in SAX mode, retaining every comment and PI until the parse ends.
      sax.comment = nullptr;
      sax.processingInstruction = nullptr;
      sax.characters = LibxmlState::characters;
      sax.cdataBlock = LibxmlState::characters;
      auto free_parser = [](xmlParserCtxt* ctx)
      {
        if (ctx->myDoc) { xmlFreeDoc(ctx->myDoc); }
        xmlFreeParserCtxt(ctx);
      };
      unique_ptr<xmlParserCtxt, decltype(free_parser)> parser(xmlCreatePushParserCtxt(&sax, nullptr, nullptr, 0, filename.c_str()), free_parser);
      if (!parser) { throw std::bad_alloc(); }
      parser->_private = &state;
      parser->userData = parser.get();
      // Keep diagnostics local to this parser, including namespace errors.
      xmlCtxtSetErrorHandler(parser.get(), LibxmlState::error, parser.get());
      xmlCtxtUseOptions(parser.get(), XML_PARSE_NONET | XML_PARSE_NO_XXE | XML_PARSE_NOENT);
      if (!encoding.empty())
      {
        if (xmlSwitchEncodingName(parser.get(), encoding.c_str()) != 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "Unsupported XML encoding: " + encoding);
        }
      }
      char buffer[65536];
      int status = 0;
      while (!state.exception)
      {
        const auto count = read(reinterpret_cast<unsigned char*>(buffer), sizeof(buffer));
        status = xmlParseChunk(parser.get(), buffer, static_cast<int>(count), count == 0);
        if (status || count == 0) { break; }
      }
      if (state.exception)
      {
        try { std::rethrow_exception(state.exception); }
        catch (const XMLHandler::EndParsingSoftly&) { return; }
      }
      if (status)
      {
        std::string message = "libxml2 rejected XML input";
        if (const auto* error = xmlCtxtGetLastError(parser.get()))
        {
          if (error->message) { message += ": " + std::string(error->message); }
          message += " (line " + std::to_string(error->line) + ", column " + std::to_string(error->int2) + ")";
        }
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, message);
      }
    }

    void XMLFile::parse_(const std::string& filename, XMLHandler* handler, bool namespace_aware)
    {
      XMLCleaner_ clean(handler);
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      std::ifstream file(std::filesystem::u8path(filename), std::ios::binary);
      if (!file)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "Cannot open XML file");
      }
      char magic[2]{};
      file.read(magic, sizeof(magic));
      const std::string header(magic, static_cast<Size>(file.gcount()));
      file.clear();
      file.seekg(0);
      if (header == "BZ" || header == "PK" || header == std::string("\x1f\x8b", 2))
      {
        file.close();
        std::unique_ptr<InputStream> stream(CompressedInputSource(filename, header).makeStream());
        if (!stream)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "Cannot open compressed XML stream");
        }
        parse([&](unsigned char* target, Size capacity) { return stream->readBytes(target, capacity); }, handler, enforced_encoding_, filename, namespace_aware);
      }
      else
      {
        parse([&](unsigned char* target, Size capacity)
        {
          file.read(reinterpret_cast<char*>(target), static_cast<std::streamsize>(capacity));
          if (file.bad())
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "Error reading XML file");
          }
          return static_cast<Size>(file.gcount());
        }, handler, enforced_encoding_, filename, namespace_aware);
      }
    }

    void XMLFile::parseBuffer_(const std::string& buffer, XMLHandler* handler, bool namespace_aware)
    {
      XMLCleaner_ clean(handler);
      Size offset = 0;
      parse([&](unsigned char* target, Size capacity)
      {
        const Size count = std::min(capacity, buffer.size() - offset);
        std::memcpy(target, buffer.data() + offset, count);
        offset += count;
        return count;
      }, handler, enforced_encoding_, "inMemory", namespace_aware);
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
