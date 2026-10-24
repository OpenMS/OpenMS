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
#include <libxml/encoding.h>
#include <cstdint>
#include <string_view>
#include <vector>
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
#include <condition_variable>
#include <deque>
#include <mutex>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace OpenMS::Internal
{

    /**
      @brief Reads a decompressing stream ahead on a helper thread

      Inflating gzip/bzip2 input costs about as much as parsing it. Reading a few
      blocks ahead lets both run concurrently. Stream exceptions are rethrown by read().
    */
    class PrefetchingReader_
    {
public:
      explicit PrefetchingReader_(InputStream& stream) :
        stream_(stream),
        worker_([this] { run_(); })
      {
      }

      ~PrefetchingReader_()
      {
        {
          std::lock_guard<std::mutex> lock(mutex_);
          stop_ = true;
        }
        changed_.notify_all();
        worker_.join();
      }

      PrefetchingReader_(const PrefetchingReader_&) = delete;
      PrefetchingReader_& operator=(const PrefetchingReader_&) = delete;

      /// Copies up to @p capacity bytes; returns 0 at the end of the stream.
      Size read(unsigned char* target, Size capacity)
      {
        Size copied = 0;
        while (copied < capacity && (offset_ < current_.size || next_()))
        {
          const Size count = std::min(capacity - copied, current_.size - offset_);
          std::memcpy(target + copied, current_.data.get() + offset_, count);
          offset_ += count;
          copied += count;
        }
        return copied;
      }

private:
      static constexpr Size block_size_ = 256 * 1024;
      static constexpr Size max_blocks_ = 4;

      struct Block
      {
        std::unique_ptr<unsigned char[]> data;
        Size size = 0;
      };

      /// Moves on to the next filled block; false at the end of the stream.
      bool next_()
      {
        std::unique_lock<std::mutex> lock(mutex_);
        if (current_.data) { spare_.push_back(std::move(current_)); }
        current_ = Block();
        offset_ = 0;
        changed_.notify_all();
        changed_.wait(lock, [this] { return !filled_.empty() || done_; });
        if (filled_.empty())
        {
          if (error_) { std::rethrow_exception(error_); }
          return false;
        }
        current_ = std::move(filled_.front());
        filled_.pop_front();
        changed_.notify_all();
        return true;
      }

      void run_()
      {
        try
        {
          while (true)
          {
            Block block;
            {
              std::unique_lock<std::mutex> lock(mutex_);
              changed_.wait(lock, [this] { return stop_ || filled_.size() < max_blocks_; });
              if (stop_) { return; }
              if (!spare_.empty())
              {
                block = std::move(spare_.back());
                spare_.pop_back();
              }
            }
            if (!block.data) { block.data.reset(new unsigned char[block_size_]); }
            block.size = stream_.readBytes(block.data.get(), block_size_);
            std::lock_guard<std::mutex> lock(mutex_);
            if (block.size == 0) { done_ = true; }
            else { filled_.push_back(std::move(block)); }
            changed_.notify_all();
            if (done_) { return; }
          }
        }
        catch (...)
        {
          std::lock_guard<std::mutex> lock(mutex_);
          error_ = std::current_exception();
          done_ = true;
          changed_.notify_all();
        }
      }

      InputStream& stream_;
      std::mutex mutex_;
      std::condition_variable changed_;
      std::deque<Block> filled_;
      std::vector<Block> spare_;
      bool done_ = false;
      bool stop_ = false;
      std::exception_ptr error_;
      // consumer side only
      Block current_;
      Size offset_ = 0;
      // started last, after all members it uses are initialized
      std::thread worker_;
    };

    /**
      @brief Leaves one core to a PrefetchingReader_ thread while it exists

      Handlers decode arrays in OpenMP regions between parsing steps. Idle OpenMP
      workers spin-wait for the next region and, using all cores, would starve the
      decompression thread. Affects only parallel regions started by this thread.
    */
    class OpenMPThreadReservation_
    {
public:
      OpenMPThreadReservation_()
      {
#ifdef _OPENMP
        previous_ = omp_get_max_threads();
        if (previous_ > 1) { omp_set_num_threads(previous_ - 1); }
#endif
      }

      ~OpenMPThreadReservation_()
      {
#ifdef _OPENMP
        if (previous_ > 1) { omp_set_num_threads(previous_); }
#endif
      }

      OpenMPThreadReservation_(const OpenMPThreadReservation_&) = delete;
      OpenMPThreadReservation_& operator=(const OpenMPThreadReservation_&) = delete;

private:
      int previous_ = 1;
    };

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

    namespace
    {
      bool isASCIIWord(const unsigned char* data)
      {
        std::uint64_t word;
        std::memcpy(&word, data, sizeof(word));
        return (word & 0x8080808080808080ull) == 0;
      }

      /// End of the ASCII run starting at @p from
      Size asciiRunEnd(const unsigned char* data, Size from, Size size)
      {
        while (from + 8 <= size && isASCIIWord(data + from)) from += 8;
        while (from < size && data[from] < 0x80) ++from;
        return from;
      }

      /**
        @brief Converts an ISO-8859-1 chunk to UTF-8

        Returns @p data itself if the chunk is ASCII (as mzML usually is), otherwise
        the converted bytes in @p storage. Latin-1 is a single-byte encoding, so
        chunks convert independently.
      */
      const unsigned char* latin1ToUTF8(const unsigned char* data, Size& size, std::vector<unsigned char>& storage)
      {
        Size in = asciiRunEnd(data, 0, size);
        if (in == size) return data;
        storage.resize(2 * size);
        unsigned char* out = storage.data();
        std::memcpy(out, data, in);
        out += in;
        while (in < size)
        {
          const unsigned char c = data[in++];
          *out++ = static_cast<unsigned char>(0xc0 | (c >> 6));
          *out++ = static_cast<unsigned char>(0x80 | (c & 0x3f));
          const Size run = asciiRunEnd(data, in, size);
          std::memcpy(out, data + in, run - in);
          out += run - in;
          in = run;
        }
        size = static_cast<Size>(out - storage.data());
        return storage.data();
      }

      /// The encoding named by the XML declaration at the start of @p data, or empty if there is none.
      std::string declaredEncoding(const unsigned char* data, Size size)
      {
        const std::string_view text(reinterpret_cast<const char*>(data), std::min<Size>(size, 1024));
        if (text.substr(0, 6) != "<?xml " && text.substr(0, 6) != "<?xml\t" && text.substr(0, 6) != "<?xml\n" && text.substr(0, 6) != "<?xml\r") return {};
        const auto declaration = text.substr(0, text.find("?>"));
        if (declaration.size() == text.size()) return {};
        auto pos = declaration.find("encoding");
        if (pos == std::string_view::npos) return {};
        pos = declaration.find_first_not_of(" \t\r\n", pos + 8);
        if (pos == std::string_view::npos || declaration[pos] != '=') return {};
        pos = declaration.find_first_not_of(" \t\r\n", pos + 1);
        if (pos == std::string_view::npos || (declaration[pos] != '"' && declaration[pos] != '\'')) return {};
        const auto close = declaration.find(declaration[pos], pos + 1);
        if (close == std::string_view::npos) return {};
        return std::string(declaration.substr(pos + 1, close - pos - 1));
      }

      /// Whether libxml2 decodes @p encoding with its built-in ISO-8859-1 decoder
      bool isLatin1(const std::string& encoding)
      {
        return !encoding.empty() && xmlParseCharEncoding(encoding.c_str()) == XML_CHAR_ENCODING_8859_1;
      }

      /// Converts parser UTF-8 to UTF-16 into @p result (reusing its capacity).
      const char16_t* utf16(std::u16string& result, const xmlChar* text, Size length)
      {
        result.clear();
        if (!text) return result.c_str();
        const unsigned char* cursor = text;
        const unsigned char* end = text + length;
        const unsigned char* ascii = text + asciiRunEnd(text, 0, length);
        result.assign(cursor, ascii);
        cursor = ascii;
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
        return result.c_str();
      }

      const char16_t* utf16(std::u16string& result, const xmlChar* text)
      {
        return utf16(result, text, text ? static_cast<Size>(xmlStrlen(text)) : 0);
      }

      std::string_view view(const xmlChar* text)
      {
        return text ? std::string_view(reinterpret_cast<const char*>(text)) : std::string_view();
      }
    }

    // C callbacks capture every exception: no C++ exception crosses libxml2 frames.
    struct LibxmlState
    {
      XMLHandler* handler;
      std::exception_ptr exception;
      // Reused across callbacks so that a parse does not allocate per element.
      XMLAttributes::NativeStore attributes;
      std::string name8;
      std::u16string name;
      std::u16string text;

      static LibxmlState& state(void* ctx)
      {
        return *static_cast<LibxmlState*>(static_cast<xmlParserCtxt*>(ctx)->_private);
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
          s.attributes.clear();
          if (attrs)
          {
            for (Size i = 0; attrs[i]; i += 2) s.attributes.add({}, view(attrs[i]), view(attrs[i + 1]));
          }
          const XMLAttributes attributes(s.attributes);
          if (s.handler->tryOnStartElementUTF8(view(name), attributes)) { return; }
          s.handler->onStartElement(utf16(s.name, name), attributes);
        });
      }
      static void endLegacy(void* ctx, const xmlChar* name)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          if (s.handler->tryOnEndElementUTF8(view(name))) { return; }
          s.handler->onEndElement(utf16(s.name, name));
        });
      }
      std::string_view qualifiedNameUTF8(const xmlChar* local, const xmlChar* prefix)
      {
        if (!prefix) return view(local);
        name8.assign(view(prefix)).append(1, ':').append(view(local));
        return name8;
      }
      const char16_t* qualifiedName(const xmlChar* local, const xmlChar* prefix)
      {
        if (!prefix) return utf16(name, local);
        utf16(text, local); // no character data is pending inside a tag
        utf16(name, prefix);
        name += u':';
        name += text;
        return name.c_str();
      }
      static void start(void* ctx, const xmlChar* local, const xmlChar* prefix,
                        const xmlChar*, int namespace_count, const xmlChar** namespaces,
                        int attribute_count, int, const xmlChar** attrs)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          s.attributes.clear();
          // Existing handlers consume qualified names and can inspect xmlns declarations.
          for (int i = 0; i < namespace_count; ++i)
          {
            if (namespaces[2 * i]) { s.attributes.add("xmlns", view(namespaces[2 * i]), view(namespaces[2 * i + 1])); }
            else { s.attributes.add({}, "xmlns", view(namespaces[2 * i + 1])); }
          }
          for (int i = 0; i < attribute_count; ++i)
          {
            const xmlChar** attr = attrs + 5 * i;
            // SAX2 values are bounded slices, not null-terminated strings.
            s.attributes.add(view(attr[1]), view(attr[0]),
                             std::string_view(reinterpret_cast<const char*>(attr[3]), static_cast<Size>(attr[4] - attr[3])));
          }
          const XMLAttributes attributes(s.attributes);
          if (s.handler->tryOnStartElementUTF8(s.qualifiedNameUTF8(local, prefix), attributes)) { return; }
          s.handler->onStartElement(s.qualifiedName(local, prefix), attributes);
        });
      }
      static void end(void* ctx, const xmlChar* local, const xmlChar* prefix, const xmlChar*)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          if (s.handler->tryOnEndElementUTF8(s.qualifiedNameUTF8(local, prefix))) { return; }
          s.handler->onEndElement(s.qualifiedName(local, prefix));
        });
      }
      static void characters(void* ctx, const xmlChar* text, int length)
      {
        callback(ctx, [&](LibxmlState& s)
        {
          if (s.handler->tryOnCharactersUTF8(reinterpret_cast<const char*>(text), static_cast<Size>(length))) { return; }
          utf16(s.text, text, static_cast<Size>(length));
          s.handler->onCharacters(s.text.c_str(), s.text.size());
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
      LibxmlState state{handler, {}, {}, {}, {}, {}};
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

      unsigned char buffer[65536];
      Size count = read(buffer, sizeof(buffer));
      // OpenMS writes ISO-8859-1 declarations. libxml2 would decode those byte by byte;
      // instead, pass ASCII chunks through unchanged and convert the others here.
      const bool latin1 = isLatin1(encoding.empty() ? declaredEncoding(buffer, count) : encoding);
      if (latin1)
      {
        // An explicit encoding takes precedence over the declaration.
        if (xmlSwitchEncodingName(parser.get(), "UTF-8") != 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "Cannot switch XML encoding");
        }
      }
      else if (!encoding.empty())
      {
        if (xmlSwitchEncodingName(parser.get(), encoding.c_str()) != 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "Unsupported XML encoding: " + encoding);
        }
      }
      std::vector<unsigned char> converted;
      int status = 0;
      while (!state.exception)
      {
        Size size = count;
        const unsigned char* chunk = latin1 ? latin1ToUTF8(buffer, size, converted) : buffer;
        status = xmlParseChunk(parser.get(), reinterpret_cast<const char*>(chunk), static_cast<int>(size), count == 0);
        if (status || count == 0) { break; }
        count = read(buffer, sizeof(buffer));
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
        OpenMPThreadReservation_ reservation;
        PrefetchingReader_ reader(*stream);
        parse([&](unsigned char* target, Size capacity) { return reader.read(target, capacity); }, handler, enforced_encoding_, filename, namespace_aware);
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
