// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/CONCEPT/Init.h>

#include <OpenMS/SYSTEM/File.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>

#include <libxml/parser.h>
#include <libxml/xmlreader.h>
#include <memory>
#include <limits>
#include <cstdint>
#include <algorithm>
#include <exception>

namespace OpenMS
{

  namespace IndexedMzMLUtils
  {
    std::streampos stringToStreampos(const std::string& s)
    {
      // Try to cast the string to a type that can hold the integer value
      // stored in the std::streampos type (which can vary from 32 to 128 bits
      // depending on the system).
      //
      // long long has a minimal size of 64 bits and can address a range up to
      // 16 Exbibit (or 2 Exbibyte), we can hopefully expect our files to be
      // smaller than an Exabyte and should be safe.
      std::streampos res;
      try
      {

        res = boost::lexical_cast< unsigned long long >(s);
        // TESTING CAST: res = (int) res; // use this only when emulating 32 bit systems

      }
      catch (boost::bad_lexical_cast&)
      {
        std::cerr << "Trying to convert corrupted / unreadable value to std::streampos : " << s << std::endl;
        std::cerr << "This can also happen if the value exceeds 63 bits, please check your input." << std::endl;
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            std::string("Could not convert string '") + s + "' to a 64 bit integer.");
      }

      // Check if the value can fit into std::streampos
      if ( std::abs( boost::lexical_cast< long double >(s) - res) > 0.1)
      {
        std::cerr << "Your system may not support addressing a file of this size,"
          << " only addresses that fit into a " << sizeof(std::streamsize)*8 <<
          " bit integer are supported on your system." << std::endl;
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            std::string("Could not convert string '") + s + "' to an integer on your system.");
      }

      return res;
    }
  }

  int IndexedMzMLDecoder::parseOffsets(const std::string& filename, std::streampos indexoffset, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open())
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    f.seekg(0, std::ios::end);
    const std::streampos length = f.tellg();
    if (!f || length < 0 || indexoffset < 0 || indexoffset > length)
    {
      std::cerr << "IndexedMzMLDecoder::parseOffsets Error: Offset was " << indexoffset
                << " (not between 0 and " << length << ").\n";
      return -1;
    }

    // Read straight into the final XML buffer: no raw allocation or full-tail copy.
    std::string xml = "<indexedmzML>";
    const Size prefix_size = xml.size();
    const std::streamoff tail_length = length - indexoffset;
    if (tail_length == 0 || tail_length > std::numeric_limits<std::streamsize>::max() ||
        static_cast<std::uintmax_t>(tail_length) > xml.max_size() - prefix_size - 1)
    {
      return -1;
    }
    const Size tail_size = static_cast<Size>(tail_length);
    xml.resize(prefix_size + tail_size + 1);
    f.seekg(indexoffset);
    if (!f) { return -1; }
    f.read(xml.data() + prefix_size, static_cast<std::streamsize>(tail_length));
    if (!f || f.gcount() != tail_length) { return -1; }
    xml.back() = '\n';
    return domParseIndexedEnd_(xml, spectra_offsets, chromatograms_offsets);
  }

  std::streampos IndexedMzMLDecoder::findIndexListOffset(const std::string& filename, int buffersize)
  {
    if (buffersize <= 0) { return -1; }
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open())
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    f.seekg(0, std::ios::end);
    const std::streampos length = f.tellg();
    if (!f || length < 0) { return -1; }
    const std::streamoff read_length = std::min(static_cast<std::streamoff>(length), static_cast<std::streamoff>(buffersize));
    std::string buffer;
    if (read_length > std::numeric_limits<std::streamsize>::max() ||
        static_cast<std::uintmax_t>(read_length) > buffer.max_size()) { return -1; }
    buffer.resize(static_cast<Size>(read_length));
    f.seekg(-read_length, std::ios::end);
    if (!f) { return -1; }
    f.read(buffer.data(), static_cast<std::streamsize>(read_length));
    if (!f || f.gcount() != read_length) { return -1; }

    // The same namespace/whitespace-tolerant pattern is shared by independent
    // readers; C++ initializes it once and regex_search never mutates it.
    static const boost::regex listoffset_rx(R"(<[^>/]*indexListOffset\s*>\s*(\d*))");
    boost::smatch matches;
    if (!boost::regex_search(buffer, matches, listoffset_rx) || matches[1].length() == 0)
    {
      std::cerr << "IndexedMzMLDecoder::findIndexListOffset Error: Could not find element indexListOffset in the last "
                << read_length << " bytes. Maybe this is not an indexedMzML.\n";
      return -1;
    }
    return OpenMS::IndexedMzMLUtils::stringToStreampos(matches[1].str());
  }

  int IndexedMzMLDecoder::domParseIndexedEnd_(const std::string& in, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {

    if (in.size() > static_cast<Size>(std::numeric_limits<int>::max())) { return -1; }
    // Keep the existing DOM path for small footers. On Linux, replacing even
    // small DOMs changes allocator reuse and can slow large-spectrum reads.
    // Bound that fast path to 2 MiB; larger indices use bounded-node traversal.
    constexpr Size small_footer_limit = 2 * 1024 * 1024;
    if (in.size() <= small_footer_limit)
    {
      std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)> document(
        xmlReadMemory(in.data(), static_cast<int>(in.size()), "index.xml", nullptr,
                      XML_PARSE_NONET | XML_PARSE_NO_XXE | XML_PARSE_NOENT | XML_PARSE_NOERROR | XML_PARSE_NOWARNING), xmlFreeDoc);
      xmlNode* root = document ? xmlDocGetRootElement(document.get()) : nullptr;
      if (!root) { return -1; }
      xmlNode* list = nullptr;
      for (xmlNode* node = root->children; node; node = node->next)
      {
        if (node->type == XML_ELEMENT_NODE && xmlStrEqual(node->name, BAD_CAST "indexList"))
        {
          if (list) { return -1; }
          list = node;
        }
      }
      if (!list) { return -1; }
      auto property = [](xmlNode* node, const char* key)
      {
        xmlChar* value = xmlGetProp(node, BAD_CAST key);
        std::string result = value ? reinterpret_cast<const char*>(value) : "";
        xmlFree(value);
        return result;
      };
      OffsetVector spectra, chromatograms;
      bool found_spectra = false, found_chromatograms = false, last_was_spectra = false;
      for (xmlNode* index = list->children; index; index = index->next)
      {
        if (index->type != XML_ELEMENT_NODE) { continue; }
        OffsetVector result;
        for (xmlNode* offset = index->children; offset; offset = offset->next)
        {
          if (offset->type != XML_ELEMENT_NODE) { continue; }
          xmlChar* content = xmlNodeGetContent(offset);
          std::string text = content ? reinterpret_cast<const char*>(content) : "";
          xmlFree(content);
          result.emplace_back(property(offset, "idRef"), IndexedMzMLUtils::stringToStreampos(text));
        }
        const auto name = property(index, "name");
        if (name == "spectrum") { spectra = std::move(result); found_spectra = true; last_was_spectra = true; }
        else if (name == "chromatogram") { chromatograms = std::move(result); found_chromatograms = true; last_was_spectra = false; }
        else { return -1; }
      }
      if (&spectra_offsets == &chromatograms_offsets && found_spectra && last_was_spectra)
      {
        spectra_offsets = std::move(spectra);
      }
      else
      {
        if (found_spectra) { spectra_offsets = std::move(spectra); }
        if (found_chromatograms) { chromatograms_offsets = std::move(chromatograms); }
      }
      return 0;
    }
    // A DOM for a large index leaves freed nodes in each worker's allocator
    // arena. Read forward and expand only one offset at a time instead.
    std::unique_ptr<xmlTextReader, decltype(&xmlFreeTextReader)> reader(
      xmlReaderForMemory(in.data(), static_cast<int>(in.size()), "index.xml", nullptr,
                         XML_PARSE_NONET | XML_PARSE_NO_XXE | XML_PARSE_NOENT | XML_PARSE_NOERROR | XML_PARSE_NOWARNING), xmlFreeTextReader);
    if (!reader) { return -1; }
    auto property = [](xmlNode* node, const char* key)
    {
      xmlChar* value = xmlGetProp(node, BAD_CAST key);
      std::string result = value ? reinterpret_cast<const char*>(value) : "";
      xmlFree(value);
      return result;
    };
    OffsetVector spectra, chromatograms, offsets;
    std::string index_name;
    bool found_list = false, in_list = false, found_spectra = false, found_chromatograms = false, last_was_spectra = false;
    auto finish_index = [&]()
    {
      if (index_name == "spectrum")
      {
        spectra = std::move(offsets);
        found_spectra = true;
        last_was_spectra = true;
      }
      else if (index_name == "chromatogram")
      {
        chromatograms = std::move(offsets);
        found_chromatograms = true;
        last_was_spectra = false;
      }
      else { return false; }
      offsets.clear();
      return true;
    };
    std::exception_ptr conversion_error;
    int status;
    while ((status = xmlTextReaderRead(reader.get())) == 1)
    {
      const int depth = xmlTextReaderDepth(reader.get());
      const int type = xmlTextReaderNodeType(reader.get());
      if (type == XML_READER_TYPE_ELEMENT)
      {
        if (depth == 1 && xmlStrEqual(xmlTextReaderConstLocalName(reader.get()), BAD_CAST "indexList"))
        {
          if (found_list) { return -1; }
          found_list = true;
          in_list = !xmlTextReaderIsEmptyElement(reader.get());
        }
        else if (conversion_error) { continue; }
        else if (in_list && depth == 2)
        {
          index_name = property(xmlTextReaderCurrentNode(reader.get()), "name");
          if (xmlTextReaderIsEmptyElement(reader.get()) && !finish_index()) { return -1; }
        }
        else if (in_list && depth == 3)
        {
          // Expand this offset's descendants to preserve xmlNodeGetContent
          // semantics (entities, CDATA and nested text), never the whole index.
          xmlNode* offset = xmlTextReaderExpand(reader.get());
          if (!offset) { return -1; }
          xmlChar* content = xmlNodeGetContent(offset);
          std::string text = content ? reinterpret_cast<const char*>(content) : "";
          xmlFree(content);
          try { offsets.emplace_back(property(offset, "idRef"), IndexedMzMLUtils::stringToStreampos(text)); }
          catch (const Exception::ConversionError&)
          {
            // As with the old whole-document parse, malformed XML takes
            // precedence over an invalid offset earlier in the document.
            conversion_error = std::current_exception();
          }
        }
      }
      else if (type == XML_READER_TYPE_END_ELEMENT && in_list && !conversion_error)
      {
        if (depth == 2 && !finish_index()) { return -1; }
        else if (depth == 1) { in_list = false; }
      }
    }
    if (status < 0 || !found_list) { return -1; }
    if (conversion_error) { std::rethrow_exception(conversion_error); }
    if (&spectra_offsets == &chromatograms_offsets && found_spectra && last_was_spectra)
    {
      spectra_offsets = std::move(spectra);
    }
    else
    {
      if (found_spectra) { spectra_offsets = std::move(spectra); }
      if (found_chromatograms) { chromatograms_offsets = std::move(chromatograms); }
    }
    return 0;
  }

}
