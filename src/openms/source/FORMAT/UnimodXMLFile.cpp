// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/UnimodXMLHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <limits>
#include <filesystem>
#include <string_view>
#include <thread>
#include <vector>

using namespace std;

namespace OpenMS
{

  UnimodXMLFile::UnimodXMLFile() :
    Internal::XMLFile()
  {

  }

  UnimodXMLFile::~UnimodXMLFile() = default;

  void UnimodXMLFile::load(const std::string& filename, vector<ResidueModification*> & modifications)
  {
    std::string file = File::find(filename);

    if (loadInSlices_(file, modifications))
    {
      return;
    }
    Internal::UnimodXMLHandler handler(modifications, file);
    parse_(file, &handler);
  }

  namespace
  {
    /// A start or end tag in the text of a file without comments, CDATA sections and processing instructions
    struct Tag
    {
      std::string_view name; ///< with its prefix, e.g. "umod:mod"
      std::string_view text; ///< from '<' to '>'
      bool closing = false;
      bool self_closing = false;
    };

    /// Reads the tag that starts at the '<' at @p p and ends before @p limit ('>' in quoted attribute values is skipped)
    bool readTag(std::string_view text, size_t p, size_t limit, Tag& tag)
    {
      size_t i = p + 1;
      tag.closing = i < limit && text[i] == '/';
      if (tag.closing) ++i;
      size_t n = i;
      while (n < limit && !std::isspace(static_cast<unsigned char>(text[n])) && text[n] != '>' && text[n] != '/') ++n;
      tag.name = text.substr(i, n - i);
      char quote = 0;
      for (size_t k = n; k < limit; ++k)
      {
        const char c = text[k];
        if (quote != 0) { if (c == quote) quote = 0; }
        else if (c == '"' || c == '\'') quote = c;
        else if (c == '<') return false;
        else if (c == '>')
        {
          tag.text = text.substr(p, k + 1 - p);
          tag.self_closing = text[k - 1] == '/';
          return !tag.name.empty();
        }
      }
      return false;
    }

    /// The value of attribute @p name of a start tag as written (entities not expanded); false if it is not there
    bool attributeValue(std::string_view tag_text, std::string_view name, std::string_view& value)
    {
      size_t i = 1;
      while (i < tag_text.size() && !std::isspace(static_cast<unsigned char>(tag_text[i])) && tag_text[i] != '>' && tag_text[i] != '/') ++i;
      while (true)
      {
        while (i < tag_text.size() && std::isspace(static_cast<unsigned char>(tag_text[i]))) ++i;
        if (i >= tag_text.size() || tag_text[i] == '>' || tag_text[i] == '/') return false;
        const size_t name_begin = i;
        while (i < tag_text.size() && !std::isspace(static_cast<unsigned char>(tag_text[i])) && tag_text[i] != '=') ++i;
        const std::string_view attribute = tag_text.substr(name_begin, i - name_begin);
        while (i < tag_text.size() && std::isspace(static_cast<unsigned char>(tag_text[i]))) ++i;
        if (i >= tag_text.size() || tag_text[i] != '=') return false;
        ++i;
        while (i < tag_text.size() && std::isspace(static_cast<unsigned char>(tag_text[i]))) ++i;
        if (i >= tag_text.size() || (tag_text[i] != '"' && tag_text[i] != '\'')) return false;
        const char quote = tag_text[i++];
        const size_t value_end = tag_text.find(quote, i);
        if (value_end == std::string_view::npos) return false;
        if (attribute == name)
        {
          value = tag_text.substr(i, value_end - i);
          return true;
        }
        i = value_end + 1;
      }
    }

    /// The elements UnimodXMLHandler acts on, with the prefix of unimod.xml; the same names without it (which the
    /// handler accepts too) make a file unsuitable for slices
    enum class HandlerTag { NONE, MOD, SPECIFICITY, NEUTRAL_LOSS, DELTA, ELEMENT, UNPREFIXED };

    HandlerTag handlerTag(std::string_view name)
    {
      const bool prefixed = name.substr(0, 5) == "umod:";
      const std::string_view local = prefixed ? name.substr(5) : name;
      HandlerTag tag = HandlerTag::NONE;
      if (local == "mod") tag = HandlerTag::MOD;
      else if (local == "specificity") tag = HandlerTag::SPECIFICITY;
      else if (local == "NeutralLoss") tag = HandlerTag::NEUTRAL_LOSS;
      else if (local == "delta") tag = HandlerTag::DELTA;
      else if (local == "element") tag = HandlerTag::ELEMENT;
      return (tag != HandlerTag::NONE && !prefixed) ? HandlerTag::UNPREFIXED : tag;
    }

    /// Follows @p tag on the stack of @p open elements; false for markup other than tags (comments, CDATA sections,
    /// processing instructions, declarations) and for an end tag that does not close the innermost open element
    bool balance(const Tag& tag, std::vector<std::string_view>& open)
    {
      if (tag.name[0] == '!' || tag.name[0] == '?') return false;
      if (tag.closing)
      {
        if (open.empty() || open.back() != tag.name) return false;
        open.pop_back();
      }
      else if (!tag.self_closing)
      {
        open.push_back(tag.name);
      }
      return true;
    }

    /**
      Is [@p from, @p to) of the \<umod:modifications\> section a sequence of complete \<umod:mod\> elements that the
      handler reads the same way with a fresh handler as after the elements before it? That is: plain tags only (no
      comments, CDATA sections, processing instructions, which could hold text like tags), \<umod:specificity\> only in
      \<umod:mod\>, \<umod:NeutralLoss\> only in \<umod:specificity\> (the neutral losses of an open specificity are what
      the handler carries from one element to the next otherwise), \<umod:delta\> and \<umod:element\> only in
      \<umod:mod\>, and every specificity with a position that the handler knows (it reports others).
    */
    bool sliceable(std::string_view text, size_t from, size_t to)
    {
      int in_mod = 0, in_specificity = 0;
      std::vector<std::string_view> open; // all elements: a slice is a sequence of complete elements
      Tag tag;
      for (size_t p = text.find('<', from); p < to; p = text.find('<', p + 1))
      {
        if (!readTag(text, p, to, tag) || !balance(tag, open)) return false;
        p += tag.text.size() - 1;
        switch (handlerTag(tag.name))
        {
          case HandlerTag::NONE: break;
          case HandlerTag::UNPREFIXED: return false;
          case HandlerTag::MOD:
            if (in_specificity != 0 || in_mod != (tag.closing ? 1 : 0)) return false;
            in_mod = (tag.closing || tag.self_closing) ? 0 : 1;
            break;
          case HandlerTag::SPECIFICITY:
            if (in_mod != 1 || in_specificity != (tag.closing ? 1 : 0)) return false;
            if (!tag.closing)
            {
              std::string_view pos;
              if (!attributeValue(tag.text, "position", pos)) return false;
              if (pos != "Anywhere" && pos != "Protein N-term" && pos != "Protein C-term" && pos != "Any C-term" && pos != "Any N-term") return false;
            }
            in_specificity = (tag.closing || tag.self_closing) ? 0 : 1;
            break;
          case HandlerTag::NEUTRAL_LOSS:
            if (in_specificity != 1) return false;
            break;
          case HandlerTag::DELTA:
          case HandlerTag::ELEMENT:
            if (in_mod != 1) return false;
            break;
        }
      }
      return in_mod == 0 && in_specificity == 0 && open.empty();
    }

    /**
      Does [@p from, @p to) hold plain tags only, none that the handler acts on except \<umod:element\> if
      @p elements_allowed (outside \<umod:delta\>/\<umod:NeutralLoss\> the handler ignores those), and does it leave
      the elements @p open (on entry: open before it) open at its end?
    */
    bool inert(std::string_view text, size_t from, size_t to, bool elements_allowed, std::vector<std::string_view>& open)
    {
      Tag tag;
      for (size_t p = text.find('<', from); p < to; p = text.find('<', p + 1))
      {
        if (!readTag(text, p, to, tag) || !balance(tag, open)) return false;
        p += tag.text.size() - 1;
        const HandlerTag h = handlerTag(tag.name);
        if (h != HandlerTag::NONE && !(elements_allowed && h == HandlerTag::ELEMENT)) return false;
      }
      return true;
    }

    /// Where the root element \<umod:unimod\> starts, if the file before it is an XML declaration, comments and
    /// white space only (no document type declaration: entities could change the meaning of the text); else npos
    size_t rootStart(std::string_view text)
    {
      size_t p = text.substr(0, 3) == "\xEF\xBB\xBF" ? 3 : 0;
      if (text.substr(p, 6) == "<?xml ")
      {
        const size_t declaration_end = text.find("?>", p);
        if (declaration_end == std::string_view::npos) return std::string_view::npos;
        // XML 1.0 only: XML 1.1 turns NEL and LSEP into line ends, which the byte scanner below would read as name
        // characters
        std::string_view version;
        if (attributeValue(text.substr(p, declaration_end - p), "version", version) && version != "1.0") return std::string_view::npos;
        // the tags are read from the bytes: UTF-8 (the default) only, not e.g. UTF-7, which can encode '<'
        std::string_view encoding;
        if (attributeValue(text.substr(p, declaration_end - p) , "encoding", encoding))
        {
          std::string e(encoding);
          std::transform(e.begin(), e.end(), e.begin(), [](unsigned char c) { return char(std::tolower(c)); });
          if (e != "utf-8" && e != "utf8") return std::string_view::npos;
        }
        p = declaration_end + 2;
      }
      while (p < text.size())
      {
        if (std::isspace(static_cast<unsigned char>(text[p]))) { ++p; continue; }
        if (text.substr(p, 4) == "<!--")
        {
          p = text.find("-->", p + 4);
          if (p == std::string_view::npos) return std::string_view::npos;
          p += 3;
          continue;
        }
        Tag tag;
        if (text[p] == '<' && readTag(text, p, text.size(), tag) && !tag.closing && !tag.self_closing && tag.name == "umod:unimod") return p;
        return std::string_view::npos;
      }
      return std::string_view::npos;
    }
  }

  bool UnimodXMLFile::loadInSlices_(const std::string & file, vector<ResidueModification*> & modifications)
  {
    // unimod.xml (2.5 MB) takes most of the start-up of a tool that needs the ModificationsDB
    const Size max_slices = 8;
    const Size min_slice_bytes = 256 * 1024;
    const Size n_threads = std::thread::hardware_concurrency();
    if (n_threads < 2 || !enforced_encoding_.empty()) return false;

    // the slices, each parsed with its own handler; deleted unless all are parsed and moved to 'modifications'
    std::vector<std::vector<ResidueModification*>> parts;
    struct PartsOwner
    {
      std::vector<std::vector<ResidueModification*>>& parts;
      ~PartsOwner() { for (auto& part : parts) for (ResidueModification* m : part) delete m; }
    } owner{parts};
    try
    {
      std::error_code ec;
      if (!std::filesystem::is_regular_file(file, ec)) return false;
      const auto file_size = std::filesystem::file_size(file, ec);
      if (ec || file_size < 2 * min_slice_bytes || file_size > uintmax_t(std::numeric_limits<std::streamsize>::max()) || file_size >= std::string().max_size()) return false;
      std::string text;
      {
        std::ifstream in(file, std::ios::binary);
        if (!in) return false;
        text.resize(size_t(file_size));
        if (!in.read(text.data(), std::streamsize(file_size)) || in.peek() != std::ifstream::traits_type::eof()) return false;
      }
      const std::string_view view(text);
      const size_t root = rootStart(view);
      if (root == std::string_view::npos) return false;
      const std::string open_tag = "<umod:modifications>", close_tag = "</umod:modifications>", mod_end = "</umod:mod>";
      const size_t open = text.find(open_tag, root);
      if (open == std::string::npos) return false;
      const size_t begin = open + open_tag.size();
      const size_t end = text.find(close_tag, begin);
      if (end == std::string::npos) return false;
      // nothing the handler acts on before or after the section: it is read again with every slice
      // ... and the section is the second level of the document: <umod:unimod><umod:modifications>, closed after it
      std::vector<std::string_view> open_elements;
      if (!inert(view, root, begin, false, open_elements) ||
          open_elements != std::vector<std::string_view>{"umod:unimod", "umod:modifications"}) return false;
      if (!inert(view, end, text.size(), true, open_elements) || !open_elements.empty()) return false;
      const Size n_slices = std::min({max_slices, n_threads, Size((end - begin) / min_slice_bytes)});
      if (n_slices < 2) return false;
      std::vector<size_t> cuts{begin};
      for (Size k = 1; k < n_slices; ++k)
      {
        const size_t p = text.find(mod_end, std::max(cuts.back(), begin + (end - begin) * k / n_slices));
        if (p == std::string::npos || p + mod_end.size() >= end) break;
        cuts.push_back(p + mod_end.size());
      }
      cuts.push_back(end);
      if (cuts.size() < 3) return false;

      // Each slice on a thread of its own, whose log streams are silenced: a slice that is not sliceable() or fails
      // to parse makes the whole file parse again as before, which reports what it reports, once and in order.
      parts.resize(cuts.size() - 1);
      std::vector<char> ok(parts.size(), 0);
      auto parse_slice = [&](Size i)
      {
        try
        {
          for (Logger::LogStream* log : {&getThreadLocalLogFatal(), &getThreadLocalLogError(), &getThreadLocalLogWarn(),
                                         &getThreadLocalLogInfo(), &getThreadLocalLogDebug()})
          {
            log->removeAllStreams();
          }
          if (!sliceable(view, cuts[i], cuts[i + 1])) return;
          // everything before the section + the slice + everything after the section
          const std::string buffer = text.substr(0, begin) + text.substr(cuts[i], cuts[i + 1] - cuts[i]) + text.substr(end);
          Internal::UnimodXMLHandler handler(parts[i], file);
          UnimodXMLFile().parseBuffer_(buffer, &handler);
          ok[i] = 1;
        }
        catch (...)
        {
        }
      };
      std::vector<std::thread> threads;
      threads.reserve(parts.size());
      try
      {
        for (Size i = 0; i < parts.size(); ++i) threads.emplace_back(parse_slice, i);
      }
      catch (...)
      {
        for (auto& t : threads) t.join();
        return false;
      }
      for (auto& t : threads) t.join();
      if (std::find(ok.begin(), ok.end(), 0) != ok.end()) return false;

      Size n_mods = modifications.size();
      for (const auto& part : parts) n_mods += part.size();
      modifications.reserve(n_mods);
    }
    catch (...)
    {
      return false;
    }
    for (auto& part : parts)
    {
      modifications.insert(modifications.end(), part.begin(), part.end()); // no allocation: reserved above
      part.clear();
    }
    return true;
  }

} // namespace OpenMS
