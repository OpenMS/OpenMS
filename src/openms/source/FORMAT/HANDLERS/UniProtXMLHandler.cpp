// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/UniProtXMLHandler.h>

#include <cctype>
#include <cstdlib>
#include <utility>

namespace OpenMS::Internal
{
  namespace
  {
    /// Forward declaration: defined below xmlchToString().
    std::string xmlchToString(const XMLCh* ch);

    /// Read a Xerces attribute by name (qname), return "" if absent.
    std::string attrValue(const xercesc::Attributes& attrs, const char* name)
    {
      XMLCh* qname = xercesc::XMLString::transcode(name);
      const XMLCh* value = attrs.getValue(qname);
      xercesc::XMLString::release(&qname);
      if (value == nullptr) return std::string();
      return xmlchToString(value);
    }

    /// Strip all whitespace from a string (pretty-printed UniProt sequences contain interior newlines/tabs).
    std::string stripWhitespace(const std::string& s)
    {
      std::string out;
      out.reserve(s.size());
      for (char c : s)
      {
        if (!std::isspace(static_cast<unsigned char>(c))) out.push_back(c);
      }
      return out;
    }

    /// Manually encode a single Unicode code point into UTF-8 and append to @p out.
    void appendUtf8(std::string& out, uint32_t cp)
    {
      if (cp < 0x80)
      {
        out.push_back(static_cast<char>(cp));
      }
      else if (cp < 0x800)
      {
        out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
      }
      else if (cp < 0x10000)
      {
        out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
        out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
      }
      else
      {
        out.push_back(static_cast<char>(0xF0 | (cp >> 18)));
        out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
      }
    }

    /// Convert XMLCh* (UTF-16) to std::string in UTF-8 encoding.
    /// We can't use Xerces's @c XMLString::transcode here because it goes through
    /// the local C-runtime encoding (often Latin-1), which silently maps non-Latin
    /// codepoints to '?' before the consumer has a chance to choose a different
    /// fallback. UniProt taxonomy names contain Greek letters (β, …) and accented
    /// Latin letters; we need both visible to the PEFF emitter's NFKD-style helper.
    std::string xmlchToString(const XMLCh* ch)
    {
      if (ch == nullptr) return std::string();
      std::string out;
      while (*ch != 0)
      {
        uint32_t cp = static_cast<uint32_t>(*ch);
        ++ch;
        // Decode UTF-16 surrogate pair if present.
        if (cp >= 0xD800 && cp <= 0xDBFF && *ch >= 0xDC00 && *ch <= 0xDFFF)
        {
          uint32_t low = static_cast<uint32_t>(*ch);
          ++ch;
          cp = 0x10000 + ((cp - 0xD800) << 10) + (low - 0xDC00);
        }
        appendUtf8(out, cp);
      }
      return out;
    }

    /// Append a Xerces UTF-16 character buffer (with explicit length) as UTF-8 to @p out.
    void appendXmlchUtf8(std::string& out, const XMLCh* ch, std::size_t length)
    {
      for (std::size_t i = 0; i < length; ++i)
      {
        uint32_t cp = static_cast<uint32_t>(ch[i]);
        if (cp >= 0xD800 && cp <= 0xDBFF && i + 1 < length && ch[i + 1] >= 0xDC00 && ch[i + 1] <= 0xDFFF)
        {
          uint32_t low = static_cast<uint32_t>(ch[++i]);
          cp = 0x10000 + ((cp - 0xD800) << 10) + (low - 0xDC00);
        }
        appendUtf8(out, cp);
      }
    }
  } // namespace

  UniProtXMLHandler::UniProtXMLHandler(const std::string& filename, EntryCallback callback)
    : XMLHandler(filename, /*version=*/""),
      callback_(std::move(callback))
  {
  }

  UniProtXMLHandler::~UniProtXMLHandler() = default;

  void UniProtXMLHandler::resetEntry_()
  {
    current_entry_ = UniProtEntry{};
    full_name_captured_ = false;
    recommended_name_depth_ = 0;
    gene_depth_ = 0;
    organism_depth_ = 0;
    sequence_depth_ = 0;
    feature_depth_ = 0;
  }

  void UniProtXMLHandler::resetFeature_()
  {
    current_feature_ = UniProtFeature{};
  }

  int UniProtXMLHandler::parsePosition_(const std::string& attr)
  {
    if (attr.empty()) return 0;
    char* end_ptr = nullptr;
    long value = std::strtol(attr.c_str(), &end_ptr, 10);
    if (end_ptr == attr.c_str() || *end_ptr != '\0' || value < 0)
    {
      return 0;  // status="unknown" -> attribute absent / non-numeric -> 0
    }
    return static_cast<int>(value);
  }

  void UniProtXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/,
                                       const XMLCh* const qname, const xercesc::Attributes& attrs)
  {
    ++depth_;
    const std::string tag = xmlchToString(qname);

    // Inside the alternative-products comment? swallow the whole subtree.
    if (alt_products_depth_ != 0) return;

    // Reset per-element capture at every start so the next end-tag only sees its own characters.
    char_buf_.clear();
    capture_ = CaptureTarget::None;

    // <entry dataset="..." version="...">
    if (tag == "entry" && entry_depth_ == 0)
    {
      resetEntry_();
      entry_depth_ = depth_;
      current_entry_.dataset = attrValue(attrs, "dataset");
      current_entry_.entry_version = attrValue(attrs, "version");
      return;
    }

    if (entry_depth_ == 0) return;  // outside any <entry>

    // <comment type="alternative products"> -> skip whole subtree
    if (tag == "comment" && attrValue(attrs, "type") == "alternative products")
    {
      alt_products_depth_ = depth_;
      return;
    }

    // <feature type="..." description="...">
    if (tag == "feature" && feature_depth_ == 0)
    {
      resetFeature_();
      feature_depth_ = depth_;
      current_feature_.type = attrValue(attrs, "type");
      current_feature_.description = attrValue(attrs, "description");
      return;
    }

    // Inside a feature: <position>, <begin>, <end>, <original>, <variation>.
    if (feature_depth_ != 0)
    {
      if (tag == "position")
      {
        current_feature_.has_position = true;
        current_feature_.position = parsePosition_(attrValue(attrs, "position"));
      }
      else if (tag == "begin")
      {
        current_feature_.has_range = true;
        current_feature_.begin = parsePosition_(attrValue(attrs, "position"));
      }
      else if (tag == "end")
      {
        current_feature_.has_range = true;
        current_feature_.end = parsePosition_(attrValue(attrs, "position"));
      }
      else if (tag == "original")
      {
        capture_ = CaptureTarget::FeatureOriginal;
      }
      else if (tag == "variation" && current_feature_.variation.empty())
      {
        // C# captures only the first <variation>.
        capture_ = CaptureTarget::FeatureVariation;
      }
      return;
    }

    // <accession> (any) — first sets primary, rest go to alt_accessions.
    if (tag == "accession")
    {
      capture_ = CaptureTarget::Accession;
      return;
    }

    // <protein><recommendedName>: gate so we capture only the recommended <fullName>.
    if (tag == "recommendedName" && recommended_name_depth_ == 0)
    {
      recommended_name_depth_ = depth_;
      return;
    }
    if (tag == "fullName" && recommended_name_depth_ != 0 && !full_name_captured_)
    {
      capture_ = CaptureTarget::FullName;
      return;
    }

    // <gene>: scoped scan for <name type="primary">.
    if (tag == "gene" && gene_depth_ == 0)
    {
      gene_depth_ = depth_;
      return;
    }
    if (gene_depth_ != 0 && tag == "name" && attrValue(attrs, "type") == "primary"
        && current_entry_.primary_gene.empty())
    {
      capture_ = CaptureTarget::PrimaryGene;
      return;
    }

    // <organism>: scoped scan for <name type="scientific"> + NCBI Taxonomy <dbReference>.
    if (tag == "organism" && organism_depth_ == 0)
    {
      organism_depth_ = depth_;
      return;
    }
    if (organism_depth_ != 0 && tag == "name" && attrValue(attrs, "type") == "scientific"
        && current_entry_.tax_name.empty())
    {
      capture_ = CaptureTarget::TaxName;
      return;
    }
    if (organism_depth_ != 0 && tag == "dbReference" && attrValue(attrs, "type") == "NCBI Taxonomy"
        && current_entry_.ncbi_tax_id.empty())
    {
      current_entry_.ncbi_tax_id = attrValue(attrs, "id");
      return;
    }

    // <proteinExistence type="...">
    if (tag == "proteinExistence" && current_entry_.protein_existence.empty())
    {
      current_entry_.protein_existence = attrValue(attrs, "type");
      return;
    }

    // <sequence length="..." version="...">: only the canonical sequence carries a "length" attr;
    // isoform <sequence> elements (inside <comment type="alternative products"> or VAR_SEQ) lack it.
    if (tag == "sequence" && depth_ == entry_depth_ + 1 && !attrValue(attrs, "length").empty())
    {
      sequence_depth_ = depth_;
      if (current_entry_.sequence_version.empty())
      {
        current_entry_.sequence_version = attrValue(attrs, "version");
      }
      capture_ = CaptureTarget::Sequence;
      return;
    }

    // Top-level <name> (entry mnemonic): direct child of <entry>, outside protein/gene/organism.
    if (tag == "name" && depth_ == entry_depth_ + 1 && gene_depth_ == 0 && organism_depth_ == 0
        && recommended_name_depth_ == 0 && current_entry_.name.empty())
    {
      capture_ = CaptureTarget::EntryName;
      return;
    }
  }

  void UniProtXMLHandler::characters(const XMLCh* const chars, const XMLSize_t length)
  {
    if (capture_ == CaptureTarget::None) return;
    // Use the length-aware UTF-16 → UTF-8 helper so we don't depend on a
    // NUL-terminated XMLCh* (Xerces guarantees the buffer length, not termination).
    appendXmlchUtf8(char_buf_, chars, length);
  }

  void UniProtXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/,
                                     const XMLCh* const qname)
  {
    const std::string tag = xmlchToString(qname);

    // Close the alt-products skip gate at its matching </comment>.
    if (alt_products_depth_ != 0)
    {
      if (depth_ == alt_products_depth_ && tag == "comment") alt_products_depth_ = 0;
      --depth_;
      return;
    }

    // Deposit buffered character data, if any, at the right destination.
    if (capture_ != CaptureTarget::None)
    {
      switch (capture_)
      {
        case CaptureTarget::Accession:
          if (current_entry_.accession.empty()) current_entry_.accession = char_buf_;
          else current_entry_.alt_accessions.push_back(char_buf_);
          break;
        case CaptureTarget::EntryName:
          current_entry_.name = char_buf_;
          break;
        case CaptureTarget::FullName:
          current_entry_.full_name = char_buf_;
          full_name_captured_ = true;
          break;
        case CaptureTarget::PrimaryGene:
          current_entry_.primary_gene = char_buf_;
          break;
        case CaptureTarget::TaxName:
          current_entry_.tax_name = char_buf_;
          break;
        case CaptureTarget::Sequence:
          current_entry_.sequence = stripWhitespace(char_buf_);
          break;
        case CaptureTarget::FeatureOriginal:
          current_feature_.original = char_buf_;
          break;
        case CaptureTarget::FeatureVariation:
          current_feature_.variation = char_buf_;
          break;
        case CaptureTarget::None:
          break;
      }
      capture_ = CaptureTarget::None;
      char_buf_.clear();
    }

    // Close subtree gates when their matching end tag arrives.
    if (sequence_depth_ != 0 && depth_ == sequence_depth_ && tag == "sequence")
    {
      sequence_depth_ = 0;
    }
    if (recommended_name_depth_ != 0 && depth_ == recommended_name_depth_ && tag == "recommendedName")
    {
      recommended_name_depth_ = 0;
    }
    if (gene_depth_ != 0 && depth_ == gene_depth_ && tag == "gene")
    {
      gene_depth_ = 0;
    }
    if (organism_depth_ != 0 && depth_ == organism_depth_ && tag == "organism")
    {
      organism_depth_ = 0;
    }
    if (feature_depth_ != 0 && depth_ == feature_depth_ && tag == "feature")
    {
      current_entry_.features.push_back(std::move(current_feature_));
      resetFeature_();
      feature_depth_ = 0;
    }
    if (entry_depth_ != 0 && depth_ == entry_depth_ && tag == "entry")
    {
      if (callback_) callback_(std::move(current_entry_));
      resetEntry_();
      entry_depth_ = 0;
    }

    --depth_;
  }

} // namespace OpenMS::Internal
