// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/UniProtXMLFile.h>

#include <functional>
#include <string>

namespace OpenMS::Internal
{
  /**
    @brief SAX handler for UniProtKB XML &lt;entry&gt; documents.

    Streams individual &lt;entry&gt; elements into UniProtEntry POD instances and
    delivers each via a callback at &lt;/entry&gt;. Implements the same scoping
    rules as the upstream C# UniPEFF tool: isoform sequences (missing the
    "length" attribute on &lt;sequence&gt;) and the &lt;comment type="alternative
    products"&gt; subtree are skipped, &lt;gene&gt; and &lt;organism&gt; are scanned via
    bounded depth tracking so their inner &lt;name&gt; elements never leak into the
    top-level "entry mnemonic" capture, and &lt;feature&gt; contents are buffered
    and classified at the closing tag rather than incrementally.
  */
  class OPENMS_DLLAPI UniProtXMLHandler :
    public XMLHandler
  {
  public:
    /// Callback invoked once per &lt;/entry&gt;: receives ownership of the populated
    /// UniProtEntry. The handler resets its working state immediately after, so the
    /// callback is the only place the entry is reachable.
    using EntryCallback = std::function<void(UniProtEntry&&)>;

    /// Where the next character-data run, if buffered, should be deposited
    /// at the closing tag of its enclosing element.
    enum class CaptureTarget
    {
      None,
      Accession,
      EntryName,
      FullName,
      PrimaryGene,
      TaxName,
      Sequence,
      FeatureOriginal,
      FeatureVariation
    };

    /// Build a handler that delivers each parsed UniProtEntry to @p callback.
    UniProtXMLHandler(const std::string& filename, EntryCallback callback);

    /// Destructor.
    ~UniProtXMLHandler() override;

    void startElement(const XMLCh* const uri, const XMLCh* const local_name,
                      const XMLCh* const qname, const xercesc::Attributes& attrs) override;
    void endElement(const XMLCh* const uri, const XMLCh* const local_name,
                    const XMLCh* const qname) override;
    void characters(const XMLCh* const chars, const XMLSize_t length) override;

  private:
    EntryCallback callback_;
    UniProtEntry current_entry_;

    /// Depth of the most recent element start, used to scope subtree handling.
    int depth_{0};

    /// Subtree gates. Each holds the depth at which the corresponding element opened
    /// so we know when its closing tag fires (depth_ == gate value during endElement).
    /// 0 means "not currently inside this subtree".
    int entry_depth_{0};                ///< &lt;entry&gt;
    int recommended_name_depth_{0};     ///< &lt;protein&gt;/&lt;recommendedName&gt;
    int gene_depth_{0};                 ///< &lt;gene&gt;
    int organism_depth_{0};             ///< &lt;organism&gt;
    int alt_products_depth_{0};         ///< &lt;comment type="alternative products"&gt; (whole subtree skipped)
    int feature_depth_{0};              ///< &lt;feature&gt;
    int sequence_depth_{0};             ///< &lt;sequence length="..."&gt; (only the canonical sequence; isoform &lt;sequence&gt; ignored)

    /// Captures the &lt;fullName&gt; the first time we encounter it inside
    /// &lt;protein&gt;/&lt;recommendedName&gt;; later &lt;fullName&gt; elements (e.g. alternative names) are skipped.
    bool full_name_captured_{false};

    /// Buffer for character data; appended-to in characters(), consumed in endElement().
    std::string char_buf_;

    /// Destination for the next buffered run of character data (None == not capturing).
    CaptureTarget capture_{CaptureTarget::None};

    /// Per-feature working state (used between startElement("feature") and endElement("feature")).
    UniProtFeature current_feature_;

    /// Whether the next &lt;name&gt; encountered inside the current &lt;gene&gt; subtree is a primary name.
    bool gene_name_is_primary_{false};
    /// Whether the next &lt;name&gt; encountered inside the current &lt;organism&gt; subtree is the scientific name.
    bool organism_name_is_scientific_{false};

    /// Clear all per-entry state so the next &lt;entry&gt; starts fresh.
    void resetEntry_();
    /// Clear all per-feature state so the next &lt;feature&gt; starts fresh.
    void resetFeature_();

    /// Parse a UniProt position attribute string; returns 0 when the attribute is
    /// absent or non-numeric (UniProt encodes "unknown" by omitting the attribute).
    static int parsePosition_(const std::string& attr);
  };

} // namespace OpenMS::Internal
