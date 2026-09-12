// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>

#include <functional>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief A single &lt;feature&gt; element from a UniProtKB XML entry.

    Fields mirror the UniProt schema after the &lt;location&gt; has been collapsed:
    a feature is either a single point (HasPosition / position) or a range
    (HasRange / begin / end). A status="unknown" coordinate is encoded as 0
    (the caller decides whether that is legal for its target PEFF key).
  */
  struct OPENMS_DLLAPI UniProtFeature
  {
    std::string type;          ///< feature/\@type, e.g. "modified residue", "disulfide bond", "sequence variant"
    std::string id;            ///< feature/\@id, e.g. "VSP_021275" (referenced by isoform definitions)
    std::string description;   ///< feature/\@description (raw, before any cleanup)
    bool has_position{false};  ///< a single &lt;position&gt; element was present
    bool has_range{false};     ///< a &lt;begin&gt; / &lt;end&gt; pair was present
    int  position{0};          ///< 1-based position (0 = unknown / absent)
    int  begin{0};             ///< 1-based range start (0 = unknown / absent)
    int  end{0};               ///< 1-based range end   (0 = unknown / absent)
    std::string original;      ///< &lt;original&gt; text (sequence variant)
    std::string variation;     ///< &lt;variation&gt; text (sequence variant; first occurrence only)
  };

  /**
    @brief One &lt;isoform&gt; definition from &lt;comment type="alternative products"&gt;.

    UniProtKB XML does not store isoform sequences verbatim: an isoform either
    "displays" the canonical sequence (@c sequence_type "displayed") or is
    "described" by the splice-variant features listed in @c sequence_ref, which
    transform the canonical sequence (a referenced range without a replacement
    is deleted, one with a &lt;variation&gt; is substituted). The remaining types
    ("external", "not described") carry no reconstructable sequence.
  */
  struct OPENMS_DLLAPI UniProtIsoform
  {
    std::string id;             ///< first &lt;id&gt; (isoform accession, e.g. "P02768-2"); later ids are secondary
    std::string name;           ///< first &lt;name&gt; (e.g. "2" or "VEGF121"); later &lt;name&gt; elements are synonyms
    std::string sequence_type;  ///< &lt;sequence type="..."&gt;: "displayed", "described", "external" or "not described"
    std::string sequence_ref;   ///< &lt;sequence ref="..."&gt;: space-separated &lt;feature type="splice variant"&gt; ids
  };

  /**
    @brief A single &lt;entry&gt; from a UniProtKB XML file, in a parser-neutral form.

    Captures the fields needed to build PEFF descriptor lines. The canonical
    sequence is stored verbatim; isoform sequences are not present in the XML
    and must be reconstructed by applying the splice-variant features listed in
    each UniProtIsoform's @c sequence_ref to the canonical sequence.
  */
  struct OPENMS_DLLAPI UniProtEntry
  {
    std::string dataset;                       ///< entry/\@dataset, e.g. "Swiss-Prot" or "TrEMBL"
    std::string accession;                     ///< first &lt;accession&gt; (primary id)
    std::vector<std::string> alt_accessions;   ///< second and subsequent &lt;accession&gt; entries
    std::string name;                          ///< first &lt;name&gt; under &lt;entry&gt; (mnemonic id, e.g. "KSINK_HUMAN")
    std::string full_name;                     ///< &lt;protein&gt;/&lt;recommendedName&gt;/&lt;fullName&gt; (first occurrence)
    std::string primary_gene;                  ///< &lt;gene&gt;/&lt;name type="primary"&gt; (first occurrence)
    std::string ncbi_tax_id;                   ///< &lt;organism&gt;/&lt;dbReference type="NCBI Taxonomy" id="..."&gt;
    std::string tax_name;                      ///< &lt;organism&gt;/&lt;name type="scientific"&gt;
    std::string protein_existence;             ///< &lt;proteinExistence type="..."&gt; (raw type string, e.g. "evidence at protein level")
    std::string sequence_version;              ///< &lt;sequence version="..."&gt;
    std::string entry_version;                 ///< &lt;entry version="..."&gt;
    std::string sequence;                      ///< canonical &lt;sequence&gt; text (whitespace stripped); isoform sequences are not stored in the XML
    std::vector<UniProtFeature> features;      ///< &lt;feature&gt; elements in document order
    std::vector<UniProtIsoform> isoforms;      ///< isoform definitions from &lt;comment type="alternative products"&gt;
  };

  /**
    @brief Reads UniProtKB XML protein databases (.xml or transparently .xml.gz).

    Streams the document via Xerces SAX. The on-disk format is the UniProtKB
    XML release schema. The &lt;comment type="alternative products"&gt; subtree is
    parsed into UniProtIsoform records (isoform sequences themselves are not
    stored in the XML — see UniProtIsoform on how to reconstruct them). Gzip is
    detected by magic-byte sniffing in the inherited XMLFile parser, so no
    special handling is required for .xml.gz inputs.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI UniProtXMLFile :
    public Internal::XMLFile
  {
  public:
    /// Default constructor.
    UniProtXMLFile();

    /// Destructor.
    ~UniProtXMLFile() override;

    /**
      @brief Load all entries from a UniProtKB XML file into memory.

      @param[in]  filename Path to a .xml or .xml.gz UniProtKB file (gzip auto-detected).
      @param[out] entries  Filled with one UniProtEntry per &lt;entry&gt; element.

      @throw Exception::FileNotFound if @p filename does not exist.
      @throw Exception::ParseError   if the file is not well-formed XML.
    */
    void load(const std::string& filename, std::vector<UniProtEntry>& entries);

    /**
      @brief Stream entries one at a time via a user-supplied callback.

      The callback receives ownership of each fully populated UniProtEntry as
      soon as its &lt;/entry&gt; is reached. This is the preferred path for
      databases too large to materialise in full (e.g. TrEMBL).

      @param[in] filename Path to a .xml or .xml.gz UniProtKB file (gzip auto-detected).
      @param[in] callback Invoked once per parsed entry; the entry is moved-from after the call.

      @throw Exception::FileNotFound if @p filename does not exist.
      @throw Exception::ParseError   if the file is not well-formed XML.
    */
    void loadStreaming(const std::string& filename, const std::function<void(UniProtEntry&&)>& callback);

  private:
    UniProtXMLFile(const UniProtXMLFile&) = delete;
    UniProtXMLFile& operator=(const UniProtXMLFile&) = delete;
  };

} // namespace OpenMS
