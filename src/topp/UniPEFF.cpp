// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/UniProtXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_UniPEFF UniPEFF

@brief Converts a UniProtKB XML protein database to PEFF 1.0 with per-entry
modification, processing, variant and disulfide-bond annotations.

This is an OpenMS-native port of the standalone C# UniPEFF tool
(David L. Tabb, UMC Groningen). For each entry it emits a PEFF descriptor
line with \\PName \\GName \\NcbiTaxId \\TaxName \\Length \\SV \\EV \\PE \\ID
\\AltAC \\ModResPsi \\ModResUnimod \\ModRes \\VariantSimple \\VariantComplex
\\Processed and \\DisulfideBond (the latter unless
@c -omit_amino_acid_modifications is set).

UniProt's <code>\<feature type="disulfide bond"\></code> entries are translated
into half-cystine modifications (PSI-MOD:00798) which are then merged into
the main modified-residue list in position-sorted order. Disulfide
connectivity is reported as <code>\\DisulfideBond=(bond:idA,idB)</code>
tuples referencing <code>id:</code> prefixes on the half-cystine tuples
inside \\ModResPsi. By default only half-cystines that participate in a
documented intrachain bond are labeled, in bond order: the k-th bond labels
its two half-cystines 2k-1 and 2k and is itself labeled k (1-based; see
issue 9829). Interchain half-cystines (single-position features) stay
unlabeled. With @c -annotation_identifiers (PEFF "Option B") every
annotation tuple instead carries a global sequential 1-based id and
\\DisulfideBond references those ids.

Modification accession lookup uses UniProt's <em>ptmlist.txt</em> (a snapshot
is bundled under @c share/OpenMS/CHEMISTRY/UniProt_ptmlist.txt); override
with @c -ptmlist. Canonical OBO names come from <em>PSI-MOD.obo</em> (bundled)
and an optional <em>unimod.obo</em>; without them, names fall back to the
UniProt @c ptmlist.txt ID and a warning is printed.

Both plain @c .xml and @c .xml.gz UniProt inputs are accepted (gzip is
auto-detected by the underlying parser).

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_UniPEFF.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_UniPEFF.html
*/

/// @cond TOPPCLASSES

namespace
{
  // ──────────────────────────────────────────────────────────────────
  // ptmlist.txt — UniProt PTM controlled vocabulary
  // ──────────────────────────────────────────────────────────────────

  struct PtmEntry
  {
    std::string id;               ///< UniProt PTM name (the lookup key from feature/@description)
    std::string accession;        ///< UniProt internal accession, e.g. "PTM-0253"
    std::string psi_mod;          ///< "MOD:xxxxx" or empty
    std::string unimod;           ///< "N" (numeric Unimod id) or empty (NOT prefixed with "UNIMOD:")
  };

  using PtmMap = std::unordered_map<std::string, PtmEntry>;

  /// Parse UniProt's @c ptmlist.txt block format:
  /// @verbatim
  ///   ID   Phosphoserine
  ///   AC   PTM-0253
  ///   DR   PSI-MOD; MOD:00046.
  ///   DR   Unimod; 21.
  ///   //
  /// @endverbatim
  /// "Half cystine" is hard-coded (it is absent from the upstream ptmlist.txt
  /// even though UniProt's @c <feature type="disulfide bond"> needs it).
  PtmMap loadPtmlist(const std::string& filename)
  {
    PtmMap result;
    // Hard-coded entry that mirrors the C# tool's setup.
    result.emplace("Half cystine", PtmEntry{"Half cystine", "", "MOD:00798", "374"});

    std::ifstream in(filename);
    if (!in)
    {
      OPENMS_LOG_WARN << "UniPEFF: cannot read ptmlist '" << filename
                      << "' — only the hard-coded 'Half cystine' entry is available." << std::endl;
      return result;
    }

    PtmEntry current;
    bool have_current = false;
    std::string line;
    auto commit = [&]() {
      if (have_current && !current.id.empty()) result[current.id] = current;
      current = PtmEntry{};
      have_current = false;
    };

    while (std::getline(in, line))
    {
      // std::getline strips the trailing LF but leaves a CR on CRLF files (Windows-edited
      // ptmlist.txt is a real concern). Trim any trailing whitespace so the fixed-offset
      // accession extraction below — `rest.size() - 10`/`-9` — doesn't off-by-one.
      while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t'))
      {
        line.pop_back();
      }
      if (line.size() < 2) continue;
      const std::string key = line.substr(0, 2);
      if (key == "ID")
      {
        commit();
        if (line.size() > 5) { current.id = line.substr(5); have_current = true; }
      }
      else if (key == "AC")
      {
        if (line.size() > 5) current.accession = line.substr(5);
      }
      else if (key == "DR")
      {
        if (line.size() <= 5) continue;
        std::string rest = line.substr(5);  // e.g. "PSI-MOD; MOD:00046." or "Unimod; 21."
        if (rest.rfind("PSI-MOD;", 0) == 0 && rest.size() > 10)
        {
          // skip "PSI-MOD; " (9 chars), drop trailing "."
          current.psi_mod = rest.substr(9, rest.size() - 10);
        }
        else if (rest.rfind("Unimod;", 0) == 0 && rest.size() > 9)
        {
          // skip "Unimod; " (8 chars), drop trailing "."
          current.unimod = rest.substr(8, rest.size() - 9);
        }
      }
      else if (key == "//")
      {
        commit();
      }
    }
    commit();  // in case the file does not end with "//"
    return result;
  }

  // ──────────────────────────────────────────────────────────────────
  // OBO name map: accession -> canonical "name:" field
  // ──────────────────────────────────────────────────────────────────

  using OboNameMap = std::unordered_map<std::string, std::string>;

  /// Parse a standard OBO file and build accession→name. Returns an empty map
  /// if the file is missing or unreadable (mirrors the C# fallback behaviour).
  /// Only [Term] stanzas are ingested. Obsolete terms are kept on purpose.
  OboNameMap loadOboNames(const std::string& filename)
  {
    OboNameMap result;
    if (filename.empty()) return result;
    std::ifstream in(filename);
    if (!in)
    {
      OPENMS_LOG_WARN << "UniPEFF: cannot read OBO file '" << filename
                      << "' — modification names will fall back to UniProt ptmlist IDs." << std::endl;
      return result;
    }
    std::string line;
    std::string current_id;
    bool in_term = false;
    while (std::getline(in, line))
    {
      if (!line.empty() && line[0] == '[')
      {
        in_term = (line.rfind("[Term]", 0) == 0);
        current_id.clear();
      }
      else if (in_term && line.rfind("id: ", 0) == 0)
      {
        current_id = line.substr(4);
        // trim
        while (!current_id.empty() && std::isspace(static_cast<unsigned char>(current_id.back()))) current_id.pop_back();
      }
      else if (in_term && line.rfind("name: ", 0) == 0 && !current_id.empty())
      {
        std::string name = line.substr(6);
        while (!name.empty() && std::isspace(static_cast<unsigned char>(name.back()))) name.pop_back();
        if (result.find(current_id) == result.end()) result[current_id] = name;
        current_id.clear();
      }
    }
    return result;
  }

  // ──────────────────────────────────────────────────────────────────
  // Small format helpers (mirror C# PEFFentry / PEFFmodel)
  // ──────────────────────────────────────────────────────────────────

  /// "TrEMBL" → "tr"; everything else (Swiss-Prot/null/unknown) → "sp".
  std::string prefixForDataset(const std::string& dataset)
  {
    return (dataset == "TrEMBL") ? "tr" : "sp";
  }

  /// UniProt's <proteinExistence type="..."> → PEFF \PE digit ("1".."5") or "".
  std::string proteinExistenceCode(const std::string& type)
  {
    if (type == "evidence at protein level")    return "1";
    if (type == "evidence at transcript level") return "2";
    if (type == "inferred from homology")       return "3";
    if (type == "predicted")                    return "4";
    if (type == "uncertain")                    return "5";
    return std::string();
  }

  /// PEFF VariantSimple newAA must be an amino-acid letter (or '*').
  bool isResidueCode(char c)
  {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '*';
  }

  /// 0 → "?", else decimal.
  std::string positionOrUnknown(int p)
  {
    return p == 0 ? std::string("?") : std::to_string(p);
  }

  /// Strip the "(Microbial infection) " prefix and any trailing "; ..." qualifier.
  /// Matches the C# CleanPtmDescription verbatim.
  std::string cleanPtmDescription(std::string description)
  {
    static const std::string kMicrobialInfection = "(Microbial infection) ";
    if (description.rfind(kMicrobialInfection, 0) == 0)
    {
      description = description.substr(kMicrobialInfection.size());
    }
    const auto semi = description.find(';');
    if (semi != std::string::npos)
    {
      description = description.substr(0, semi);
    }
    return description;
  }

  // ──────────────────────────────────────────────────────────────────
  // PEFF text escaping (UniPEFF semantics)
  // ──────────────────────────────────────────────────────────────────

  /// PEFF allows only ASCII. NFKD-normalise so that combining diacritics
  /// (é, ö, …) collapse to their base letter; replace any remaining non-ASCII
  /// byte with '?'. Without ICU we cover the Latin-1 supplement (the only
  /// non-ASCII characters that occur with any frequency in UniProtKB taxonomy
  /// and protein names). Greek/Cyrillic/etc. fall back to '?', matching the
  /// C# tool's behaviour for characters without an ASCII decomposition.
  std::string toAscii(const std::string& value)
  {
    bool already_ascii = true;
    for (unsigned char c : value) { if (c > 0x7F) { already_ascii = false; break; } }
    if (already_ascii) return value;

    // Decode UTF-8 to code points, map known Latin-1 characters to ASCII,
    // emit '?' for everything else. This is a minimal transliterator that
    // covers UniProtKB's common cases (NFKD of é/ö/ß/… via lookup); for the
    // few remaining cases the C# NFKD would also fall back to '?'.
    static const std::unordered_map<uint32_t, const char*> kMap = {
      // ── Latin-1 Supplement (U+00C0..U+00FF) ──────────────────────────────
      {0x00C0, "A"}, {0x00C1, "A"}, {0x00C2, "A"}, {0x00C3, "A"}, {0x00C4, "A"}, {0x00C5, "A"},
      {0x00C6, "AE"}, {0x00C7, "C"},
      {0x00C8, "E"}, {0x00C9, "E"}, {0x00CA, "E"}, {0x00CB, "E"},
      {0x00CC, "I"}, {0x00CD, "I"}, {0x00CE, "I"}, {0x00CF, "I"},
      {0x00D0, "D"}, {0x00D1, "N"},
      {0x00D2, "O"}, {0x00D3, "O"}, {0x00D4, "O"}, {0x00D5, "O"}, {0x00D6, "O"}, {0x00D8, "O"},
      {0x00D9, "U"}, {0x00DA, "U"}, {0x00DB, "U"}, {0x00DC, "U"},
      {0x00DD, "Y"}, {0x00DE, "Th"}, {0x00DF, "ss"},
      {0x00E0, "a"}, {0x00E1, "a"}, {0x00E2, "a"}, {0x00E3, "a"}, {0x00E4, "a"}, {0x00E5, "a"},
      {0x00E6, "ae"}, {0x00E7, "c"},
      {0x00E8, "e"}, {0x00E9, "e"}, {0x00EA, "e"}, {0x00EB, "e"},
      {0x00EC, "i"}, {0x00ED, "i"}, {0x00EE, "i"}, {0x00EF, "i"},
      {0x00F0, "d"}, {0x00F1, "n"},
      {0x00F2, "o"}, {0x00F3, "o"}, {0x00F4, "o"}, {0x00F5, "o"}, {0x00F6, "o"}, {0x00F8, "o"},
      {0x00F9, "u"}, {0x00FA, "u"}, {0x00FB, "u"}, {0x00FC, "u"},
      {0x00FD, "y"}, {0x00FE, "th"}, {0x00FF, "y"},
      // ── Latin Extended-A (U+0100..U+017F) — full block ──────────────────
      // NFKD of every character in this block decomposes to ASCII (plus the
      // few digraphs IJ/ij/OE/oe).
      {0x0100, "A"}, {0x0101, "a"}, {0x0102, "A"}, {0x0103, "a"}, {0x0104, "A"}, {0x0105, "a"},
      {0x0106, "C"}, {0x0107, "c"}, {0x0108, "C"}, {0x0109, "c"}, {0x010A, "C"}, {0x010B, "c"},
      {0x010C, "C"}, {0x010D, "c"}, {0x010E, "D"}, {0x010F, "d"}, {0x0110, "D"}, {0x0111, "d"},
      {0x0112, "E"}, {0x0113, "e"}, {0x0114, "E"}, {0x0115, "e"}, {0x0116, "E"}, {0x0117, "e"},
      {0x0118, "E"}, {0x0119, "e"}, {0x011A, "E"}, {0x011B, "e"}, {0x011C, "G"}, {0x011D, "g"},
      {0x011E, "G"}, {0x011F, "g"}, {0x0120, "G"}, {0x0121, "g"}, {0x0122, "G"}, {0x0123, "g"},
      {0x0124, "H"}, {0x0125, "h"}, {0x0126, "H"}, {0x0127, "h"},
      {0x0128, "I"}, {0x0129, "i"}, {0x012A, "I"}, {0x012B, "i"}, {0x012C, "I"}, {0x012D, "i"},
      {0x012E, "I"}, {0x012F, "i"}, {0x0130, "I"}, {0x0131, "i"},
      {0x0132, "IJ"}, {0x0133, "ij"}, {0x0134, "J"}, {0x0135, "j"},
      {0x0136, "K"}, {0x0137, "k"}, {0x0138, "k"},
      {0x0139, "L"}, {0x013A, "l"}, {0x013B, "L"}, {0x013C, "l"}, {0x013D, "L"}, {0x013E, "l"},
      {0x013F, "L"}, {0x0140, "l"}, {0x0141, "L"}, {0x0142, "l"},
      {0x0143, "N"}, {0x0144, "n"}, {0x0145, "N"}, {0x0146, "n"}, {0x0147, "N"}, {0x0148, "n"},
      {0x0149, "n"}, {0x014A, "N"}, {0x014B, "n"},
      {0x014C, "O"}, {0x014D, "o"}, {0x014E, "O"}, {0x014F, "o"}, {0x0150, "O"}, {0x0151, "o"},
      {0x0152, "OE"}, {0x0153, "oe"},
      {0x0154, "R"}, {0x0155, "r"}, {0x0156, "R"}, {0x0157, "r"}, {0x0158, "R"}, {0x0159, "r"},
      {0x015A, "S"}, {0x015B, "s"}, {0x015C, "S"}, {0x015D, "s"}, {0x015E, "S"}, {0x015F, "s"},
      {0x0160, "S"}, {0x0161, "s"},
      {0x0162, "T"}, {0x0163, "t"}, {0x0164, "T"}, {0x0165, "t"}, {0x0166, "T"}, {0x0167, "t"},
      {0x0168, "U"}, {0x0169, "u"}, {0x016A, "U"}, {0x016B, "u"}, {0x016C, "U"}, {0x016D, "u"},
      {0x016E, "U"}, {0x016F, "u"}, {0x0170, "U"}, {0x0171, "u"}, {0x0172, "U"}, {0x0173, "u"},
      {0x0174, "W"}, {0x0175, "w"}, {0x0176, "Y"}, {0x0177, "y"}, {0x0178, "Y"},
      {0x0179, "Z"}, {0x017A, "z"}, {0x017B, "Z"}, {0x017C, "z"}, {0x017D, "Z"}, {0x017E, "z"},
      {0x017F, "s"},
    };

    std::string out;
    out.reserve(value.size());
    size_t i = 0;
    while (i < value.size())
    {
      unsigned char b0 = static_cast<unsigned char>(value[i]);
      if (b0 < 0x80) { out.push_back(static_cast<char>(b0)); ++i; continue; }
      // Decode one UTF-8 code point
      uint32_t cp = 0;
      int len = 0;
      if ((b0 & 0xE0) == 0xC0) { cp = b0 & 0x1F; len = 2; }
      else if ((b0 & 0xF0) == 0xE0) { cp = b0 & 0x0F; len = 3; }
      else if ((b0 & 0xF8) == 0xF0) { cp = b0 & 0x07; len = 4; }
      else { out.push_back('?'); ++i; continue; }  // invalid start byte
      if (i + len > value.size()) { out.push_back('?'); ++i; continue; }
      bool ok = true;
      for (int k = 1; k < len; ++k)
      {
        unsigned char bk = static_cast<unsigned char>(value[i + k]);
        if ((bk & 0xC0) != 0x80) { ok = false; break; }
        cp = (cp << 6) | (bk & 0x3F);
      }
      if (!ok) { out.push_back('?'); ++i; continue; }
      i += len;
      auto it = kMap.find(cp);
      if (it != kMap.end()) out.append(it->second);
      else out.push_back('?');
    }
    return out;
  }

  /// PEFF reserved-character escaping. Per UniPEFF: \\ and | are always
  /// escaped; only UNPAIRED parens get a backslash; balanced pairs are
  /// left intact so that descriptions like "N-linked (GlcNAc...)" survive.
  std::string escapePeff(const std::string& raw)
  {
    if (raw.empty()) return raw;
    const std::string ascii = toAscii(raw);
    // Mark indices of unpaired '(' / ')'.
    std::vector<bool> unpaired(ascii.size(), false);
    std::vector<size_t> open_stack;
    for (size_t i = 0; i < ascii.size(); ++i)
    {
      if (ascii[i] == '(') open_stack.push_back(i);
      else if (ascii[i] == ')')
      {
        if (!open_stack.empty()) open_stack.pop_back();
        else unpaired[i] = true;
      }
    }
    for (size_t idx : open_stack) unpaired[idx] = true;

    std::string out;
    out.reserve(ascii.size() + 4);
    for (size_t i = 0; i < ascii.size(); ++i)
    {
      char c = ascii[i];
      if (c == '\\' || c == '|' || ((c == '(' || c == ')') && unpaired[i])) out.push_back('\\');
      out.push_back(c);
    }
    return out;
  }

  // ──────────────────────────────────────────────────────────────────
  // Intermediate annotation structures (built per UniProtEntry)
  // ──────────────────────────────────────────────────────────────────

  /// Sentinel for "no annotation identifier assigned": the tuple is emitted without an id: prefix.
  constexpr uint32_t kNoId = std::numeric_limits<uint32_t>::max();

  struct ModResItem
  {
    int  position{0};        ///< 0 = unknown ('?')
    bool is_halfcys{false};  ///< for stable tie-break on merge
    const PtmEntry* ptm{nullptr};  ///< nullptr = generic (no accession), description in name
    std::string name;        ///< display name (UniProt PTM ID or fallback description)
    uint32_t annotation_id{kNoId};  ///< 1-based; all tuples in Option B, paired half-cystines in default mode
  };

  struct VariantSimpleItem
  {
    int  position{0};
    char new_aa{'?'};
    uint32_t annotation_id{kNoId};  ///< assigned in Option B only
  };

  struct VariantComplexItem
  {
    int begin{0};
    int end{0};
    std::string new_seq;
    uint32_t annotation_id{kNoId};  ///< assigned in Option B only
  };

  struct ProcessedItem
  {
    int begin{0};
    int end{0};
    std::string cv;          ///< PEFF CV accession, e.g. "PEFF:0001021"
    std::string type_name;   ///< human name, e.g. "signal peptide"
    uint32_t annotation_id{kNoId};  ///< assigned in Option B only
  };

  struct DisulfidePairItem
  {
    /// Indices into EntryAnnotations::mods after merge — resolved to annotation_ids during emission.
    size_t idx_a{0};
    size_t idx_b{0};
    bool valid{false};       ///< true only when both endpoints come from a <begin>/<end> half-cystine pair
    uint32_t annotation_id{kNoId};  ///< id of the \DisulfideBond tuple itself (kNoId = not emitted)
  };

  struct EntryAnnotations
  {
    std::vector<ModResItem> regular_mods;
    std::vector<ModResItem> halfcys;
    std::vector<DisulfidePairItem> disulfides;   ///< indices reference @c halfcys before merge, @c mods after merge
    std::vector<VariantSimpleItem> simple_variants;
    std::vector<VariantComplexItem> complex_variants;
    std::vector<ProcessedItem> processed;

    /// Final merged + sorted list (built by classifyAndMerge).
    std::vector<ModResItem> mods;
  };

  // ──────────────────────────────────────────────────────────────────
  // Feature classifier — mirrors C# FromUniProtXML feature switch
  // ──────────────────────────────────────────────────────────────────

  /// Lookup PTM by description in @p ptms, returns nullptr on miss.
  const PtmEntry* findPtm(const PtmMap& ptms, const std::string& description)
  {
    auto it = ptms.find(description);
    return (it == ptms.end()) ? nullptr : &it->second;
  }

  /// Classify a single <feature> and append into the right annotation list(s).
  /// @p record_processing / @p record_aa_mods / @p record_variants gate the three
  /// categories of -Omit* flags. Returns nothing; updates @p annotations in place.
  void classifyAndAppend(const UniProtFeature& f, const PtmMap& ptms, EntryAnnotations& a,
                         bool record_processing, bool record_aa_mods, bool record_variants,
                         const std::string& accession)
  {
    // Single-residue site: prefer <position>, else <begin>; unknown -> 0.
    const int mod_position = f.has_position ? f.position : (f.has_range ? f.begin : 0);

    if (record_processing)
    {
      if (f.type == "chain")
      {
        ProcessedItem p{f.has_range ? f.begin : 0, f.has_range ? f.end : 0,
                        "PEFF:0001020", "mature protein"};
        a.processed.push_back(std::move(p));
        return;
      }
      if (f.type == "initiator methionine")
      {
        const int pos = f.has_position ? f.position : 0;
        a.processed.push_back(ProcessedItem{pos, pos, "PEFF:0001035", "initiator methionine"});
        return;
      }
      if (f.type == "propeptide")
      {
        a.processed.push_back(ProcessedItem{f.has_range ? f.begin : 0, f.has_range ? f.end : 0,
                                            "PEFF:0001034", "propeptide"});
        return;
      }
      if (f.type == "signal peptide")
      {
        a.processed.push_back(ProcessedItem{f.has_range ? f.begin : 0, f.has_range ? f.end : 0,
                                            "PEFF:0001021", "signal peptide"});
        return;
      }
      if (f.type == "transit peptide")
      {
        a.processed.push_back(ProcessedItem{f.has_range ? f.begin : 0, f.has_range ? f.end : 0,
                                            "PEFF:0001022", "transit peptide"});
        return;
      }
    }

    if (record_aa_mods)
    {
      if (f.type == "cross-link")
      {
        std::string desc = cleanPtmDescription(f.description);
        if (desc.empty()) desc = "cross-link";
        // Generic ModRes (no CV accession) per UniPEFF policy.
        if (f.has_range)
        {
          a.regular_mods.push_back(ModResItem{f.begin, false, nullptr, desc});
          a.regular_mods.push_back(ModResItem{f.end,   false, nullptr, desc});
        }
        else
        {
          a.regular_mods.push_back(ModResItem{mod_position, false, nullptr, desc});
        }
        return;
      }
      if (f.type == "disulfide bond")
      {
        const PtmEntry* hc = findPtm(ptms, "Half cystine");
        if (f.has_range)
        {
          a.halfcys.push_back(ModResItem{f.begin, true, hc, hc ? hc->id : "half cystine"});
          const size_t first_idx = a.halfcys.size() - 1;
          a.halfcys.push_back(ModResItem{f.end,   true, hc, hc ? hc->id : "half cystine"});
          const size_t second_idx = a.halfcys.size() - 1;
          a.disulfides.push_back(DisulfidePairItem{first_idx, second_idx, true});
        }
        else
        {
          a.halfcys.push_back(ModResItem{mod_position, true, hc, hc ? hc->id : "half cystine"});
        }
        return;
      }
      if (f.type == "glycosylation site")
      {
        std::string desc = cleanPtmDescription(f.description);
        if (desc.empty()) desc = "glycosylation site";
        a.regular_mods.push_back(ModResItem{mod_position, false, nullptr, desc});
        return;
      }
      if (f.type == "lipid moiety-binding region")
      {
        std::string desc = cleanPtmDescription(f.description);
        if (desc.empty()) desc = "lipid moiety-binding region";
        const PtmEntry* p = findPtm(ptms, desc);
        a.regular_mods.push_back(ModResItem{mod_position, false, p, p ? p->id : desc});
        return;
      }
      if (f.type == "modified residue")
      {
        std::string desc = cleanPtmDescription(f.description);
        const PtmEntry* p = findPtm(ptms, desc);
        if (p != nullptr)
        {
          a.regular_mods.push_back(ModResItem{mod_position, false, p, p->id});
        }
        else
        {
          // Mirror C# stdout messages so users can sync their ptmlist.txt.
          OPENMS_LOG_INFO << "UniPEFF: failed to find this modified residue in ptmlist.txt: " << desc << std::endl;
        }
        return;
      }
    }

    if (record_variants && f.type == "sequence variant")
    {
      const std::string new_seq = f.variation;  // may be empty (a deletion)
      if (f.has_range)
      {
        if (f.begin == 0 || f.end == 0)
        {
          OPENMS_LOG_WARN << "UniPEFF: " << accession << " sequence variant with unknown range; omitted." << std::endl;
          return;
        }
        a.complex_variants.push_back(VariantComplexItem{f.begin, f.end, new_seq});
      }
      else if (f.has_position)
      {
        if (f.position == 0)
        {
          OPENMS_LOG_WARN << "UniPEFF: " << accession << " sequence variant with unknown position; omitted." << std::endl;
          return;
        }
        if (new_seq.size() == 1 && isResidueCode(new_seq[0]))
        {
          a.simple_variants.push_back(VariantSimpleItem{f.position, new_seq[0]});
        }
        else
        {
          a.complex_variants.push_back(VariantComplexItem{f.position, f.position, new_seq});
        }
      }
    }
  }

  /// Merge halfcys into regular_mods producing @c a.mods in UniPEFF sort order:
  /// (position == 0 sorts LAST, then by position ascending, then regular before halfcys on ties).
  /// Updates @c a.disulfides to reference indices into the merged list.
  void mergeHalfCystines(EntryAnnotations& a)
  {
    struct Indexed { ModResItem m; bool from_halfcys; size_t orig_idx; };
    std::vector<Indexed> all;
    all.reserve(a.regular_mods.size() + a.halfcys.size());
    for (size_t i = 0; i < a.regular_mods.size(); ++i)
    {
      all.push_back(Indexed{a.regular_mods[i], false, i});
    }
    for (size_t i = 0; i < a.halfcys.size(); ++i)
    {
      all.push_back(Indexed{a.halfcys[i], true, i});
    }
    // Stable sort with the UniPEFF tie-break.
    std::stable_sort(all.begin(), all.end(), [](const Indexed& x, const Indexed& y) {
      const bool xu = (x.m.position == 0);
      const bool yu = (y.m.position == 0);
      if (xu != yu) return !xu;  // !unknown < unknown -> known first
      if (x.m.position != y.m.position) return x.m.position < y.m.position;
      // Regular mods (is_halfcys=false) sort before halfcys (true).
      return !x.m.is_halfcys && y.m.is_halfcys;
    });

    // Build new merged index for each (was_halfcys, orig_idx) so disulfide pairs can be re-pointed.
    std::vector<size_t> halfcys_new_idx(a.halfcys.size(), 0);
    a.mods.clear();
    a.mods.reserve(all.size());
    for (size_t i = 0; i < all.size(); ++i)
    {
      if (all[i].from_halfcys) halfcys_new_idx[all[i].orig_idx] = a.mods.size();
      a.mods.push_back(all[i].m);
    }
    for (auto& d : a.disulfides)
    {
      d.idx_a = halfcys_new_idx[d.idx_a];
      d.idx_b = halfcys_new_idx[d.idx_b];
    }
    // After merge the source lists are no longer authoritative.
    a.regular_mods.clear();
    a.halfcys.clear();
  }

  // ──────────────────────────────────────────────────────────────────
  // Annotation-ID assignment
  // ──────────────────────────────────────────────────────────────────

  /// Option B (-annotation_identifiers): walk @p a in UniPEFF's emit order assigning
  /// monotonically increasing 1-based IDs to every annotation tuple (ModRes split into
  /// PSI/Unimod/generic buckets, then VariantSimple, VariantComplex, Processed, finally
  /// DisulfideBond). The disulfide pair items receive IDs of their own; @c idx_a/@c idx_b
  /// already point to half-cystine ModRes items whose IDs have just been stamped.
  /// Tuples skipped here (unknown/out-of-bounds positions) keep @c kNoId — the writer
  /// skips exactly the same tuples, so every emitted tuple carries an id.
  void assignAnnotationIds(EntryAnnotations& a, const std::string& base_sequence,
                           bool emit_processed, bool emit_aa_mods, bool emit_variants)
  {
    uint32_t next_id = 1;
    if (emit_aa_mods)
    {
      // PSI-MOD bucket: ModRes whose PTM has a PSI-MOD accession.
      for (auto& m : a.mods)
      {
        if (m.ptm != nullptr && !m.ptm->psi_mod.empty()) m.annotation_id = next_id++;
      }
      // Unimod bucket: ModRes with Unimod but no PSI-MOD.
      for (auto& m : a.mods)
      {
        if (m.ptm != nullptr && m.ptm->psi_mod.empty() && !m.ptm->unimod.empty()) m.annotation_id = next_id++;
      }
      // Generic bucket: ModRes with no CV accession (or no PTM at all).
      for (auto& m : a.mods)
      {
        if (m.ptm == nullptr || (m.ptm->psi_mod.empty() && m.ptm->unimod.empty())) m.annotation_id = next_id++;
      }
    }

    if (emit_variants)
    {
      // VariantSimple: real simple list, then complex entries that get demoted to simple.
      const int seq_len = static_cast<int>(base_sequence.size());
      for (auto& v : a.simple_variants)
      {
        if (v.position == 0) continue;
        if (seq_len > 0 && v.position > seq_len) continue;
        v.annotation_id = next_id++;
      }
      for (auto& v : a.complex_variants)
      {
        if (v.begin != 0 && v.begin == v.end && v.new_seq.size() == 1 && isResidueCode(v.new_seq[0])
            && !(seq_len > 0 && v.begin > seq_len))
        {
          v.annotation_id = next_id++;
        }
      }
      for (auto& v : a.complex_variants)
      {
        if (v.begin == 0 || v.end == 0) continue;
        if (v.begin > v.end || (seq_len > 0 && (v.begin > seq_len || v.end > seq_len))) continue;
        if (v.begin == v.end && v.new_seq.size() == 1 && isResidueCode(v.new_seq[0])) continue;  // demoted
        v.annotation_id = next_id++;
      }
    }

    if (emit_processed)
    {
      const int seq_len = static_cast<int>(base_sequence.size());
      for (auto& p : a.processed)
      {
        if (p.begin == 0 || p.end == 0) continue;
        if (p.begin > p.end || (seq_len > 0 && (p.begin > seq_len || p.end > seq_len))) continue;
        p.annotation_id = next_id++;
      }
    }

    if (emit_aa_mods)
    {
      // Disulfide tuples continue the id sequence; only valid pairs whose endpoints
      // exist in the merged mod list get one (same guards as the writer).
      for (auto& d : a.disulfides)
      {
        if (!d.valid) continue;
        if (d.idx_a >= a.mods.size() || d.idx_b >= a.mods.size()) continue;
        d.annotation_id = next_id++;
      }
    }
  }

  /// Default mode (no -annotation_identifiers): implement the selective labeling scheme
  /// from issue #9829 — the k-th (1-based) documented intrachain disulfide labels its
  /// begin half-cystine 2k-1 and its end half-cystine 2k, and the \DisulfideBond tuple
  /// itself is labeled k, referencing those two ids. All other annotations (including
  /// lone interchain half-cystines from single-<position> features) stay unlabeled.
  void assignDisulfideLabels(EntryAnnotations& a)
  {
    uint32_t next_bond = 1;
    for (auto& d : a.disulfides)
    {
      if (!d.valid) continue;
      if (d.idx_a >= a.mods.size() || d.idx_b >= a.mods.size()) continue;
      d.annotation_id = next_bond++;
      a.mods[d.idx_a].annotation_id = 2 * d.annotation_id - 1;
      a.mods[d.idx_b].annotation_id = 2 * d.annotation_id;
    }
  }

  // ──────────────────────────────────────────────────────────────────
  // PEFF text writer (byte-exact UniPEFF format)
  // ──────────────────────────────────────────────────────────────────

  struct PeffHeader
  {
    std::string db_name{"UniProtKB"};
    std::string prefix{"sp"};
    std::string db_description;
    std::string db_source{"https://www.uniprot.org"};
    std::string db_version{"unknown"};
    std::string db_date;
    int number_of_entries{0};
    bool has_annotation_identifiers{false};
  };

  void writeFileDescriptionBlock(std::ostream& out)
  {
    out << "# PEFF 1.0\n"
        << "# GeneralComment=Generated by UniPEFF\n"
        << "# //\n";
  }

  void writeDbDescriptionBlock(std::ostream& out, const PeffHeader& h)
  {
    out << "# DbName=" << h.db_name << "\n"
        << "# Prefix=" << h.prefix << "\n";
    if (!h.db_description.empty()) out << "# DbDescription=" << h.db_description << "\n";
    out << "# Decoy=false\n"
        << "# DbSource=" << h.db_source << "\n"
        << "# DbVersion=" << h.db_version << "\n";
    if (!h.db_date.empty()) out << "# DbDate=" << h.db_date << "\n";
    out << "# NumberOfEntries=" << h.number_of_entries << "\n"
        << "# SequenceType=AA\n";
    if (h.has_annotation_identifiers) out << "# HasAnnotationIdentifiers=true\n";
    out << "# //\n";
  }

  /// Tracks how often the OBO-name lookup fell back to the UniProt ptmlist ID.
  /// Total occurrences vs. distinct missed accessions — "2886 occurrences across
  /// 30 distinct accessions" tells the user the OBO is ~30 terms out of date
  /// (actionable refresh); "2886 across 2886" tells them the OBO is broken.
  struct OboFallbackTracker
  {
    size_t total{0};
    std::set<std::string> distinct;
  };

  /// Resolve the canonical display name for one ModRes via OBO maps (PEFF
  /// requires the OBO "name:" field for ModResPsi/ModResUnimod; PEFF allows
  /// a generic name for plain ModRes). On a miss we fall back to the UniProt
  /// ptmlist ID and record both the occurrence and the distinct accession.
  std::string resolveDisplayName(const ModResItem& m, const std::string& accession,
                                 const OboNameMap& psi, const OboNameMap& uni,
                                 OboFallbackTracker& tracker)
  {
    if (!accession.empty())
    {
      const OboNameMap* maps[2] = { &psi, &uni };
      for (const OboNameMap* om : maps)
      {
        auto it = om->find(accession);
        if (it != om->end()) return it->second;
      }
      ++tracker.total;
      tracker.distinct.insert(accession);
    }
    return m.name;  // UniProt PTM ID
  }

  /// Write the descriptor + sequence for one entry. OBO-name fallbacks taken during
  /// emission are accumulated into @p tracker so the final report can distinguish
  /// "N occurrences across M distinct accessions".
  void writePeffEntry(std::ostream& out, const UniProtEntry& e, EntryAnnotations& a,
                      const std::string& prefix,
                      const OboNameMap& psi_obo, const OboNameMap& unimod_obo,
                      OboFallbackTracker& tracker,
                      bool emit_processed, bool emit_aa_mods, bool emit_variants)
  {
    const int seq_len = static_cast<int>(e.sequence.size());

    out << ">" << prefix << ":" << (e.accession.empty() ? std::string() : e.accession);

    auto tag_kv = [&](const std::string& key, const std::string& value) {
      out << " \\" << key << "=" << value;
    };

    if (!e.full_name.empty())
    {
      // Emit \PName as a one-element list: (value). The PEFF spec permits either form
      // (scalar OR list) for single-value keys, but the list form is unambiguous to parse
      // (a scalar value containing a balanced parenthetical can otherwise be silently
      // truncated by a paren-list-greedy reader). This diverges from the upstream C#
      // UniPEFF's scalar emission but stays valid PEFF 1.0.
      tag_kv("PName", std::string("(") + escapePeff(e.full_name) + ")");
    }
    if (!e.primary_gene.empty()) tag_kv("GName", escapePeff(e.primary_gene));
    if (!e.ncbi_tax_id.empty())  tag_kv("NcbiTaxId", escapePeff(e.ncbi_tax_id));
    if (!e.tax_name.empty())     tag_kv("TaxName", escapePeff(e.tax_name));
    if (!e.sequence.empty())     tag_kv("Length", std::to_string(seq_len));
    if (!e.sequence_version.empty()) tag_kv("SV", escapePeff(e.sequence_version));
    if (!e.entry_version.empty())    tag_kv("EV", escapePeff(e.entry_version));
    const std::string pe = proteinExistenceCode(e.protein_existence);
    if (!pe.empty()) tag_kv("PE", pe);
    if (!e.name.empty()) tag_kv("ID", escapePeff(e.name));
    if (!e.alt_accessions.empty())
    {
      std::string val;
      for (const auto& ac : e.alt_accessions)
      {
        val.push_back('(');
        val += escapePeff(ac);
        val.push_back(')');
      }
      tag_kv("AltAC", val);
    }

    // Modifications grouped by accession kind. IDs were stamped beforehand —
    // by assignAnnotationIds() (Option B, every tuple, in this exact bucket
    // order) or by assignDisulfideLabels() (default mode, paired half-cystines
    // only) — and a tuple prints an "id:" prefix iff it carries one.
    auto writeMods = [&](const std::string& key,
                         std::function<bool(const ModResItem&)> belongs,
                         std::function<std::string(const ModResItem&)> accession_of)
    {
      std::string val;
      for (const ModResItem& m : a.mods)
      {
        if (!belongs(m)) continue;
        val.push_back('(');
        if (m.annotation_id != kNoId) { val += std::to_string(m.annotation_id); val.push_back(':'); }
        val += positionOrUnknown(m.position);
        const std::string acc = accession_of(m);
        val.push_back('|');
        val += escapePeff(acc);
        const std::string display = resolveDisplayName(m, acc, psi_obo, unimod_obo, tracker);
        if (!display.empty())
        {
          val.push_back('|');
          val += escapePeff(display);
        }
        val.push_back(')');
      }
      if (!val.empty()) tag_kv(key, val);
    };

    if (emit_aa_mods)
    {
      writeMods("ModResPsi",
                [](const ModResItem& m) { return m.ptm != nullptr && !m.ptm->psi_mod.empty(); },
                [](const ModResItem& m) { return m.ptm->psi_mod; });
      writeMods("ModResUnimod",
                [](const ModResItem& m) { return m.ptm != nullptr && m.ptm->psi_mod.empty() && !m.ptm->unimod.empty(); },
                [](const ModResItem& m) { return std::string("UNIMOD:") + m.ptm->unimod; });
      writeMods("ModRes",
                [](const ModResItem& m) { return m.ptm == nullptr || (m.ptm->psi_mod.empty() && m.ptm->unimod.empty()); },
                [](const ModResItem&) { return std::string(); });
    }

    if (emit_variants)
    {
      // VariantSimple: real simples + demoted complex-as-simple.
      {
        std::string val;
        for (const auto& v : a.simple_variants)
        {
          if (v.position == 0) continue;
          if (seq_len > 0 && v.position > seq_len)
          {
            OPENMS_LOG_WARN << "UniPEFF: " << e.accession << " VariantSimple position " << v.position
                            << " exceeds sequence length " << seq_len << "; omitted." << std::endl;
            continue;
          }
          val.push_back('(');
          if (v.annotation_id != kNoId) { val += std::to_string(v.annotation_id); val.push_back(':'); }
          val += std::to_string(v.position);
          val.push_back('|');
          val += escapePeff(std::string(1, v.new_aa));
          val.push_back(')');
        }
        for (const auto& v : a.complex_variants)
        {
          if (v.begin != 0 && v.begin == v.end && v.new_seq.size() == 1 && isResidueCode(v.new_seq[0])
              && !(seq_len > 0 && v.begin > seq_len))
          {
            val.push_back('(');
            if (v.annotation_id != kNoId) { val += std::to_string(v.annotation_id); val.push_back(':'); }
            val += std::to_string(v.begin);
            val.push_back('|');
            val += escapePeff(v.new_seq);
            val.push_back(')');
          }
        }
        if (!val.empty()) tag_kv("VariantSimple", val);
      }
      // VariantComplex: range substitutions/deletions/insertions (excluding demoted).
      {
        std::string val;
        for (const auto& v : a.complex_variants)
        {
          if (v.begin == 0 || v.end == 0) continue;
          if (v.begin > v.end || (seq_len > 0 && (v.begin > seq_len || v.end > seq_len)))
          {
            OPENMS_LOG_WARN << "UniPEFF: " << e.accession << " VariantComplex " << v.begin << "-" << v.end
                            << " is out of bounds (length " << seq_len << "); omitted." << std::endl;
            continue;
          }
          if (v.begin == v.end && v.new_seq.size() == 1 && isResidueCode(v.new_seq[0])) continue;
          val.push_back('(');
          if (v.annotation_id != kNoId) { val += std::to_string(v.annotation_id); val.push_back(':'); }
          val += std::to_string(v.begin);
          val.push_back('|');
          val += std::to_string(v.end);
          val.push_back('|');
          val += escapePeff(v.new_seq);
          val.push_back(')');
        }
        if (!val.empty()) tag_kv("VariantComplex", val);
      }
    }

    if (emit_processed)
    {
      std::string val;
      for (const auto& p : a.processed)
      {
        if (p.begin == 0 || p.end == 0) continue;
        if (p.begin > p.end || (seq_len > 0 && (p.begin > seq_len || p.end > seq_len)))
        {
          OPENMS_LOG_WARN << "UniPEFF: " << e.accession << " Processed " << p.begin << "-" << p.end
                          << " is out of bounds (length " << seq_len << "); omitted." << std::endl;
          continue;
        }
        val.push_back('(');
        if (p.annotation_id != kNoId) { val += std::to_string(p.annotation_id); val.push_back(':'); }
        val += std::to_string(p.begin);
        val.push_back('|');
        val += std::to_string(p.end);
        val.push_back('|');
        val += escapePeff(p.cv);
        val.push_back('|');
        val += escapePeff(p.type_name);
        val.push_back(')');
      }
      if (!val.empty()) tag_kv("Processed", val);
    }

    if (emit_aa_mods)
    {
      // \DisulfideBond=(bond_id:idA,idB) — one tuple per documented intrachain bond,
      // referencing the ids its two half-cystines carry inside \ModResPsi. The ids were
      // stamped by whichever assignment pass ran (issue #9829 selective labels in the
      // default mode; global annotation identifiers in Option B); a pair without an id
      // (interchain/invalid) is not emitted.
      std::string val;
      for (const auto& d : a.disulfides)
      {
        if (d.annotation_id == kNoId) continue;
        val.push_back('(');
        val += std::to_string(d.annotation_id);
        val.push_back(':');
        val += std::to_string(a.mods[d.idx_a].annotation_id);
        val.push_back(',');
        val += std::to_string(a.mods[d.idx_b].annotation_id);
        val.push_back(')');
      }
      if (!val.empty()) tag_kv("DisulfideBond", val);
    }

    out << "\n";

    // Sequence, 60 chars per line.
    constexpr int kWidth = 60;
    const std::string& seq = e.sequence;
    for (int offset = 0; offset < seq_len; offset += kWidth)
    {
      const int take = std::min(kWidth, seq_len - offset);
      out.write(seq.data() + offset, take);
      out << "\n";
    }
  }

  // ──────────────────────────────────────────────────────────────────
  // Whole-file emission orchestration
  // ──────────────────────────────────────────────────────────────────

  struct PreparedEntry
  {
    UniProtEntry source;
    EntryAnnotations annotations;
    std::string prefix;       ///< sp/tr or user override
  };

  /// Prepare one entry: classify features, merge halfcys, (Option B) stamp IDs.
  PreparedEntry prepareEntry(UniProtEntry&& e, const PtmMap& ptms, const std::string& prefix_override,
                             bool option_b, bool record_processing, bool record_aa_mods,
                             bool record_variants)
  {
    PreparedEntry pe;
    pe.source = std::move(e);
    pe.prefix = prefix_override.empty() ? prefixForDataset(pe.source.dataset) : prefix_override;

    for (const auto& f : pe.source.features)
    {
      classifyAndAppend(f, ptms, pe.annotations, record_processing, record_aa_mods,
                        record_variants, pe.source.accession);
    }
    mergeHalfCystines(pe.annotations);

    if (option_b)
    {
      assignAnnotationIds(pe.annotations, pe.source.sequence,
                          record_processing, record_aa_mods, record_variants);
    }
    else
    {
      // No-op when the entry has no documented disulfide pairs
      // (including under -omit_amino_acid_modifications).
      assignDisulfideLabels(pe.annotations);
    }
    return pe;
  }

  bool isWritableEntry(const UniProtEntry& e)
  {
    return !e.accession.empty() && !e.sequence.empty();
  }

} // namespace

// ──────────────────────────────────────────────────────────────────
// TOPP wrapper
// ──────────────────────────────────────────────────────────────────

class TOPPUniPEFF :
  public TOPPBase
{
public:
  TOPPUniPEFF() :
    TOPPBase("UniPEFF", "Convert a UniProtKB XML protein database to PEFF 1.0 with rich annotations.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input UniProtKB XML file (plain or gzip; gzip is auto-detected).");
    setValidFormats_("in", {"xml"});  // .xml.gz is accepted: FileHandler strips .gz and re-checks

    registerOutputFile_("out", "<file>", "", "Output PEFF 1.0 file.");
    setValidFormats_("out", {"peff"});

    registerInputFile_("ptmlist", "<file>", "", "UniProt ptmlist.txt; defaults to the bundled snapshot.", false);
    setValidFormats_("ptmlist", {"txt"});

    registerInputFile_("psimod_obo", "<file>", "", "PSI-MOD.obo for canonical modification names; defaults to the bundled OpenMS PSI-MOD.obo.", false);
    setValidFormats_("psimod_obo", {"obo"});

    registerInputFile_("unimod_obo", "<file>", "", "Optional unimod.obo for canonical Unimod names; if absent, names fall back to the UniProt ptmlist IDs.", false);
    setValidFormats_("unimod_obo", {"obo"});

    registerStringOption_("prefix", "<string>", "", "Force a single PEFF prefix for every entry (e.g. 'sp'); if empty, sp/tr is derived from the UniProt dataset.", false);
    registerStringOption_("dbversion", "<string>", "unknown", "Value for the mandatory '# DbVersion=' PEFF header line.", false);

    registerFlag_("annotation_identifiers", "Emit PEFF Option B: assign a global sequential id: prefix to every annotation tuple, referenced by \\DisulfideBond. By default only half-cystines in documented disulfide bonds get (1-based, bond-ordered) ids.");
    registerFlag_("omit_molecular_processing", "Skip the \\Processed annotations (initiator methionine, signal/transit peptide, propeptide, chain).");
    registerFlag_("omit_amino_acid_modifications", "Skip \\ModResPsi / \\ModResUnimod / \\ModRes and \\DisulfideBond; ptmlist is not read.");
    registerFlag_("omit_sequence_variations", "Skip \\VariantSimple and \\VariantComplex annotations.");
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string in_file        = getStringOption_("in");
    const std::string out_file       = getStringOption_("out");
    const std::string prefix_override = getStringOption_("prefix");
    const std::string dbversion      = getStringOption_("dbversion");
    const bool option_b              = getFlag_("annotation_identifiers");
    const bool omit_proc             = getFlag_("omit_molecular_processing");
    const bool omit_aa               = getFlag_("omit_amino_acid_modifications");
    const bool omit_var              = getFlag_("omit_sequence_variations");

    const bool record_processing = !omit_proc;
    const bool record_aa_mods    = !omit_aa;
    const bool record_variants   = !omit_var;

    // PEFF is a plain-text format; we do not (yet) write directly into a compressed
    // container, and TOPPBase's format validation accepts `.peff.gz` etc. via its
    // recursive suffix stripping. Reject compressed suffixes explicitly rather than
    // silently writing uncompressed bytes under a compressed file name.
    if (out_file.ends_with(".gz") || out_file.ends_with(".bz2") || out_file.ends_with(".zip"))
    {
      writeLogError_("UniPEFF: compressed PEFF outputs are not supported; pass a plain '.peff' filename.");
      return ILLEGAL_PARAMETERS;
    }

    // Resolve auxiliary files.
    std::string ptmlist_file = getStringOption_("ptmlist");
    if (record_aa_mods && ptmlist_file.empty())
    {
      try { ptmlist_file = File::find("CHEMISTRY/UniProt_ptmlist.txt"); }
      catch (...) {
        writeLogError_("UniPEFF: cannot locate bundled CHEMISTRY/UniProt_ptmlist.txt and -ptmlist not given.");
        return INPUT_FILE_NOT_FOUND;
      }
    }

    std::string psimod_obo_file = getStringOption_("psimod_obo");
    if (record_aa_mods && psimod_obo_file.empty())
    {
      try { psimod_obo_file = File::find("CHEMISTRY/PSI-MOD.obo"); }
      catch (...) { /* tolerated: name lookup falls back to UniProt IDs */ }
    }
    const std::string unimod_obo_file = getStringOption_("unimod_obo");  // empty by default

    // Load auxiliary data.
    PtmMap ptms;
    OboNameMap psi_obo;
    OboNameMap unimod_obo;
    if (record_aa_mods)
    {
      OPENMS_LOG_INFO << "UniPEFF: loading ptmlist from " << ptmlist_file << std::endl;
      ptms = loadPtmlist(ptmlist_file);
      if (!psimod_obo_file.empty())
      {
        OPENMS_LOG_INFO << "UniPEFF: loading PSI-MOD OBO from " << psimod_obo_file << std::endl;
        psi_obo = loadOboNames(psimod_obo_file);
      }
      if (!unimod_obo_file.empty())
      {
        OPENMS_LOG_INFO << "UniPEFF: loading Unimod OBO from " << unimod_obo_file << std::endl;
        unimod_obo = loadOboNames(unimod_obo_file);
      }
    }

    // Streaming pipeline (bounded memory: one entry in flight at a time):
    //   * For each parsed UniProtEntry, write its serialized PEFF text into a
    //     per-prefix temporary file as we go. Per-prefix entry counts are
    //     accumulated alongside so we can emit the right `# NumberOfEntries=`
    //     in each DB-description block.
    //   * After the parse finishes, assemble the final file as
    //         <file header> + <DB blocks per prefix> + <concatenated per-prefix temp files>
    //     and remove the temp files.
    //
    // PEFF requires the entire header section to precede the entries, and
    // `# NumberOfEntries=` cannot be known until parsing is complete; the
    // per-prefix spool keeps memory use independent of input size while
    // preserving entry ordering within each prefix (which the byte-exact
    // golden tests depend on).
    struct PrefixSpool
    {
      std::string path;
      std::ofstream out;
      int count{0};
    };
    std::map<std::string, PrefixSpool> spools;
    std::vector<std::string> prefixes;  // first-seen order
    size_t skipped_no_accession = 0;
    size_t skipped_no_sequence = 0;

    auto cleanup_spools = [&]() {
      for (auto& [_, s] : spools)
      {
        if (s.out.is_open()) s.out.close();
        std::remove(s.path.c_str());
      }
    };

    bool spool_open_failed = false;
    std::string spool_open_failed_path;
    OboFallbackTracker fallback_tracker;
    {
      UniProtXMLFile xml;
      xml.loadStreaming(in_file, [&](UniProtEntry&& entry) {
        if (spool_open_failed) return;  // abort downstream work once any spool open fails
        if (!isWritableEntry(entry))
        {
          if (entry.accession.empty()) ++skipped_no_accession;
          else
          {
            OPENMS_LOG_WARN << "UniPEFF: skipping entry " << entry.accession << " with no sequence." << std::endl;
            ++skipped_no_sequence;
          }
          return;
        }
        PreparedEntry pe = prepareEntry(std::move(entry), ptms, prefix_override, option_b,
                                        record_processing, record_aa_mods, record_variants);
        auto it = spools.find(pe.prefix);
        if (it == spools.end())
        {
          // Open the spool BEFORE registering the prefix — on failure we must not leave a
          // dangling prefix that would later default-construct a 0-entry DB block from
          // spools[p] and silently drop everything we tried to write for it.
          PrefixSpool s;
          s.path = out_file + "." + pe.prefix + ".tmp";
          s.out.open(s.path, std::ios::binary | std::ios::trunc);
          if (!s.out)
          {
            spool_open_failed = true;
            spool_open_failed_path = s.path;
            return;
          }
          it = spools.emplace(pe.prefix, std::move(s)).first;
          prefixes.push_back(pe.prefix);
        }
        writePeffEntry(it->second.out, pe.source, pe.annotations, pe.prefix,
                       psi_obo, unimod_obo, fallback_tracker,
                       record_processing, record_aa_mods, record_variants);
        ++it->second.count;
      });
    }
    if (spool_open_failed)
    {
      cleanup_spools();
      writeLogError_("UniPEFF: cannot open spool file '" + spool_open_failed_path +
                     "'. Ensure the output directory exists and is writable.");
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    if (skipped_no_accession > 0)
    {
      OPENMS_LOG_WARN << "UniPEFF: skipped " << skipped_no_accession << " entries with no accession." << std::endl;
    }

    int writable = 0;
    for (const auto& [_, s] : spools) writable += s.count;
    OPENMS_LOG_INFO << "UniPEFF: parsed " << (writable + static_cast<int>(skipped_no_accession + skipped_no_sequence))
                    << " UniProtKB entries (" << writable << " writable)." << std::endl;
    if (writable == 0)
    {
      cleanup_spools();
      writeLogError_("UniPEFF: no writable entries; refusing to write an empty PEFF file.");
      return PARSE_ERROR;
    }

    // Close spools so their data is flushed before we read them back.
    for (auto& [_, s] : spools)
    {
      if (s.out.is_open()) s.out.close();
    }

    std::ofstream out(out_file, std::ios::binary | std::ios::trunc);
    if (!out)
    {
      cleanup_spools();
      writeLogError_("UniPEFF: cannot open output file '" + out_file + "'.");
      return CANNOT_WRITE_OUTPUT_FILE;
    }

    writeFileDescriptionBlock(out);

    // PEFF requires every database-description block to appear before any entry.
    for (const std::string& p : prefixes)
    {
      PeffHeader h;
      h.prefix = p;
      h.db_version = dbversion;
      h.has_annotation_identifiers = option_b;
      h.number_of_entries = spools[p].count;
      writeDbDescriptionBlock(out, h);
    }

    // Now stream the per-prefix entry text in first-seen prefix order.
    for (const std::string& p : prefixes)
    {
      std::ifstream tmp(spools[p].path, std::ios::binary);
      if (!tmp)
      {
        cleanup_spools();
        writeLogError_("UniPEFF: cannot re-read spool file '" + spools[p].path + "'.");
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      out << tmp.rdbuf();
    }
    out.close();
    cleanup_spools();

    if (fallback_tracker.total > 0)
    {
      OPENMS_LOG_WARN << "UniPEFF: " << fallback_tracker.total << " modification name(s) across "
                      << fallback_tracker.distinct.size() << " distinct accession(s) fell back to UniProt ptmlist IDs "
                         "(no OBO 'name:' entry found). Provide an updated -psimod_obo / -unimod_obo for strict PEFF "
                         "conformance." << std::endl;
    }
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPUniPEFF tool;
  return tool.main(argc, argv);
}

/// @endcond
