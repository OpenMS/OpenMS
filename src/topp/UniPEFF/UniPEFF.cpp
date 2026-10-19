// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include "UniPEFFEntries.h"
#include "UniProtXMLFile.h"
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <regex>
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
inside \\ModResPsi, for every bond whose two half-cystines UniProt locates in
this sequence (<code>\<begin\></code>/<code>\<end\></code>). These are mostly
intrachain bonds, but also bonds between chains cleaved from the same
precursor (e.g. insulin's "Interchain (between B and A chains)"). A bond given
as a single <code>\<position\></code> (partner cysteine in another molecule) or
with an endpoint beyond the sequence end gets no connectivity; its in-range
half-cystines are written as plain modifications (see below for the rest). By
default only the reported bonds and their half-cystines are labeled, in bond
order: of K reported bonds, bond k (counting from 0) labels its two
half-cystines 2k and 2k+1 and is itself labeled 2K+k, so the ids 0..3K-1 are
unique within the entry, as PEFF requires (see issue 9829). With
@c -annotation_identifiers (PEFF "Option B") every
annotation tuple instead carries a global sequential id (0, 1, 2, ...) and
\\DisulfideBond references those ids. Each database block with ids declares
<code># HasAnnotationIdentifiers=true</code>.

Annotations positioned beyond the end of the sequence (malformed input), i.e.
modifications, variants and processed regions, are omitted with a warning,
since PEFF treats such a position as an error in the file.

Features whose <code>\<location sequence="..."\></code> names another isoform
carry coordinates of that isoform, not of the entry's canonical sequence. They
are never applied to the canonical entry. If that isoform is written as an entry
of its own (see below), they become annotations of that entry, in its
coordinates; otherwise they are not applied. Summary log lines count both kinds
among the features that would otherwise have become annotations. A location
naming the "displayed" isoform refers to the canonical sequence.

Modification accession lookup uses UniProt's <em>ptmlist.txt</em> (a snapshot
is bundled under @c share/OpenMS/CHEMISTRY/UniProt_ptmlist.txt); override
with @c -ptmlist. Canonical OBO names come from <em>PSI-MOD.obo</em> (bundled)
and an optional <em>unimod.obo</em>; without them, names fall back to the
UniProt @c ptmlist.txt ID and a warning is printed.

UniProt alternative products (isoforms) are expanded into their own PEFF
entries: every isoform whose sequence is "described" by
<code>\<feature type="splice variant"\></code> records is reconstructed from
the canonical sequence (a referenced region without a replacement is deleted,
one with a <code>\<variation\></code> is substituted) and emitted directly
after its parent entry, e.g. <code>\>sp:P02768-2 \\PName=(Isoform 2 of
Albumin) ...</code>, with gene/taxonomy/mnemonic metadata inherited from the
parent (mirroring UniProt's isoform FASTA convention, which carries no
per-isoform PE/SV). The isoform flagged "displayed" is identical to the
canonical sequence and is not emitted again; isoforms typed "external" or
"not described" have no reconstructable sequence and are reported in the log.
An isoform entry carries the features UniProt annotates on that isoform
(<code>\<location sequence="..."\></code>, see above), with ids and
\\DisulfideBond connectivity as for canonical entries; the canonical entry's
annotations are not mapped onto the spliced sequence. Disable isoform expansion
with @c -omit_isoforms.

Both plain @c .xml and @c .xml.gz UniProt inputs are accepted (gzip is
auto-detected by the underlying parser).

The mandatory <code># DbVersion=</code> header line takes the value of
@c -dbversion. If that is empty, it is the download date UniProt puts into the
names of its files, e.g. @c 2025_07_05 for
<em>uniprotkb_proteome_UP000005640_2025_07_05.xml.gz</em>, or @c unknown if
the input file name contains no such date.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_UniPEFF.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_UniPEFF.html
*/

/// @cond TOPPCLASSES

namespace
{
  using namespace OpenMS::UniPEFFEntries;

  // ──────────────────────────────────────────────────────────────────
  // ptmlist.txt — UniProt PTM controlled vocabulary
  // ──────────────────────────────────────────────────────────────────

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
  // Isoform reconstruction (UniProt "alternative products")
  // ──────────────────────────────────────────────────────────────────

  // ──────────────────────────────────────────────────────────────────
  // Intermediate annotation structures (built per UniProtEntry)
  // ──────────────────────────────────────────────────────────────────

  /// True when any tuple of @p a carries an annotation identifier; the entry's database
  /// block must then declare HasAnnotationIdentifiers=true (PEFF 1.0, section 3.4.2).
  bool hasAnnotationIds(const EntryAnnotations& a)
  {
    auto has_id = [](const auto& item) { return item.annotation_id != kNoId; };
    return std::any_of(a.mods.begin(), a.mods.end(), has_id)
        || std::any_of(a.simple_variants.begin(), a.simple_variants.end(), has_id)
        || std::any_of(a.complex_variants.begin(), a.complex_variants.end(), has_id)
        || std::any_of(a.processed.begin(), a.processed.end(), has_id)
        || std::any_of(a.disulfides.begin(), a.disulfides.end(), has_id);
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

  /// Write @p seq wrapped at 60 characters per line.
  void writeWrappedSequence(std::ostream& out, const std::string& seq)
  {
    constexpr int kWidth = 60;
    const int seq_len = static_cast<int>(seq.size());
    for (int offset = 0; offset < seq_len; offset += kWidth)
    {
      const int take = std::min(kWidth, seq_len - offset);
      out.write(seq.data() + offset, take);
      out << "\n";
    }
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

  /// Write the annotation tags (\\ModResPsi ... \\DisulfideBond) of one entry, canonical or isoform.
  /// Positions beyond @p seq_len are omitted with a warning naming @p accession. OBO-name
  /// fallbacks taken here are accumulated into @p tracker.
  void writeAnnotationTags(std::ostream& out, const std::string& accession, int seq_len,
                           const EntryAnnotations& a,
                           const OboNameMap& psi_obo, const OboNameMap& unimod_obo,
                           OboFallbackTracker& tracker,
                           bool emit_processed, bool emit_aa_mods, bool emit_variants)
  {
    auto tag_kv = [&](const std::string& key, const std::string& value) {
      out << " \\" << key << "=" << value;
    };

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
        if (seq_len > 0 && m.position > seq_len)
        {
          OPENMS_LOG_WARN << "UniPEFF: " << accession << " " << key << " position " << m.position
                          << " exceeds sequence length " << seq_len << "; omitted.\n";
          continue;
        }
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
            OPENMS_LOG_WARN << "UniPEFF: " << accession << " VariantSimple position " << v.position
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
            OPENMS_LOG_WARN << "UniPEFF: " << accession << " VariantComplex " << v.begin << "-" << v.end
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
          OPENMS_LOG_WARN << "UniPEFF: " << accession << " Processed " << p.begin << "-" << p.end
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
      // \DisulfideBond=(bond_id:idA,idB) — one tuple per valid <begin>/<end> disulfide pair,
      // referencing the ids its two half-cystines carry inside \ModResPsi. The ids were
      // stamped by whichever assignment pass ran (issue #9829 selective labels in the
      // default mode; global annotation identifiers in Option B); a pair without an id
      // (endpoint beyond the sequence end) is not emitted. Single-<position> disulfide
      // features never form a pair.
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

    writeAnnotationTags(out, e.accession, seq_len, a, psi_obo, unimod_obo, tracker,
                        emit_processed, emit_aa_mods, emit_variants);

    out << "\n";
    writeWrappedSequence(out, e.sequence);
  }

  /// Write the descriptor + sequence for one reconstructed isoform entry.
  /// Mirrors UniProt's isoform FASTA convention: "Isoform <name> of <parent name>"
  /// with gene / taxonomy / mnemonic inherited from the parent and no SV/EV/PE
  /// (UniProt versions isoform sequences with the parent entry, not separately).
  /// Its annotations are the features UniProt annotates on this isoform
  /// (&lt;location sequence="..."&gt;), in the isoform's coordinates; the canonical
  /// entry's annotations are not mapped onto the spliced sequence.
  void writeIsoformPeffEntry(std::ostream& out, const UniProtEntry& parent, const IsoformEntry& iso,
                             const std::string& prefix,
                             const OboNameMap& psi_obo, const OboNameMap& unimod_obo,
                             OboFallbackTracker& tracker,
                             bool emit_processed, bool emit_aa_mods, bool emit_variants)
  {
    out << ">" << prefix << ":" << iso.accession;

    auto tag_kv = [&](const std::string& key, const std::string& value) {
      out << " \\" << key << "=" << value;
    };

    std::string pname = "Isoform " + (iso.name.empty() ? iso.accession : iso.name);
    if (!parent.full_name.empty()) pname += " of " + parent.full_name;
    tag_kv("PName", std::string("(") + escapePeff(pname) + ")");
    if (!parent.primary_gene.empty()) tag_kv("GName", escapePeff(parent.primary_gene));
    if (!parent.ncbi_tax_id.empty())  tag_kv("NcbiTaxId", escapePeff(parent.ncbi_tax_id));
    if (!parent.tax_name.empty())     tag_kv("TaxName", escapePeff(parent.tax_name));
    tag_kv("Length", std::to_string(iso.sequence.size()));
    if (!parent.name.empty()) tag_kv("ID", escapePeff(parent.name));

    writeAnnotationTags(out, iso.accession, static_cast<int>(iso.sequence.size()), iso.annotations,
                        psi_obo, unimod_obo, tracker, emit_processed, emit_aa_mods, emit_variants);

    out << "\n";
    writeWrappedSequence(out, iso.sequence);
  }

  // ──────────────────────────────────────────────────────────────────
  // Whole-file emission orchestration
  // ──────────────────────────────────────────────────────────────────

  bool isWritableEntry(const UniProtEntry& e)
  {
    return !e.accession.empty() && !e.sequence.empty();
  }

  /// UniProt names downloads after their date, e.g. "uniprotkb_proteome_UP000005640_2025_07_05.xml.gz".
  /// Returns the last date token (YYYY_MM_DD, delimited by '_', '-', '.' or the ends of the name) in the
  /// file name of @p path, or an empty string if it has none.
  std::string dateFromFileName(const std::string& path)
  {
    static const std::regex date_token(R"((?:^|[_.-])((?:19|20)\d{2}_(?:0[1-9]|1[0-2])_(?:0[1-9]|[12]\d|3[01]))(?=[_.-]|$))");
    const std::string name = File::basename(path);
    std::string date;
    for (std::sregex_iterator it(name.begin(), name.end(), date_token), end; it != end; ++it)
    {
      date = (*it)[1].str();
    }
    return date;
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
    registerStringOption_("dbversion", "<string>", "", "Value for the mandatory '# DbVersion=' PEFF header line. If empty, the download date in the input file name is used (e.g. '2025_07_05' for 'uniprotkb_proteome_UP000005640_2025_07_05.xml.gz'), or 'unknown' if it has none.", false);

    registerFlag_("annotation_identifiers", "Emit PEFF Option B: assign a global sequential id: prefix to every annotation tuple, referenced by \\DisulfideBond. By default only the \\DisulfideBond tuples and the half-cystines they reference get ids: of K bonds, bond k (counting from 0) labels its half-cystines 2k and 2k+1 and is itself labeled 2K+k.");
    registerFlag_("omit_molecular_processing", "Skip the \\Processed annotations (initiator methionine, signal/transit peptide, propeptide, chain).");
    registerFlag_("omit_amino_acid_modifications", "Skip \\ModResPsi / \\ModResUnimod / \\ModRes and \\DisulfideBond; ptmlist is not read.");
    registerFlag_("omit_sequence_variations", "Skip \\VariantSimple and \\VariantComplex annotations.");
    registerFlag_("omit_isoforms", "Skip the additional PEFF entries for UniProt isoforms (alternative products whose sequences are reconstructed from splice-variant features).");
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string in_file        = getStringOption_("in");
    const std::string out_file       = getStringOption_("out");
    const std::string prefix_override = getStringOption_("prefix");
    std::string dbversion            = getStringOption_("dbversion");
    const bool option_b              = getFlag_("annotation_identifiers");
    const bool omit_proc             = getFlag_("omit_molecular_processing");
    const bool omit_aa               = getFlag_("omit_amino_acid_modifications");
    const bool omit_var              = getFlag_("omit_sequence_variations");
    const bool omit_iso              = getFlag_("omit_isoforms");

    const bool record_processing = !omit_proc;
    const bool record_aa_mods    = !omit_aa;
    const bool record_variants   = !omit_var;
    const bool record_isoforms   = !omit_iso;

    if (dbversion.empty())
    {
      dbversion = dateFromFileName(in_file);
      if (dbversion.empty())
      {
        dbversion = "unknown";
        OPENMS_LOG_INFO << "UniPEFF: no -dbversion given and no date (YYYY_MM_DD) in the input file name; "
                           "writing '# DbVersion=unknown'." << std::endl;
      }
      else
      {
        OPENMS_LOG_INFO << "UniPEFF: no -dbversion given; using the date in the input file name: '# DbVersion="
                        << dbversion << "'." << std::endl;
      }
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
      bool has_annotation_ids{false};  ///< some entry carries ids -> HasAnnotationIdentifiers=true
    };
    std::map<std::string, PrefixSpool> spools;
    std::vector<std::string> prefixes;  // first-seen order
    size_t skipped_no_accession = 0;
    size_t skipped_no_sequence = 0;
    size_t isoforms_written = 0;
    SkippedIsoformCounts skipped_isoforms;

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
    IsoformScopedTracker isoform_scoped_tracker;
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
                                        record_processing, record_aa_mods, record_variants,
                                        record_isoforms, skipped_isoforms);
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
        isoform_scoped_tracker.add(pe.isoform_scoped_applied, pe.isoform_scoped_skipped);
        writePeffEntry(it->second.out, pe.source, pe.annotations, pe.prefix,
                       psi_obo, unimod_obo, fallback_tracker,
                       record_processing, record_aa_mods, record_variants);
        ++it->second.count;
        if (hasAnnotationIds(pe.annotations)) it->second.has_annotation_ids = true;
        // Isoform entries follow their parent immediately and count towards
        // this prefix's `# NumberOfEntries=` like any other entry.
        for (const IsoformEntry& iso : pe.isoforms)
        {
          writeIsoformPeffEntry(it->second.out, pe.source, iso, pe.prefix,
                                psi_obo, unimod_obo, fallback_tracker,
                                record_processing, record_aa_mods, record_variants);
          ++it->second.count;
          ++isoforms_written;
          if (hasAnnotationIds(iso.annotations)) it->second.has_annotation_ids = true;
        }
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

    int total_written = 0;
    for (const auto& [_, s] : spools) total_written += s.count;
    const int writable = total_written - static_cast<int>(isoforms_written);
    OPENMS_LOG_INFO << "UniPEFF: parsed " << (writable + static_cast<int>(skipped_no_accession + skipped_no_sequence))
                    << " UniProtKB entries (" << writable << " writable)." << std::endl;
    if (isoforms_written > 0)
    {
      OPENMS_LOG_INFO << "UniPEFF: wrote " << isoforms_written
                      << " isoform entries reconstructed from alternative products.\n";
    }
    if (!skipped_isoforms.empty())
    {
      size_t skipped_total = 0;
      std::string by_type;
      for (const auto& [type, n] : skipped_isoforms)
      {
        skipped_total += n;
        if (!by_type.empty()) by_type += ", ";
        by_type += std::to_string(n) + " '" + type + "'";
      }
      OPENMS_LOG_INFO << "UniPEFF: " << skipped_total << " isoform(s) without a reconstructable sequence were not emitted ("
                      << by_type << ").\n";
    }
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
      // Option B declares ids for every block; the default mode only for blocks whose
      // entries carry disulfide labels (PEFF allows the flag per database).
      h.has_annotation_identifiers = option_b || spools[p].has_annotation_ids;
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
    if (isoform_scoped_tracker.applied > 0)
    {
      OPENMS_LOG_INFO << "UniPEFF: " << isoform_scoped_tracker.applied << " feature(s) annotated on isoforms "
                         "(location/@sequence) were written to their isoform entries ("
                      << isoform_scoped_tracker.applied_entries << " entries affected).\n";
    }
    if (isoform_scoped_tracker.skipped > 0)
    {
      OPENMS_LOG_INFO << "UniPEFF: " << isoform_scoped_tracker.skipped << " feature(s) annotated on isoforms "
                         "that are not written (location/@sequence) were not applied ("
                      << isoform_scoped_tracker.skipped_entries << " entries affected).\n";
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
