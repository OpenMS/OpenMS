// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/VariantKnockoffs.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <algorithm>
#include <fstream>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_VariantDatabase VariantDatabase

@brief Adds single-residue protein variants, and knockoff entries that estimate their error, to a protein database.

<CENTER>
  <table>
    <tr>
      <th ALIGN = "center"> potential predecessor tools </th>
      <td VALIGN="middle" ROWSPAN=2> &rarr; VariantDatabase &rarr;</td>
      <th ALIGN = "center"> potential successor tools </th>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_DecoyDatabase, search engines </td>
    </tr>
  </table>
</CENTER>

Each variant of @p -variants becomes an entry of its own: the protein with the substitution, named
<tt>&lt;accession&gt;_VAR_&lt;ref&gt;&lt;position&gt;&lt;alt&gt;</tt>. With @p -knockoffs, up to that many <b>knockoff</b>
entries per variant are added as well: the same substitution at another position of the protein that holds the same
residue, in a matched sequence context and far from every listed variant (see VariantKnockoffs). A knockoff is almost
certainly absent from the sample, but a modification of the substitution's mass or an isotope error fits it as well as
it fits the real variant. Search the output (plus decoys) as usual, then pass it to @ref TOPP_FalseDiscoveryRate as
@p -PSM_variant_database: the knockoff PSMs estimate the false variant PSMs, and knockoff PSMs are never reported.

The variants are a tab-separated file with a header naming the columns @p accession, @p position (1-based), @p ref
and @p alt (one residue each), and optionally @p id (written to the description) and @p tier (the variant's prior:
@p sample, @p common, @p curated or @p catalogue; variants of different tiers get separate error estimates). The
accession is a FASTA identifier, or the middle field of a UniProt identifier (<tt>sp|P02769|ALBU_BOVIN</tt>).

Ported from bigbio/andes #173 and #180.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_VariantDatabase.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_VariantDatabase.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPVariantDatabase :
  public TOPPBase
{
public:
  TOPPVariantDatabase() :
    TOPPBase("VariantDatabase", "Adds single-residue protein variants and their knockoff entries to a protein database.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Protein database");
    setValidFormats_("in", {"fasta"});
    registerInputFile_("variants", "<file>", "", "Variants (tab-separated, header: accession, position, ref, alt[, id][, tier])");
    setValidFormats_("variants", {"tsv"});
    registerOutputFile_("out", "<file>", "", "Protein database with the variant and knockoff entries");
    setValidFormats_("out", {"fasta"});
    registerIntOption_("knockoffs", "<number>", 4, "Knockoff entries per variant, at most (0: none)", false);
    setMinInt_("knockoffs", 0);
    StringList all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<enzyme>", "Trypsin", "Enzyme of the search: a knockoff lies in another enzymatic peptide of a similar length", false);
    setValidStrings_("enzyme", all_enzymes);
    registerFlag_("only_variants", "Write only the variant and knockoff entries, not the proteins of 'in'");
  }

  /// read the variants; returns false (after logging) on a malformed file
  bool loadVariants_(const std::string& filename, std::vector<VariantKnockoffs::Variant>& variants)
  {
    std::ifstream in(filename);
    std::string line;
    if (!std::getline(in, line))
    {
      OPENMS_LOG_ERROR << "Error: " << filename << " is empty." << endl;
      return false;
    }
    std::map<std::string, Size> col;
    const auto header = ListUtils::create<std::string>(StringUtils::trimmed(line), '\t');
    for (Size i = 0; i < header.size(); ++i) col[StringUtils::toLowered(StringUtils::trimmed(header[i]))] = i;
    for (const char* required : {"accession", "position", "ref", "alt"})
    {
      if (col.count(required) == 0)
      {
        OPENMS_LOG_ERROR << "Error: the header of " << filename << " has no column '" << required << "'." << endl;
        return false;
      }
    }
    const std::vector<std::string> tiers = {"", "sample", "common", "curated", "catalogue"};
    Size line_no = 1;
    while (std::getline(in, line))
    {
      ++line_no;
      line = StringUtils::trimmed(line);
      if (line.empty() || line[0] == '#') continue;
      const auto fields = ListUtils::create<std::string>(line, '\t');
      auto field = [&](const std::string& name) -> std::string
      {
        auto it = col.find(name);
        return it == col.end() || it->second >= fields.size() ? std::string() : StringUtils::trimmed(fields[it->second]);
      };
      VariantKnockoffs::Variant v;
      v.accession = field("accession");
      const std::string ref = StringUtils::toUppered(field("ref"));
      const std::string alt = StringUtils::toUppered(field("alt"));
      v.id = field("id");
      v.tier = StringUtils::toLowered(field("tier"));
      int position = 0;
      try
      {
        position = std::stoi(field("position"));
      }
      catch (const std::exception&)
      {
        position = 0;
      }
      if (v.accession.empty() || position < 1 || ref.size() != 1 || alt.size() != 1)
      {
        OPENMS_LOG_WARN << "Line " << line_no << " of " << filename << ": not a single-residue substitution, skipped." << endl;
        continue;
      }
      if (std::find(tiers.begin(), tiers.end(), v.tier) == tiers.end())
      {
        OPENMS_LOG_ERROR << "Error: line " << line_no << " of " << filename << ": unknown tier '" << v.tier
                         << "' (sample, common, curated or catalogue)." << endl;
        return false;
      }
      v.position = static_cast<Size>(position - 1);
      v.ref = ref[0];
      v.alt = alt[0];
      variants.push_back(v);
    }
    return true;
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string in = getStringOption_("in");
    const std::string out = getStringOption_("out");

    std::vector<VariantKnockoffs::Variant> variants;
    if (!loadVariants_(getStringOption_("variants"), variants)) return INPUT_FILE_CORRUPT;

    std::vector<FASTAFile::FASTAEntry> db;
    FASTAFile().load(in, db);

    VariantKnockoffs::DatabaseStats stats;
    std::vector<FASTAFile::FASTAEntry> entries;
    try
    {
      entries = VariantKnockoffs::buildDatabase(db, variants, getStringOption_("enzyme"), static_cast<Size>(getIntOption_("knockoffs")), stats);
    }
    catch (Exception::InvalidValue& e)
    {
      OPENMS_LOG_ERROR << "Error: a variant substitutes an unknown residue: " << e.what() << endl;
      return INCOMPATIBLE_INPUT_DATA;
    }

    OPENMS_LOG_INFO << variants.size() << " variants read: " << stats.variant_entries << " variant entries, "
                    << stats.knockoff_entries << " knockoff entries (" << stats.variants_without_knockoff
                    << " variants with no position for a knockoff)." << endl;
    if (stats.unknown_protein + stats.reference_mismatch + stats.duplicates > 0)
    {
      OPENMS_LOG_WARN << "Skipped: " << stats.unknown_protein << " variants of a protein not in the database, "
                      << stats.reference_mismatch << " whose reference residue is not the protein's, "
                      << stats.duplicates << " listed more than once." << endl;
    }

    if (getFlag_("only_variants")) db.clear();
    db.insert(db.end(), entries.begin(), entries.end());
    FASTAFile().store(out, db);
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPVariantDatabase tool;
  return tool.main(argc, argv);
}

/// @endcond
