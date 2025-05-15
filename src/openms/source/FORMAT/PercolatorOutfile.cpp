// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PercolatorOutfile.h>

#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/FORMAT/CsvFile.h>

namespace OpenMS
{
  using namespace std;

  // initialize static variable:
  const std::string PercolatorOutfile::score_type_names[] =
    {"qvalue", "PEP", "score"};


  PercolatorOutfile::PercolatorOutfile() = default;


  enum PercolatorOutfile::ScoreType PercolatorOutfile::getScoreType(
    String score_type_name)
  {
    score_type_name.toLower();
    if ((score_type_name == "q-value") || (score_type_name == "qvalue") ||
        (score_type_name == "q value"))
    {
      return QVALUE;
    }
    if ((score_type_name == "pep") ||
        (score_type_name == "posterior error probability"))
    {
      return POSTERRPROB;
    }
    if (score_type_name == "score")
    {
      return SCORE;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Not a valid Percolator score type",
                                  score_type_name);
  }


  void PercolatorOutfile::resolveMisassignedNTermMods_(String& peptide) const
  {
    boost::regex re(R"(^[A-Z]\[(?<MOD1>-?\d+(\.\d+)?)\](\[(?<MOD2>-?\d+(\.\d+)?)\])?)");
    boost::smatch match;
    bool found = boost::regex_search(peptide, match, re);
    if (found && match["MOD1"].matched)
    {
      const ResidueModification* null = nullptr;
      vector<const ResidueModification*> maybe_nterm(2, null);
      String residue = peptide[0];
      String mod1 = match["MOD1"].str();
      double mass1 = mod1.toDouble();
      maybe_nterm[0] = ModificationsDB::getInstance()->
        getBestModificationByDiffMonoMass(mass1, 0.01, residue,
                                          ResidueModification::N_TERM);
      if (maybe_nterm[0] && !match["MOD2"].matched &&
          ((maybe_nterm[0]->getId() != "Carbamidomethyl") || (residue != "C")))
      { // only 1 mod, may be terminal -> assume terminal (unless it's CAM!):
        String replacement = ".(" + maybe_nterm[0]->getId() + ")" + residue;
        peptide = boost::regex_replace(peptide, re, replacement);
      }
      // only 1 mod, may not be terminal -> nothing to do
      else if (match["MOD2"].matched) // two mods
      {
        String mod2 = match["MOD2"].str();
        double mass2 = mod2.toDouble();
        maybe_nterm[1] = ModificationsDB::getInstance()->
          getBestModificationByDiffMonoMass(mass2, 0.01, residue,
                                            ResidueModification::N_TERM);
        if (maybe_nterm[0] && !maybe_nterm[1])
        { // first mod is terminal:
          String replacement = "(" + maybe_nterm[0]->getId() + ")" + residue +
            "[" + mod2 + "]";
          peptide = boost::regex_replace(peptide, re, replacement);
        }
        else if (maybe_nterm[1] && !maybe_nterm[0])
        { // second mod is terminal:
          String replacement = "(" + maybe_nterm[1]->getId() + ")" + residue +
            "[" + mod1 + "]";
          peptide = boost::regex_replace(peptide, re, replacement);
        }
        else // ambiguous cases
        {
          vector<const ResidueModification*> maybe_residue(2, null);
          maybe_residue[0] = ModificationsDB::getInstance()->
            getBestModificationByDiffMonoMass(mass1, 0.01, residue,
                                              ResidueModification::ANYWHERE);
          maybe_residue[1] = ModificationsDB::getInstance()->
            getBestModificationByDiffMonoMass(mass2, 0.01, residue,
                                              ResidueModification::ANYWHERE);
          if (maybe_nterm[0] && maybe_nterm[1]) // both mods may be terminal
          {
            if (maybe_residue[0] && !maybe_residue[1])
            { // first mod must be non-terminal -> second mod is terminal:
              String replacement = "(" + maybe_nterm[1]->getId() + ")" +
                residue + "[" + mod1 + "]";
              peptide = boost::regex_replace(peptide, re, replacement);
            }
            else if (maybe_residue[1] && !maybe_residue[0])
            { // second mod must be non-terminal -> first mod is terminal:
              String replacement = "(" + maybe_nterm[0]->getId() + ")" +
                residue + "[" + mod2 + "]";
              peptide = boost::regex_replace(peptide, re, replacement);
            }
            else // both mods may be terminal or non-terminal :-(
            { // arbitrarily assume first mod is terminal
              String replacement = "(" + maybe_nterm[0]->getId() + ")" +
                residue + "[" + mod2 + "]";
              peptide = boost::regex_replace(peptide, re, replacement);
            }
          }
          // if neither mod can be terminal, something is wrong -> let
          // AASequence deal with it
        }
      }
    }
  }


  void PercolatorOutfile::getPeptideSequence_(String peptide, AASequence& seq)
    const
  {
    // 'peptide' includes neighboring amino acids, e.g.: K.AAAR.A
    // but unclear to which protein neighboring AAs belong, so we ignore them:
    size_t len = peptide.size(), start = 0, count = std::string::npos;
    if (peptide[1] == '.')
    {
      start = 2;
    }
    if (peptide[len - 2] == '.')
    {
      count = len - start - 2;
    }
    peptide = peptide.substr(start, count);

    // re-format modifications:
    String unknown_mod = "[unknown]";
    if (peptide.hasSubstring(unknown_mod))
    {
      OPENMS_LOG_WARN << "Removing unknown modification(s) from peptide '" << peptide
               << "'" << endl;
      peptide.substitute(unknown_mod, "");
    }
    boost::regex re(R"(\[UNIMOD:(\d+)\])");
    std::string replacement = "(UniMod:$1)";
    peptide = boost::regex_replace(peptide, re, replacement);
    // search results from X! Tandem:
    // N-terminal mods may be wrongly assigned to the first residue; there may
    // be up to two mass shifts (one terminal, one residue) in random order!
    resolveMisassignedNTermMods_(peptide);
    // positive mass shifts are missing the "+":
    re.assign("\\[(\\d)");
    replacement = "[+$1";
    peptide = boost::regex_replace(peptide, re, replacement);

    seq = AASequence::fromString(peptide);
  }


  void PercolatorOutfile::load(const String& filename,
                               ProteinIdentification& proteins,
                               vector<PeptideIdentification>& peptides,
                               SpectrumMetaDataLookup& lookup,
                               enum ScoreType output_score)
  {
    SpectrumMetaDataLookup::MetaDataFlags lookup_flags =
      (SpectrumMetaDataLookup::MDF_RT |
       SpectrumMetaDataLookup::MDF_PRECURSORMZ |
       SpectrumMetaDataLookup::MDF_PRECURSORCHARGE);

    if (lookup.reference_formats.empty())
    {
      // MS-GF+ Percolator (mzid?) format:
      lookup.addReferenceFormat(R"(_SII_(?<INDEX1>\d+)_\d+_\d+_(?<CHARGE>\d+)_\d+)");
      // Mascot Percolator format (RT may be missing, e.g. for searches via
      // ProteomeDiscoverer):
      lookup.addReferenceFormat(R"(spectrum:[^;]+[(scans:)(scan=)(spectrum=)](?<INDEX0>\d+)[^;]+;rt:(?<RT>\d*(\.\d+)?);mz:(?<MZ>\d+(\.\d+)?);charge:(?<CHARGE>-?\d+))");
      // X! Tandem Percolator format:
      lookup.addReferenceFormat(R"(_(?<INDEX0>\d+)_(?<CHARGE>\d+)_\d+$)");
    }

    vector<String> items;
    CsvFile source(filename, '\t');
    source.getRow(0, items);
    String header = ListUtils::concatenate<String>(items, '\t');
    if (header != 
        "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds")
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  header, "Not a valid header for Percolator (PSM level) output");
    }

    set<String> accessions;
    Size no_charge = 0, no_rt = 0, no_mz = 0; // counters for missing data
    peptides.clear();

    for (Size row = 1; row < source.rowCount(); ++row)
    {
      source.getRow(row, items);

      SpectrumMetaDataLookup::SpectrumMetaData meta_data;
      try
      {
        lookup.getSpectrumMetaData(items[0], meta_data, lookup_flags);
      }
      catch (...)
      {
        String msg = "Error: Could not extract data for spectrum reference '" +
          items[0] + "' from row " + String(row);
        OPENMS_LOG_ERROR << msg << endl;
      }

      PeptideHit hit;
      if (meta_data.precursor_charge != 0)
      {
        hit.setCharge(meta_data.precursor_charge);
      }
      else
      {
        ++no_charge; // maybe TODO: calculate charge from m/z and peptide mass?
      }

      PeptideIdentification peptide;
      peptide.setIdentifier("id");
      if (!std::isnan(meta_data.rt))
      {
        peptide.setRT(meta_data.rt);
      }
      else
      {
        ++no_rt;
      }
      if (!std::isnan(meta_data.precursor_mz))
      {
        peptide.setMZ(meta_data.precursor_mz);
      }
      else
      {
        ++no_mz;
      }

      double score = items[1].toDouble();
      double qvalue = items[2].toDouble();
      double posterrprob = items[3].toDouble();
      hit.setMetaValue("Percolator_score", score);
      hit.setMetaValue("Percolator_qvalue", qvalue);
      hit.setMetaValue("Percolator_PEP", posterrprob);
      switch (output_score)
      {
        case SCORE:
          hit.setScore(score);
          peptide.setScoreType("Percolator_score");
          peptide.setHigherScoreBetter(true);
          break;
        case QVALUE:
          hit.setScore(qvalue);
          peptide.setScoreType("q-value");
          peptide.setHigherScoreBetter(false);
          break;
        case POSTERRPROB:
          hit.setScore(posterrprob);
          peptide.setScoreType("Posterior Error Probability");
          peptide.setHigherScoreBetter(false);
          break;
        case SIZE_OF_SCORETYPE:
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'output_score' must not be 'SIZE_OF_SCORETYPE'!");
      }

      AASequence seq;
      getPeptideSequence_(items[4], seq);
      hit.setSequence(seq);

      for (Size pos = 5; pos < items.size(); ++pos)
      {
        accessions.insert(items[pos]);
        PeptideEvidence evidence;
        evidence.setProteinAccession(items[pos]);
        hit.addPeptideEvidence(evidence);
      }

      peptide.insertHit(hit);
      peptides.push_back(peptide);
    }

    proteins = ProteinIdentification();
    proteins.setIdentifier("id");
    proteins.setDateTime(DateTime::now());
    proteins.setSearchEngine("Percolator");

    for (set<String>::const_iterator it = accessions.begin();
         it != accessions.end(); ++it)
    {
      ProteinHit hit;
      hit.setAccession(*it);
      proteins.insertHit(hit);
    }

    // add info about allowed modifications:
    ModificationDefinitionsSet mod_defs;
    mod_defs.inferFromPeptides(peptides);
    ProteinIdentification::SearchParameters params;
    mod_defs.getModificationNames(params.fixed_modifications,
                                  params.variable_modifications);
    proteins.setSearchParameters(params);

    OPENMS_LOG_INFO << "Created " << proteins.getHits().size() << " protein hits.\n"
             << "Created " << peptides.size() << " peptide hits (PSMs)."
             << endl;
    if (no_charge > 0)
    {
      OPENMS_LOG_WARN << no_charge << " peptide hits without charge state information."
               << endl;
    }
    if (no_rt > 0)
    {
      OPENMS_LOG_WARN << no_rt << " peptide hits without retention time information." 
               << endl;
    }
    if (no_mz > 0)
    {
      OPENMS_LOG_WARN << no_mz << " peptide hits without mass-to-charge information." 
               << endl;
    }
  }

} // namespace OpenMS
