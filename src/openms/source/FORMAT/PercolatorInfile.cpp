// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Johannes von Kleist $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PercolatorInfile.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>

#include <regex>
#include <functional>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
  using namespace std;

  void PercolatorInfile::store(const String& pin_file,
    const PeptideIdentificationList& peptide_ids, 
    const StringList& feature_set, 
    const std::string& enz, 
    int min_charge, int max_charge)
  {
    TextFile txt = preparePin_(peptide_ids, feature_set, enz, min_charge, max_charge);
    txt.store(pin_file);
  }

  PeptideIdentificationList PercolatorInfile::load(const String& pin_file)
  {
    PeptideIdentificationList peptide_ids;
    CsvFile csv(pin_file, '\t');
    
    // map ScanId (from percolator file) to PeptideHit
    map<String, PeptideHit> scan_id_to_peptide_hit;

    // determine column indices
    StringList header;
    csv.getRow(0, header);
    Size scan_id_col = -1;
    Size score_col = -1;
    Size q_value_col = -1;
    Size pep_col = -1;
    Size peptide_col = -1;
    Size protein_col = -1;
    
    for (Size i = 0; i < header.size(); ++i)
    {
      if (header[i] == "ScanId") scan_id_col = i;
      else if (header[i] == "score") score_col = i;
      else if (header[i] == "q-value") q_value_col = i;
      else if (header[i] == "posterior_error_prob") pep_col = i;
      else if (header[i] == "peptide") peptide_col = i;
      else if (header[i] == "proteinIds") protein_col = i;
    }

    if (scan_id_col == -1 || score_col == -1 || q_value_col == -1 || pep_col == -1 || peptide_col == -1 || protein_col == -1)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Percolator input file is missing required columns (ScanId, score, q-value, posterior_error_prob, peptide, proteinIds).", "");
    }

    // parse rows
    StringList row;
    for (Size i = 1; i < csv.rowCount(); ++i)
    {
      csv.getRow(i, row);
      PeptideHit hit;
      hit.setScore(row[score_col].toDouble());
      hit.setMetaValue("spectral_q_value", row[q_value_col].toDouble()); // using generic name for q-value
      hit.setMetaValue("Posterior Error Probability", row[pep_col].toDouble()); // using generic name for PEP
      
      // parse peptide sequence (format: X.SEQUENCE.X)
      String seq = row[peptide_col];
      if (seq.size() > 4 && seq[1] == '.' && seq[seq.size() - 2] == '.')
      {
        seq = seq.substr(2, seq.size() - 4);
      }
      hit.setSequence(AASequence::fromString(seq));

      // parse protein accessions
      String proteins = row[protein_col];
      // remove leading/trailing whitespace
      proteins.trim();
      vector<String> protein_accessions;
      proteins.split('\t', protein_accessions); // not sure about delimiter here?
      
      for (const auto& acc : protein_accessions)
      {
        PeptideEvidence pe;
        pe.setProteinAccession(acc);
        hit.addPeptideEvidence(pe);
      }

      // store using ScanId as key to reconstruct structure later if needed
      // but for now we just return a flat list
      peptide_ids.insertHit(hit);
    }

    return peptide_ids;
  }

  String PercolatorInfile::getScanIdentifier(const String& path_to_spectrum_file, const String& spectrum_native_id)
  {
    // check if we have a special case where we want to use the native ID as scan identifier
    // this is the case for MzML files where we can just use the native ID
    // for MzXML etc we might need to parse the scan number
    return File::basename(path_to_spectrum_file) + ":" + spectrum_native_id;
  }

  TextFile PercolatorInfile::preparePin_(const PeptideIdentificationList& peptide_ids,
                                       const StringList& feature_set,
                                       const std::string& enz,
                                       int min_charge,
                                       int max_charge)
  {
    ProteaseDigestion digest;
    if (!enz.empty() && enz != "no_enzyme")
    {
      try
      {
        digest.setEnzyme(enz);
      }
      catch (Exception::ElementNotFound&)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown enzyme name passed to PercolatorInfile", enz);
      }
    }
    else if (enz == "no_enzyme")
    {
      digest.setEnzyme("unspecific cleavage");
    }

    TextFile txt;
    txt.addLine(ListUtils::concatenate(feature_set, '\t'));
    if (peptide_ids.empty())
    {
      OPENMS_LOG_WARN << "No identifications provided. Creating empty percolator input." << std::endl;
      return txt;
    }

    // extract native id (usually in spectrum_reference)
    const String sid = getScanIdentifier(peptide_ids[0], 0);

    // iterate over identifications
    for (Size i = 0; i < peptide_ids.getHits().size(); ++i)
    {
      const PeptideHit& hit = peptide_ids.getHits()[i];
      int label = 0;
      if (hit.getMetaValue("target_decoy").toString() == "target") label = 1;
      else if (hit.getMetaValue("target_decoy").toString() == "decoy") label = -1;
      else
      {
        OPENMS_LOG_WARN << "PSM without target/decoy information found. "
                        << "This may indicate incomplete mapping during PeptideIndexing (e.g., wrong decoy prefix settings). " 
                        << "Will skip this PSM." << std::endl;
        continue;
      }

      if (hit.getPeptideEvidences().empty())
      {
        OPENMS_LOG_WARN << "PSM (PeptideHit) without protein reference found. "
                        << "This may indicate incomplete mapping during PeptideIndexing (e.g., wrong enzyme settings). " 
                        << "Will skip this PSM." << std::endl;
        continue;
      }

      String unmodified_sequence = hit.getSequence().toUnmodifiedString();

      String line = sid + "\t" + String(label) + "\t" + sid + "\t";

      double exp_mass = 0;
      if (hit.metaValueExists("mass"))
      {
        exp_mass = hit.getMetaValue("mass");
      }
      else
      {
        // TODO: calculate from MZ and charge?
      }

      double calc_mass = hit.getSequence().getMonoWeight();
      double delta_mass = exp_mass - calc_mass;

      int charge = hit.getCharge();

      double rt = 0;
      if (hit.metaValueExists("retentiontime"))
      {
        rt = hit.getMetaValue("retentiontime");
      }

      // Enzymatic features
      // just first peptide evidence
      char aa_before = hit.getPeptideEvidences().front().getAABefore();
      char aa_after = hit.getPeptideEvidences().front().getAAAfter();

      bool enzN = isEnz_(aa_before, unmodified_sequence.prefix(1)[0], digest);
      bool enzC = isEnz_(unmodified_sequence.suffix(1)[0], aa_after, digest);
      int enzInt = countEnzymatic_(unmodified_sequence, digest);

      // build feature line
      StringList features;
      for (const auto& f : feature_set)
      {
        if (f == "ExpMass") features.push_back(String(exp_mass));
        else if (f == "CalcMass") features.push_back(String(calc_mass));
        else if (f == "lnrSp") features.push_back(String(hit.getMetaValue("lnrSp", 0.0)));
        else if (f == "deltLCn") features.push_back(String(hit.getMetaValue("deltLCn", 0.0)));
        else if (f == "deltCn") features.push_back(String(hit.getMetaValue("deltCn", 0.0)));
        else if (f == "Xcorr") features.push_back(String(hit.getMetaValue("Xcorr", 0.0)));
        else if (f == "Sp") features.push_back(String(hit.getMetaValue("Sp", 0.0)));
        else if (f == "IonFrac") features.push_back(String(hit.getMetaValue("IonFrac", 0.0)));
        else if (f == "Mass") features.push_back(String(exp_mass));
        else if (f == "PepLen") features.push_back(String(unmodified_sequence.size()));
        else if (f == "Charge1") features.push_back(String(charge == 1 ? 1 : 0));
        else if (f == "Charge2") features.push_back(String(charge == 2 ? 1 : 0));
        else if (f == "Charge3") features.push_back(String(charge == 3 ? 1 : 0));
        else if (f == "Charge4") features.push_back(String(charge == 4 ? 1 : 0));
        else if (f == "EnzN") features.push_back(String(enzN ? 1 : 0));
        else if (f == "EnzC") features.push_back(String(enzC ? 1 : 0));
        else if (f == "EnzInt") features.push_back(String(enzInt));
        else if (f == "LnNumSP") features.push_back(String(hit.getMetaValue("LnNumSP", 0.0)));
        else if (f == "DM") features.push_back(String(delta_mass));
        else if (f == "AbsDM") features.push_back(String(std::abs(delta_mass)));
        else
        {
          features.push_back(String(hit.getMetaValue(f, 0.0)));
        }
      }

      line += ListUtils::concatenate(features, '\t');

      // Peptide and proteins
      line += "\t-." + unmodified_sequence + ".-"; // Percolator peptide format
      
      String proteins;
      for (const auto& pe : hit.getPeptideEvidences())
      {
        if (!proteins.empty()) proteins += "\t";
        proteins += pe.getProteinAccession();
      }
      line += "\t" + proteins;

      txt.addLine(line);
    }

    return txt;
  }

  bool PercolatorInfile::isEnz_(const char& n, const char& c, const ProteaseDigestion& digest)
  {
    // Terminal positions (protein N/C-terminus) are always considered enzymatic
    if (n == '-' || c == '-' || n == '[' || c == ']')
    {
      return true;
    }

    // Construct a minimal 2-residue "protein" and check if the
    // single-residue peptide at position 0 is a valid digestion product
    // (i.e., there is a valid cleavage site between n and c)
    const String mini_protein = String(1, n) + String(1, c);
    return digest.isValidProduct(mini_protein, 0, 1, true);
  }

  Size PercolatorInfile::countEnzymatic_(const String& peptide, const ProteaseDigestion& digest)
  {
    return digest.countInternalCleavageSites(peptide);
  }

} // namespace OpenMS
