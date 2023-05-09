// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PercolatorInfile.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/FORMAT/CsvFile.h>

#include <regex>
#include <functional>

namespace OpenMS
{
  using namespace std;

  void PercolatorInfile::store(const String& pin_file,
    const vector<PeptideIdentification>& peptide_ids, 
    const StringList& feature_set, 
    const std::string& enz, 
    int min_charge, 
    int max_charge)
  {
    TextFile txt = preparePin_(peptide_ids, feature_set, enz, min_charge, max_charge);
    txt.store(pin_file);
  }

  // uses spectrum_reference, if empty uses spectrum_id, if also empty fall back to using index
  String PercolatorInfile::getScanIdentifier(const PeptideIdentification& pid, size_t index)
  {
    // MSGF+ uses this field, is empty if not specified
    String scan_identifier = pid.getMetaValue("spectrum_reference");

    if (scan_identifier.empty())
    {
      // XTandem uses this (integer) field
      // these ids are 1-based in contrast to the index which is 0-based. This might be problematic to use for merging
      if (pid.metaValueExists("spectrum_id") && !pid.getMetaValue("spectrum_id").toString().empty())
      {
        scan_identifier = "scan=" + pid.getMetaValue("spectrum_id").toString();
      }
      else
      {
        scan_identifier = "index=" + String(index); // fall back
        OPENMS_LOG_WARN << "no known spectrum identifiers, using index [1,n] - use at own risk." << endl;
      }
    }
    return scan_identifier.removeWhitespaces();
  }

  vector<PeptideIdentification> PercolatorInfile::load(const String& pin_file, bool higher_score_better, const String& score_name, String decoy_prefix)
  {
    CsvFile csv(pin_file, '\t');
    StringList header;
    csv.getRow(0, header);

    unordered_map<String, size_t> to_idx; // map column name to column index
    {
      size_t idx{};
      for (const auto& h : header) { to_idx[h] = idx++; }
    }

    // charge columns are not standardized so we check for the format and create hash to lookup column name to charge mapping
    std::regex charge_one_hot_pattern("^charge\\d+$");
    std::regex sage_one_hot_pattern("^z=\\d+$");
    String charge_prefix;
    unordered_map<String, int> col_name_to_charge;
    for (const String& c : header)
    {
      if (std::regex_match(c, charge_one_hot_pattern))
      {
        col_name_to_charge[c] = c.substr(6).toInt();
        charge_prefix = "charge";
      }
      else if (std::regex_match(c, sage_one_hot_pattern))
      {
        col_name_to_charge[c] = c.substr(2).toInt();
        charge_prefix = "z=";
      }
      else if (c == "z=other") // SAGE
      {
        col_name_to_charge[c] = 0;
      }
    }

    auto n_rows = csv.rowCount();
    
    vector<PeptideIdentification> pids;
    pids.reserve(n_rows);
    String spec_id;
    for (size_t i = 1; i != n_rows; ++i)
    {
      StringList row;      
      csv.getRow(i, row);

      if (row.size() != header.size())
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: line " + String(i) + " of file '" + pin_file + "' does not have the same number of columns as the header!", String(i));
      }

      const String& sSpecId = row[to_idx.at("SpecId")];
      if (sSpecId != spec_id)
      {
        pids.resize(pids.size() + 1);
        pids.back().setHigherScoreBetter(higher_score_better);
        pids.back().setScoreType(score_name);
      }

      int sScanNr = row[to_idx.at("ScanNr")].toInt();

      String sPeptide = row[to_idx.at("Peptide")];
      const double score = row[to_idx.at(score_name)].toDouble();
      String target_decoy = row[to_idx.at("Label")].toInt() == 1 ? "target" : "decoy";
      const String& sProteins = row[to_idx.at("Proteins")];
      int rank = to_idx.count("rank") ? row[to_idx.at("rank")].toInt() : 1;
      StringList accessions;

      int charge = 0;
      for (const auto&[name, z] : col_name_to_charge)
      {
        if (row[to_idx.at(name)] == "1")
        {
          charge = z;
          break;
        }
      }

      sProteins.split(';', accessions);

      // deduce decoy state from accessions if decoy_prefix is set
      if (!decoy_prefix.empty())
      {
        target_decoy = std::all_of(accessions.begin(), accessions.end(), [&decoy_prefix](const String& acc) { return acc.hasPrefix(decoy_prefix); }) ? "decoy" : "target" ;
      }          

      // needs to handle strings like: [+42]-MVLVQDLLHPTAASEAR, [+304.207]-ETC[+57.0215]RQLGLGTNIYNAER
      sPeptide.substitute("]-", "]."); // we can parse [+42].MVLVQDLLHPTAASEAR
      AASequence aa_seq = AASequence::fromString(sPeptide);
      PeptideHit ph(score, rank, charge, std::move(aa_seq));
      ph.setMetaValue("SpecId", sSpecId);
      ph.setMetaValue("ScanNr", sScanNr);
      ph.setMetaValue("target_decoy", target_decoy);
      ph.setRank(rank);

      // add link to protein (we only know the accession but not start/end, aa_before/after in protein at this point)
      for (const String& accession : accessions)
      {
        ph.addPeptideEvidence(PeptideEvidence(accession));
      }
      
      pids.back().insertHit(std::move(ph));
    }
    return pids;
  }


  TextFile PercolatorInfile::preparePin_(
    const vector<PeptideIdentification>& peptide_ids, 
    const StringList& feature_set, 
    const std::string& enz, 
    int min_charge, 
    int max_charge)
  {
    TextFile txt;  
    txt.addLine(ListUtils::concatenate(feature_set, '\t'));
    if (peptide_ids.empty()) 
    {
      OPENMS_LOG_WARN << "No identifications provided. Creating empty percolator input." << endl;
      return txt;
    }

    // extract native id (usually in spectrum_reference)
    const String sid = getScanIdentifier(peptide_ids[0], 0);

    // determine RegEx to extract scan/index number
    boost::regex scan_regex = boost::regex(SpectrumLookup::getRegExFromNativeID(sid));

    // keep track of errors
    size_t missing_meta_value_count{};
    set<String> missing_meta_values;

    size_t index = 0;
    for (const PeptideIdentification& pep_id : peptide_ids)
    {
      index++;
      // try to make a file and scan unique identifier
      String scan_identifier = getScanIdentifier(pep_id, index);
      String file_identifier = pep_id.getMetaValue("file_origin", String());

      file_identifier += (String)pep_id.getMetaValue("id_merge_index", String());

      Int scan_number = SpectrumLookup::extractScanNumber(scan_identifier, scan_regex, true);
      
      double exp_mass = pep_id.getMZ();
      double retention_time = pep_id.getRT();
      for (const PeptideHit& psm : pep_id.getHits())
      {
        if (psm.getPeptideEvidences().empty())
        {
          OPENMS_LOG_WARN << "PSM (PeptideHit) without protein reference found. "
                  << "This may indicate incomplete mapping during PeptideIndexing (e.g., wrong enzyme settings)." 
                  << "Will skip this PSM." << endl;
          continue;
        }
        PeptideHit hit(psm); // make a copy of the hit to store temporary features
        hit.setMetaValue("SpecId", file_identifier + scan_identifier);
        hit.setMetaValue("ScanNr", scan_number);
        
        if (!hit.metaValueExists("target_decoy") 
          || hit.getMetaValue("target_decoy").toString().empty()) 
        {
          continue;
        }
        
        int label = 1;
        if (hit.getMetaValue("target_decoy") == "decoy")
        {
          label = -1;
        }
        hit.setMetaValue("Label", label);
        
        int charge = hit.getCharge();
        String unmodified_sequence = hit.getSequence().toUnmodifiedString();
      
        double calc_mass; 
        if (!hit.metaValueExists("CalcMass"))
        {
          calc_mass = hit.getSequence().getMZ(charge);
          hit.setMetaValue("CalcMass", calc_mass); // Percolator calls is CalcMass instead of m/z
        }
        else
        {
          calc_mass = hit.getMetaValue("CalcMass");
        }

        if (hit.metaValueExists("IsotopeError"))  // for backwards compatibility (generated by MSGFPlusAdaper OpenMS < 2.6)
        {
          float isoErr = hit.getMetaValue("IsotopeError").toString().toFloat();
          exp_mass = exp_mass - (isoErr * Constants::C13C12_MASSDIFF_U) / charge;
        }
        else if (hit.metaValueExists(Constants::UserParam::ISOTOPE_ERROR)) // OpenMS user param name for isotope error
        {
          float isoErr = hit.getMetaValue(Constants::UserParam::ISOTOPE_ERROR).toString().toFloat();
          exp_mass = exp_mass - (isoErr * Constants::C13C12_MASSDIFF_U) / charge;
        }
                
        hit.setMetaValue("ExpMass", exp_mass);

        // needed in case "description of correct" option is used
        double delta_mass = exp_mass - calc_mass;
        hit.setMetaValue("deltamass", delta_mass);
        hit.setMetaValue("retentiontime", retention_time);

        hit.setMetaValue("mass", exp_mass);
        
        double score = hit.getScore();
        // TODO better to use log scores for E-value based scores
        hit.setMetaValue("score", score);
        
        int peptide_length = unmodified_sequence.size();
        hit.setMetaValue("peplen", peptide_length);
        
        for (int i = min_charge; i <= max_charge; ++i)
        {
          hit.setMetaValue("charge" + String(i), charge == i);
        }

        // just first peptide evidence
        char aa_before = hit.getPeptideEvidences().front().getAABefore();
        char aa_after = hit.getPeptideEvidences().front().getAAAfter();

        bool enzN = isEnz_(aa_before, unmodified_sequence.prefix(1)[0], enz);
        hit.setMetaValue("enzN", enzN);
        bool enzC = isEnz_(unmodified_sequence.suffix(1)[0], aa_after, enz);
        hit.setMetaValue("enzC", enzC);
        int enzInt = countEnzymatic_(unmodified_sequence, enz);
        hit.setMetaValue("enzInt", enzInt);

        hit.setMetaValue("dm", delta_mass);
        
        double abs_delta_mass = abs(delta_mass);
        hit.setMetaValue("absdm", abs_delta_mass);
        
        //peptide
        String sequence = "";

        aa_before = aa_before == '[' ? '-' : aa_before;
        aa_after = aa_after == ']' ? '-' : aa_after;

        sequence += aa_before;
        sequence += "."; 
        // Percolator uses square brackets to indicate PTMs
        sequence += hit.getSequence().toBracketString(false, true);
        sequence += "."; 
        sequence += aa_after;
        
        hit.setMetaValue("Peptide", sequence);
        
        //proteinId1
        StringList proteins;
        for (const PeptideEvidence& pep : hit.getPeptideEvidences())
        {
          proteins.push_back(pep.getProteinAccession());
        }
        hit.setMetaValue("Proteins", ListUtils::concatenate(proteins, '\t'));
        
        StringList feats;
        for (const String& feat : feature_set)
        {
        // Some Hits have no NumMatchedMainIons, and MeanError, etc. values. Have to ignore them!
          if (hit.metaValueExists(feat))
          {
            feats.push_back(hit.getMetaValue(feat).toString());
          }
        }
        // here: feats (metavalues in peptide hits) and feature_set are equal if they have same size (if no metavalue is missing)

        if (feats.size() == feature_set.size())
        { // only if all feats were present add
          txt.addLine(ListUtils::concatenate(feats, '\t'));
        }        
        else
        { // at least one feature is missing in the current peptide hit
          missing_meta_value_count++;        
          for (const auto& f : feature_set)
          {
            if (std::find(feats.begin(), feats.end(), f) == feats.end()) missing_meta_values.insert(f);
          }
        }
      }
    }

    // print warnings
    if (missing_meta_value_count != 0)
    {
      OPENMS_LOG_WARN << "There were peptide hits with missing features/meta values. Skipped peptide hits: " << missing_meta_value_count << endl;
      OPENMS_LOG_WARN << "Names of missing meta values: " << endl;
      for (const auto& f : missing_meta_values)
      {
        OPENMS_LOG_WARN << f << endl;
      }
    }

    return txt;
  }


  bool PercolatorInfile::isEnz_(const char& n, const char& c, const std::string& enz)
  {
    if (enz == "trypsin")
    {
      return ((n == 'K' || n == 'R') && c != 'P') || n == '-' || c == '-';
    }
    else if (enz == "trypsinp")
    {
      return (n == 'K' || n == 'R') || n == '-' || c == '-';
    }
    else if (enz == "chymotrypsin")
    {
      return ((n == 'F' || n == 'W' || n == 'Y' || n == 'L') && c != 'P') || n == '-' || c == '-';
    }
    else if (enz == "thermolysin")
    {
      return ((c == 'A' || c == 'F' || c == 'I' || c == 'L' || c == 'M'
              || c == 'V' || (n == 'R' && c == 'G')) && n != 'D' && n != 'E') || n == '-' || c == '-';
    }
    else if (enz == "proteinasek")
    {
      return (n == 'A' || n == 'E' || n == 'F' || n == 'I' || n == 'L'
             || n == 'T' || n == 'V' || n == 'W' || n == 'Y') || n == '-' || c == '-';
    }
    else if (enz == "pepsin")
    {
      return ((c == 'F' || c == 'L' || c == 'W' || c == 'Y' || n == 'F'
              || n == 'L' || n == 'W' || n == 'Y') && n != 'R') || n == '-' || c == '-';
    }
    else if (enz == "elastase")
    {
      return ((n == 'L' || n == 'V' || n == 'A' || n == 'G') && c != 'P')
             || n == '-' || c == '-';
    }
    else if (enz == "lys-n")
    {
      return (c == 'K')
             || n == '-' || c == '-';
    }
    else if (enz == "lys-c")
    {
      return ((n == 'K') && c != 'P')
             || n == '-' || c == '-';
    }
    else if (enz == "arg-c")
    {
      return ((n == 'R') && c != 'P')
             || n == '-' || c == '-';
    }
    else if (enz == "asp-n")
    {
      return (c == 'D')
             || n == '-' || c == '-';
    }
    else if (enz == "glu-c")
    {
      return ((n == 'E') && (c != 'P'))
             || n == '-' || c == '-';
    }
    else
    {
      return true;
    }
  }

  // Function adapted from Enzyme.h in Percolator converter
  // TODO: Use existing OpenMS functionality.
  Size PercolatorInfile::countEnzymatic_(const String& peptide, const string& enz)
  {
    Size count = 0;
    for (Size ix = 1; ix < peptide.size(); ++ix)
    {
      if (isEnz_(peptide[ix - 1], peptide[ix], enz))
      {
        ++count;
      }
    }
    return count;
  }

}

