// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Johannes von Kleist $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PercolatorInfile.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <regex>
#include <functional>
#include <unordered_set>

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
    String scan_identifier = pid.getSpectrumReference();

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

  vector<PeptideIdentification> PercolatorInfile::load(
    const String& pin_file,
    bool higher_score_better,
    const String& score_name,
    const StringList& extra_scores,
    StringList& filenames,
    String decoy_prefix, 
    double threshold, 
    bool SageAnnotation)
  {
    CsvFile csv(pin_file, '\t');
    
    //Sage Variables, initialized in the following block if SageAnnotation is set 
    map<int, vector<PeptideHit::PeakAnnotation>> anno_mapping; 
    CsvFile tsv; 
    CsvFile annos; 
    unordered_map<String, size_t> to_idx_t; 

    if (SageAnnotation) // Block for special treatment of sage 
    {
      String tsv_file_path = pin_file.substr(0, pin_file.size()-3);
      tsv_file_path = tsv_file_path + "tsv"; 
      tsv = CsvFile(tsv_file_path,'\t'); 

      String temp_diff = "results.sage.pin"; 
      String anno_file_path = pin_file.substr(0, pin_file.size()-temp_diff.length());
      anno_file_path = anno_file_path + "matched_fragments.sage.tsv"; 
      annos = CsvFile(anno_file_path, '\t'); 
      //map PSMID to vec of PeakAnnotation 
      StringList sage_tsv_header;
      tsv.getRow(0, sage_tsv_header); 
      {
        int idx_t{};
        for (const auto& h : sage_tsv_header) { to_idx_t[h] = idx_t++; }
      }

      // processs annotation file
      StringList sage_annotation_header; 
      annos.getRow(0, sage_annotation_header);
      unordered_map<String, size_t> to_idx_a; // map column name to column index, for full annotation file file 
      {
        int idx_a{};
        for (const auto& h : sage_annotation_header) { to_idx_a[h] = idx_a++; }
      }
      // map PSMs -> PeakAnnotation vector
      auto num_rows = annos.rowCount(); 

      for (size_t i = 1; i < num_rows; ++i)
      {
        StringList row;
        annos.getRow(i, row);
        
        //Check if mapping already has PSM, if it does add 
        if (anno_mapping.find(row[to_idx_a.at("psm_id")].toInt()) == anno_mapping.end())
        {                 
          //Make a new vector of annotations 
          PeptideHit::PeakAnnotation peak_temp; 

          peak_temp.annotation = row[to_idx_a.at("fragment_type")] + row[to_idx_a.at("fragment_ordinals")]; 
          peak_temp.charge = row[to_idx_a.at("fragment_charge")].toInt(); 
          peak_temp.intensity = row[to_idx_a.at("fragment_intensity")].toDouble(); 
          peak_temp.mz = row[to_idx_a.at("fragment_mz_experimental")].toDouble(); 

          vector<PeptideHit::PeakAnnotation> temp_anno_vec; 
          temp_anno_vec.push_back(peak_temp); 
          anno_mapping[ row[to_idx_a.at("psm_id")].toInt() ] = temp_anno_vec;
        }
        else
        {
          //Add values to exisiting vector 
          PeptideHit::PeakAnnotation peak_temp; 

          peak_temp.annotation = row[to_idx_a.at("fragment_type")] + row[to_idx_a.at("fragment_ordinals")]; 
          peak_temp.charge = row[to_idx_a.at("fragment_charge")].toInt(); 
          peak_temp.intensity = row[to_idx_a.at("fragment_intensity")].toDouble(); 
          peak_temp.mz = row[to_idx_a.at("fragment_mz_experimental")].toDouble(); 

          anno_mapping[ row[to_idx_a.at("psm_id")].toInt() ].push_back(peak_temp);         
        }
      }      
    }
    
    StringList pin_header;

    csv.getRow(0, pin_header);

    unordered_map<String, size_t> to_idx; // map column name to column index
    {
      int idx{};
      for (const auto& h : pin_header) { to_idx[h] = idx++; }
    }

    // determine file name column index in percolator in file
    int file_name_column_index{-1};
    if (auto it = std::find(pin_header.begin(), pin_header.end(), "FileName"); it != pin_header.end())
    {
      file_name_column_index = it - pin_header.begin();
    }
 
    // determine extra scores and store column indices
    std::set<String> found_extra_scores; // additional (non-main) scores that should be stored in the PeptideHit, order important for comparable idXML  
    for (const String& s : extra_scores)
    {
      if (auto it = std::find(pin_header.begin(), pin_header.end(), s); it != pin_header.end())
      {
        found_extra_scores.insert(s);
      }
      else
      {
        OPENMS_LOG_WARN << "Extra score: " << s << " not found in Percolator input file." << endl;
      }
    }
      
    // charge columns are not standardized, so we check for the format and create hash to lookup column name to charge mapping
    std::regex charge_one_hot_pattern("^charge\\d+$");
    std::regex sage_one_hot_pattern("^z=\\d+$");
    String charge_prefix;
    unordered_map<String, int> col_name_to_charge;

    // Special handling for sage: sage produces an additional column for PSMs outside of the "suggested" charge search range (e.g., charge 2-5).
    // The reason is that sage searches always for the charge annotated in the spectrum raw file. Only if the annotation is missing it will search
    // the suggested charge range.
    bool found_sage_otherz_charge_column{false}; 
    for (const String& c : pin_header)
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
        found_sage_otherz_charge_column = true;
        OPENMS_LOG_DEBUG << "Found SAGE charge column 'z=other'. Will extract charge from this column if charge was not set in the one-hot encoded charge columns." << endl;
      }
    }

    auto n_rows = csv.rowCount();

    vector<PeptideIdentification> pids;
    pids.reserve(n_rows);
    String spec_id;
    String raw_file_name("UNKNOWN");
    unordered_map<String, size_t> map_filename_to_idx; // fast lookup of filename to index in filenames vector

    for (size_t i = 1; i != n_rows; ++i)
    {
      StringList row;
      csv.getRow(i, row);

      StringList t_row; 

      if (SageAnnotation)
      {
        tsv.getRow(i, t_row); 
        // skip if spectrum_q is above threshold
        if (t_row[to_idx_t.at("spectrum_q")].toDouble() > threshold ) continue;
      }

      if (row.size() != pin_header.size())
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: line " + String(i) + " of file '" + pin_file + "' does not have the same number of columns as the pin_header!", String(i));
      }


      if (file_name_column_index >= 0)
      {
        raw_file_name = row[file_name_column_index];
        if (map_filename_to_idx.find(raw_file_name) == map_filename_to_idx.end())
        {
          filenames.push_back(raw_file_name);
          map_filename_to_idx[raw_file_name] = filenames.size() - 1;
        }
      }
      // NOTE: In our pin files that we WRITE, SpecID will be filename + vendor spectrum native ID
      // However, many search engines (e.g. Sage) choose arbitrary IDs, which is unfortunately allowed
      //  by this loosely defined format.
      const String& sSpecId = row[to_idx.at("SpecId")];

      if (auto it = to_idx.find("ion_mobility"); it != to_idx.end())
      {
        const String& sIM = row[it->second];
        const double IM = sIM.toDouble();   
        if (!pids.empty())  pids.back().setMetaValue(Constants::UserParam::IM, IM);
      }
      // In theory, this should be an integer, but Sage currently cannot extract the number from all vendor spectrum IDs,
      //  so it writes the full ID as string
      String sScanNr = row[to_idx.at("ScanNr")];
      if (sSpecId != spec_id)
      {
        pids.resize(pids.size() + 1);
        pids.back().setHigherScoreBetter(higher_score_better);
        pids.back().setScoreType(score_name);
        pids.back().setMetaValue(Constants::UserParam::ID_MERGE_INDEX, map_filename_to_idx.at(raw_file_name));
        pids.back().setRT(row[to_idx.at("retentiontime")].toDouble() * 60.0); // search engines typically write minutes (e.g., sage)
        pids.back().setMetaValue("PinSpecId", sSpecId);  
        // Since ScanNr is the closest to help in identifying the spectrum in the file later on,
        // we use it as spectrum_reference. Since it can be integer only or the complete
        // vendor ID, you will need a lookup in case of number only later!!
        pids.back().setSpectrumReference(sScanNr);     
      }
      String sPeptide = row[to_idx.at("Peptide")];
      const double score = row[to_idx.at(score_name)].toDouble();
      String target_decoy = row[to_idx.at("Label")].toInt() == 1 ? "target" : "decoy";
      const String& sProteins = row[to_idx.at("Proteins")];
      int rank = to_idx.count("rank") ? row[to_idx.at("rank")].toInt() : 1;
      StringList accessions;

      int charge = 0;
      for (const auto& [name, z] : col_name_to_charge)
      {        
        if (row[to_idx.at(name)] == "1")
        {
          charge = z;
          break;
        }
      }

      // all one-hot encoded charge columns are zero. Use value in the sage "z=other" column if it exists.
      if (charge == 0 && found_sage_otherz_charge_column)
      {
        charge = row[to_idx.at("z=other")].toInt();
      }
      
      if (charge != 0)
      {
        pids.back().setMZ(row[to_idx.at("ExpMass")].toDouble() / std::fabs(charge) + Constants::PROTON_MASS_U);
      }

      sProteins.split(';', accessions);

      // deduce decoy state from accessions if decoy_prefix is set
      if (!decoy_prefix.empty())
      {
        bool dec = false;
        bool tgt = false;
        for (const auto& acc : accessions)
        {
          if (!(dec && tgt))
          {
            if (acc.hasPrefix(decoy_prefix))
            {
              dec = true;
            }
            else
            {
              tgt = true;
            }
          }
          else
          {
            break;
          }
        }
        if (tgt && dec)
        {
          target_decoy = "target+decoy";
        }
        else if (tgt)
        {
          target_decoy = "target";
        }
        else {
          target_decoy = "decoy";
        }
      }

      // needs to handle strings like: [+42]-MVLVQDLLHPTAASEAR, [+304.207]-ETC[+57.0215]RQLGLGTNIYNAER etc.
      sPeptide.substitute("]-", "]."); // we can parse [+42].MVLVQDLLHPTAASEAR
      sPeptide.substitute("-[", ".["); // we can parse MVLVQDLLHPTAASEAR.[+111]
      AASequence aa_seq = AASequence::fromString(sPeptide);
      PeptideHit ph(score, rank, charge, std::move(aa_seq));
      ph.setMetaValue("target_decoy", target_decoy);

      for (const auto& name : found_extra_scores)
      {
        ph.setMetaValue(name, row[to_idx.at(name)]);
      }
      ph.setRank(rank);

      // adding own meta values 
      if (SageAnnotation)
      {
        ph.setMetaValue("spectrum_q", t_row[to_idx_t.at("spectrum_q")].toDouble());  //TODO: check if column exists / SAGE specific treatment
      }
      ph.setMetaValue("DeltaMass", ( row[to_idx.at("ExpMass")].toDouble() - row[to_idx.at("CalcMass")].toDouble()) ); 
      // add annotations
      if (SageAnnotation)
      {
        if (anno_mapping.find(sSpecId.toInt()) != anno_mapping.end())
        {
          // copy annotations from mapping to PeptideHit
          vector<PeptideHit::PeakAnnotation> pep_vec;
          for (const PeptideHit::PeakAnnotation& pep : anno_mapping[sSpecId.toInt()])
          {
            pep_vec.push_back(pep) ; 
          }        
          ph.setPeakAnnotations(pep_vec);        
        } 
      }
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