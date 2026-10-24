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
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <optional>
#include <tuple>
#include <type_traits>
#include <regex>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
  using namespace std;

  void PercolatorInfile::store(const std::string& pin_file,
    const PeptideIdentificationList& peptide_ids,
    const StringList& feature_set,
    const std::string& enz,
    int min_charge,
    int max_charge)
  {
    TextFile txt = preparePin_(peptide_ids, feature_set, enz, min_charge, max_charge);
    txt.store(pin_file);
  }

  StringList PercolatorInfile::getStandardFeatureSet(int min_charge, int max_charge)
  {
    StringList fs = {
      "SpecId", "Label", "ScanNr",
      "ExpMass", "CalcMass", "mass", "peplen"
    };
    for (int c = min_charge; c <= max_charge; ++c)
    {
      fs.push_back("charge" + StringUtils::toStr(c));
    }
    fs.insert(fs.end(), {"enzN", "enzC", "enzInt", "dm", "absdm"});
    return fs;
  }

  // uses spectrum_reference, if empty uses spectrum_id, if also empty fall back to using index
  std::string PercolatorInfile::getScanIdentifier(
    const PeptideIdentification& pid, size_t index)
  {
    // MSGF+ uses this field, is empty if not specified
    std::string scan_identifier = pid.getSpectrumReference();

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
        scan_identifier = "index=" + StringUtils::toStr(index); // fall back
        OPENMS_LOG_WARN << "no known spectrum identifiers, using index [1,n] - use at own risk." << endl;
      }
    }
    return StringUtils::removeWhitespaces(scan_identifier);
  }

  std::string PercolatorInfile::getFileIdentifier(const PeptideIdentification& pid)
  {
    std::string identifier = pid.getMetaValue("file_origin", std::string()).toString();
    // the index is a number: after the last '|', it cannot run into the file origin ("run1" vs. "run" + 1)
    if (pid.metaValueExists("id_merge_index")) identifier += "|" + pid.getMetaValue("id_merge_index").toString();
    return identifier;
  }

  PeptideIdentificationList PercolatorInfile::load(
    const std::string& pin_file,
    bool higher_score_better,
    const std::string& score_name,
    const StringList& extra_scores,
    StringList& filenames,
    std::string decoy_prefix, 
    double threshold, 
    bool SageAnnotation)
  {
    CsvFile csv(pin_file, '\t');
    
    //Sage Variables, initialized in the following block if SageAnnotation is set 
    map<int, vector<PeptideHit::PeakAnnotation>> anno_mapping; 
    CsvFile tsv; 
    CsvFile annos; 
    unordered_map<std::string, size_t> to_idx_t; 

    if (SageAnnotation) // Block for special treatment of sage 
    {
      std::string tsv_file_path = StringUtils::substr(pin_file, 0, pin_file.size()-3);
      tsv_file_path = tsv_file_path + "tsv"; 
      tsv = CsvFile(tsv_file_path,'\t'); 

      std::string temp_diff = "results.sage.pin"; 
      std::string anno_file_path = StringUtils::substr(pin_file, 0, pin_file.size()-temp_diff.length());
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
      unordered_map<std::string, size_t> to_idx_a; // map column name to column index, for full annotation file file 
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
        if (!anno_mapping.contains(StringUtils::toInt32(row[to_idx_a.at("psm_id")])))
        {                 
          //Make a new vector of annotations 
          PeptideHit::PeakAnnotation peak_temp; 

          peak_temp.annotation = row[to_idx_a.at("fragment_type")] + row[to_idx_a.at("fragment_ordinals")]; 
          peak_temp.charge = StringUtils::toInt32(row[to_idx_a.at("fragment_charge")]); 
          peak_temp.intensity = StringUtils::toDouble(row[to_idx_a.at("fragment_intensity")]); 
          peak_temp.mz = StringUtils::toDouble(row[to_idx_a.at("fragment_mz_experimental")]); 

          vector<PeptideHit::PeakAnnotation> temp_anno_vec; 
          temp_anno_vec.push_back(peak_temp); 
          anno_mapping[ StringUtils::toInt32(row[to_idx_a.at("psm_id")]) ] = temp_anno_vec;
        }
        else
        {
          //Add values to exisiting vector 
          PeptideHit::PeakAnnotation peak_temp; 

          peak_temp.annotation = row[to_idx_a.at("fragment_type")] + row[to_idx_a.at("fragment_ordinals")]; 
          peak_temp.charge = StringUtils::toInt32(row[to_idx_a.at("fragment_charge")]); 
          peak_temp.intensity = StringUtils::toDouble(row[to_idx_a.at("fragment_intensity")]); 
          peak_temp.mz = StringUtils::toDouble(row[to_idx_a.at("fragment_mz_experimental")]); 

          anno_mapping[ StringUtils::toInt32(row[to_idx_a.at("psm_id")]) ].push_back(peak_temp);         
        }
      }      
    }
    
    StringList pin_header;

    csv.getRow(0, pin_header);

    unordered_map<std::string, size_t> to_idx; // map column name to column index
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
    std::set<std::string> found_extra_scores; // additional (non-main) scores that should be stored in the PeptideHit, order important for comparable idXML  
    for (const std::string& s : extra_scores)
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
    std::string charge_prefix;
    unordered_map<std::string, int> col_name_to_charge;

    // Special handling for sage: sage produces an additional column for PSMs outside of the "suggested" charge search range (e.g., charge 2-5).
    // The reason is that sage searches always for the charge annotated in the spectrum raw file. Only if the annotation is missing it will search
    // the suggested charge range.
    bool found_sage_otherz_charge_column{false}; 
    for (const std::string& c : pin_header)
    {
      if (std::regex_match(c, charge_one_hot_pattern))
      {
        col_name_to_charge[c] = StringUtils::toInt32(StringUtils::substr(c, 6));
        charge_prefix = "charge";
      }
      else if (std::regex_match(c, sage_one_hot_pattern))
      {
        col_name_to_charge[c] = StringUtils::toInt32(StringUtils::substr(c, 2));
        charge_prefix = "z=";
      }      
      else if (c == "z=other") // SAGE
      {
        found_sage_otherz_charge_column = true;
        OPENMS_LOG_DEBUG << "Found SAGE charge column 'z=other'. Will extract charge from this column if charge was not set in the one-hot encoded charge columns." << endl;
      }
    }

    auto n_rows = csv.rowCount();

    PeptideIdentificationList pids;
    pids.reserve(n_rows);
    std::string spec_id;
    std::string raw_file_name("UNKNOWN");
    unordered_map<std::string, size_t> map_filename_to_idx; // fast lookup of filename to index in filenames vector

    for (size_t i = 1; i != n_rows; ++i)
    {
      StringList row;
      csv.getRow(i, row);

      StringList t_row; 

      if (SageAnnotation)
      {
        tsv.getRow(i, t_row); 
        // skip if spectrum_q is above threshold
        if (StringUtils::toDouble(t_row[to_idx_t.at("spectrum_q")]) > threshold ) continue;
      }

      if (row.size() != pin_header.size())
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: line " + StringUtils::toStr(i) + " of file '" + pin_file + "' does not have the same number of columns as the pin_header!",StringUtils::toStr(i));
      }


      if (file_name_column_index >= 0)
      {
        raw_file_name = row[file_name_column_index];
        if (!map_filename_to_idx.contains(raw_file_name))
        {
          filenames.push_back(raw_file_name);
          map_filename_to_idx[raw_file_name] = filenames.size() - 1;
        }
      }
      // NOTE: In our pin files that we WRITE, SpecID will be filename + vendor spectrum native ID
      // However, many search engines (e.g. Sage) choose arbitrary IDs, which is unfortunately allowed
      //  by this loosely defined format.
      const std::string& sSpecId = row[to_idx.at("SpecId")];

      double IM = 0.0;
      if (auto it = to_idx.find("ion_mobility"); it != to_idx.end())
      {
        const std::string& sIM = row[it->second]; // read IM from e.g. Sage output
        IM = StringUtils::toDouble(sIM);  
      }
      // In theory, this should be an integer, but Sage currently cannot extract the number from all vendor spectrum IDs,
      //  so it writes the full ID as string
      std::string sScanNr = row[to_idx.at("ScanNr")];
      if (sSpecId != spec_id)
      {
        pids.resize(pids.size() + 1);
        pids.back().setHigherScoreBetter(higher_score_better);
        pids.back().setScoreType(score_name);
        pids.back().setMetaValue(Constants::UserParam::ID_MERGE_INDEX, map_filename_to_idx.at(raw_file_name));
        pids.back().setRT(StringUtils::toDouble(row[to_idx.at("retentiontime")]) * 60.0); // search engines typically write minutes (e.g., sage)
        pids.back().setMetaValue("PinSpecId", sSpecId);
        if (IM > 0.0) // Sage might annotate 0.0 if no IM is present
        {
          pids.back().setMetaValue(Constants::UserParam::IM, IM);
        }
        // Since ScanNr is the closest to help in identifying the spectrum in the file later on,
        // we use it as spectrum_reference. Since it can be integer only or the complete
        // vendor ID, you will need a lookup in case of number only later!!
        pids.back().setSpectrumReference(sScanNr);     
      }
      std::string sPeptide = row[to_idx.at("Peptide")];
      const double score = StringUtils::toDouble(row[to_idx.at(score_name)]);
      std::string target_decoy = StringUtils::toInt32(row[to_idx.at("Label")]) == 1 ? "target" : "decoy";
      const std::string& sProteins = row[to_idx.at("Proteins")];
      int rank = to_idx.count("rank") ? StringUtils::toInt32(row[to_idx.at("rank")]) : 1;
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
        charge = StringUtils::toInt32(row[to_idx.at("z=other")]);
      }
      
      if (charge != 0)
      {
        pids.back().setMZ(StringUtils::toDouble(row[to_idx.at("ExpMass")]) / std::fabs(charge) + Constants::PROTON_MASS_U);
      }

      StringUtils::split(sProteins, ';', accessions);

      // deduce decoy state from accessions if decoy_prefix is set
      if (!decoy_prefix.empty())
      {
        bool dec = false;
        bool tgt = false;
        for (const auto& acc : accessions)
        {
          if (!(dec && tgt))
          {
            if (StringUtils::hasPrefix(acc, decoy_prefix))
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
      StringUtils::substitute(sPeptide, "]-", "]."); // we can parse [+42].MVLVQDLLHPTAASEAR
      StringUtils::substitute(sPeptide, "-[", ".["); // we can parse MVLVQDLLHPTAASEAR.[+111]
      AASequence aa_seq = AASequence::fromString(sPeptide);
      PeptideHit ph(score, rank - 1, charge, std::move(aa_seq));
      ph.setMetaValue("target_decoy", target_decoy);

      for (const auto& name : found_extra_scores)
      {
        std::string value = row[to_idx.at(name)];
        if (name == "ln(-poisson)" && value == "inf") value = "3.5"; // workaround for Sage
        ph.setMetaValue(name, value);
      }

      // adding own meta values 
      if (SageAnnotation)
      {
        ph.setMetaValue("spectrum_q", StringUtils::toDouble(t_row[to_idx_t.at("spectrum_q")]));  //TODO: check if column exists / SAGE specific treatment
      }
      ph.setMetaValue("DeltaMass", ( StringUtils::toDouble(row[to_idx.at("ExpMass")]) - StringUtils::toDouble(row[to_idx.at("CalcMass")])) ); 
      // add annotations
      if (SageAnnotation)
      {
        if (anno_mapping.contains(StringUtils::toInt32(sSpecId)))
        {
          // copy annotations from mapping to PeptideHit
          vector<PeptideHit::PeakAnnotation> pep_vec;
          for (const PeptideHit::PeakAnnotation& pep : anno_mapping[StringUtils::toInt32(sSpecId)])
          {
            pep_vec.push_back(pep) ; 
          }        
          ph.setPeakAnnotations(pep_vec);        
        } 
      }
      // add link to protein (we only know the accession but not start/end, aa_before/after in protein at this point)
      for (const std::string& accession : accessions)
      {
        ph.addPeptideEvidence(PeptideEvidence(accession));
      }

      pids.back().insertHit(std::move(ph));
    }

    return pids;
  }


  template <bool with_strings, typename Exists, typename Get, typename Put>
  double PercolatorInfile::pinHitColumns_(const PeptideHit& hit, bool decoy, const std::string& spec_id, Int scan_number,
    double exp_mass, double retention_time, const std::string& enz, int min_charge, int max_charge,
    Exists&& exists, Get&& get, Put&& put)
  {
    if constexpr (with_strings) put(PIN_SPEC_ID, spec_id);
    put(PIN_SCAN_NR, scan_number);
    put(PIN_LABEL, decoy ? -1 : 1);

    const int charge = hit.getCharge();
    const std::string unmodified_sequence = hit.getSequence().toUnmodifiedString();

    double calc_mass;
    if (!exists(LOOKUP_CALC_MASS))
    {
      calc_mass = hit.getSequence().getMZ(charge);
      put(PIN_CALC_MASS, calc_mass);
    }
    else
    {
      calc_mass = getFeatureValue(get(LOOKUP_CALC_MASS), "CalcMass");
    }

    // The precursor isotope error (observed minus theoretical, in 13C spacings) is removed from the mass
    // difference only. ExpMass and mass stay the observed precursor m/z, the same for every PSM of the spectrum:
    // Percolator takes PSMs with the same ScanNr but different ExpMass for different spectra, so PSMs at
    // different isotope errors would not compete with each other.
    double isotope_corrected_mass = exp_mass;
    if (exists(LOOKUP_ISOTOPE_ERROR_LEGACY))  // legacy MSGF+ adapter meta (<2.6)
    {
      const float isoErr = StringUtils::toFloat(get(LOOKUP_ISOTOPE_ERROR_LEGACY).toString());
      isotope_corrected_mass -= (isoErr * Constants::C13C12_MASSDIFF_U) / charge;
    }
    else if (exists(LOOKUP_ISOTOPE_ERROR))
    {
      const float isoErr = StringUtils::toFloat(get(LOOKUP_ISOTOPE_ERROR).toString());
      isotope_corrected_mass -= (isoErr * Constants::C13C12_MASSDIFF_U) / charge;
    }

    put(PIN_EXP_MASS, exp_mass);

    const double delta_mass = isotope_corrected_mass - calc_mass;
    put(PIN_DELTAMASS, delta_mass);
    put(PIN_RETENTIONTIME, retention_time);
    put(PIN_MASS, exp_mass);

    put(PIN_SCORE, hit.getScore());
    put(PIN_PEPLEN, static_cast<int>(unmodified_sequence.size()));

    for (int i = min_charge; i <= max_charge; ++i)
    {
      put(PIN_CHARGE + (i - min_charge), charge == i);
    }

    // use first peptide evidence for flanking AAs
    char aa_before = hit.getPeptideEvidences().front().getAABefore();
    char aa_after  = hit.getPeptideEvidences().front().getAAAfter();

    const bool enzN = isEnz_(aa_before, StringUtils::prefix(unmodified_sequence, 1)[0], enz);
    put(PIN_ENZ_N, enzN);
    const bool enzC = isEnz_(StringUtils::suffix(unmodified_sequence, 1)[0], aa_after, enz);
    put(PIN_ENZ_C, enzC);
    put(PIN_ENZ_INT, static_cast<int>(countEnzymatic_(unmodified_sequence, enz)));

    put(PIN_DM, delta_mass);
    put(PIN_ABSDM, std::abs(delta_mass));

    if constexpr (with_strings)
    {
      // Percolator-formatted peptide string with PTM brackets
      aa_before = aa_before == '[' ? '-' : aa_before;
      aa_after  = aa_after  == ']' ? '-' : aa_after;
      std::string sequence;
      sequence += aa_before;
      sequence += ".";
      sequence += hit.getSequence().toBracketString(false, true);
      sequence += ".";
      sequence += aa_after;
      put(PIN_PEPTIDE, sequence);

      StringList proteins;
      for (const PeptideEvidence& pep : hit.getPeptideEvidences())
      {
        proteins.push_back(pep.getProteinAccession());
      }
      put(PIN_PROTEINS, ListUtils::concatenate(proteins, "\t"));
    }
    return calc_mass;
  }

  namespace
  {
    // The names of the PIN columns (PercolatorInfile::PinColumn_ order; the charge{N} columns from min_charge on)
    std::vector<std::string> pinColumnNames(int min_charge, int max_charge)
    {
      std::vector<std::string> names{"SpecId", "ScanNr", "Label", "CalcMass", "ExpMass", "deltamass", "retentiontime",
                                     "mass", "score", "peplen", "enzN", "enzC", "enzInt", "dm", "absdm", "Peptide",
                                     "Proteins"};
      for (int i = min_charge; i <= max_charge; ++i)
      {
        names.push_back("charge" + StringUtils::toStr(i));
      }
      return names;
    }

    // The names of the PIN columns in the order stamping registers them on its first hit
    std::vector<std::string> pinColumnNamesInStampingOrder(int min_charge, int max_charge)
    {
      std::vector<std::string> names{"SpecId", "ScanNr", "Label", "CalcMass", "ExpMass", "deltamass", "retentiontime",
                                     "mass", "score", "peplen"};
      for (int i = min_charge; i <= max_charge; ++i)
      {
        names.push_back("charge" + StringUtils::toStr(i));
      }
      names.insert(names.end(), {"enzN", "enzC", "enzInt", "dm", "absdm", "Peptide", "Proteins"});
      return names;
    }

    // PeptideHit::getTargetDecoyType() with the meta value name looked up once (index): the lookup by name takes a
    // lock on the registry, which threads would contend for. An unregistered name (index -1) exists on no hit.
    PeptideHit::TargetDecoyType targetDecoyType(const PeptideHit& hit, UInt td_index)
    {
      if (td_index == UInt(-1) || !hit.metaValueExists(td_index))
      {
        return PeptideHit::TargetDecoyType::UNKNOWN;
      }
      std::string td = StringUtils::toLowered(hit.getMetaValue(td_index).toString());
      if (td == "decoy") return PeptideHit::TargetDecoyType::DECOY;
      if (td == "target+decoy") return PeptideHit::TargetDecoyType::TARGET_DECOY;
      if (td == "target") return PeptideHit::TargetDecoyType::TARGET;

      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown value of meta value 'target_decoy'", td);
    }
  }

  std::set<std::pair<size_t, size_t>> PercolatorInfile::stampPinFeaturesOnHits(
    PeptideIdentificationList& peptide_ids,
    const std::string& enz,
    int min_charge,
    int max_charge)
  {
    PinFeatureMetaValueMap added_meta_values;
    return stampPinFeaturesOnHits(peptide_ids, enz, min_charge, max_charge, added_meta_values);
  }

  std::set<std::pair<size_t, size_t>> PercolatorInfile::stampPinFeaturesOnHits(
    PeptideIdentificationList& peptide_ids,
    const std::string& enz,
    int min_charge,
    int max_charge,
    PinFeatureMetaValueMap& added_meta_values)
  {
    std::set<std::pair<size_t, size_t>> skipped;
    added_meta_values.clear();
    if (peptide_ids.empty()) return skipped;

    // Scan-number regex is derived from the first pid's scan identifier.
    const std::string sid = getScanIdentifier(peptide_ids[0], 0);
    RegularExpression scan_regex = RegularExpression(SpectrumLookup::getRegExFromNativeID(sid));

    // Per spectrum (PeptideIdentification): the meta values it added (hit index, keys), the values it overwrote (hit
    // index, key, previous value), the hits it skipped (hit index, reason) and the error it stopped at, collected in
    // spectrum order below. On an error, the spectra after the failing one are restored, as if the loop had stopped there.
    enum SkipReason { NO_EVIDENCE, UNKNOWN_TARGET_DECOY };
    std::vector<std::vector<std::pair<size_t, std::set<std::string>>>> added(peptide_ids.size());
    std::vector<std::vector<std::tuple<size_t, std::string, DataValue>>> overwritten(peptide_ids.size());
    std::vector<std::vector<std::pair<size_t, SkipReason>>> skipped_hits(peptide_ids.size());
    std::vector<std::exception_ptr> errors(peptide_ids.size());

    // The meta values by name index once all names below are registered (by the first stamped hit): looking names up
    // takes a lock on the registry, which the threads would contend for. Unregistered names are looked up by name.
    std::unordered_map<std::string, UInt> name_index;
    auto exists = [&](const PeptideHit& hit, const std::string& key)
    {
      const auto it = name_index.find(key);
      return it == name_index.end() ? hit.metaValueExists(key) : hit.metaValueExists(it->second);
    };
    auto get = [&](const PeptideHit& hit, const std::string& key) -> const DataValue&
    {
      const auto it = name_index.find(key);
      return it == name_index.end() ? hit.getMetaValue(key) : hit.getMetaValue(it->second);
    };

    const std::vector<std::string> column_names = pinColumnNames(min_charge, max_charge);
    const std::string lookup_names[] = {"CalcMass", "IsotopeError", Constants::UserParam::ISOTOPE_ERROR};

    // Stamps the hits of one spectrum; returns whether it stamped any.
    auto stamp_spectrum = [&](size_t pid_index)
    {
      PeptideIdentification& pep_id = peptide_ids[pid_index];
      const size_t pid_index_1based = pid_index + 1;

      const std::string scan_identifier = getScanIdentifier(pep_id, pid_index_1based);
      const std::string file_identifier = getFileIdentifier(pep_id);

      const Int scan_number = SpectrumLookup::extractScanNumber(
        scan_identifier, scan_regex, /*no_error=*/true);

      double exp_mass = pep_id.getMZ();
      const double retention_time = pep_id.getRT();

      bool stamped = false;
      auto& hits = pep_id.getHits();
      for (size_t hit_index = 0; hit_index < hits.size(); ++hit_index)
      {
        PeptideHit& hit = hits[hit_index];
        std::set<std::string> added_keys;

        // Remember only keys introduced by this call. Existing values are not
        // cleanup candidates, even when the PIN calculation refreshes them.
        auto stamp_meta_value = [&](const std::string& key, const auto& value)
        {
          if (!exists(hit, key))
          {
            added_keys.insert(key);
          }
          else
          {
            overwritten[pid_index].emplace_back(hit_index, key, get(hit, key));
          }
          const auto it = name_index.find(key);
          if (it == name_index.end()) hit.setMetaValue(key, value);
          else hit.setMetaValue(it->second, value);
        };

        try
        {
          if (hit.getPeptideEvidences().empty())
          {
            skipped_hits[pid_index].emplace_back(hit_index, NO_EVIDENCE);
            continue;
          }
          // looked up once (by name, under the lock of the registry): hit.isDecoy() is target_decoy == DECOY
          const PeptideHit::TargetDecoyType target_decoy = hit.getTargetDecoyType();
          if (target_decoy == PeptideHit::TargetDecoyType::UNKNOWN)
          {
            skipped_hits[pid_index].emplace_back(hit_index, UNKNOWN_TARGET_DECOY);
            continue;
          }
          auto put = [&](int column, const auto& value) { stamp_meta_value(column_names[column], value); };
          auto exists_lookup = [&](int lookup) { return exists(hit, lookup_names[lookup]); };
          auto get_lookup = [&](int lookup) -> const DataValue& { return get(hit, lookup_names[lookup]); };
          pinHitColumns_<true>(hit, target_decoy == PeptideHit::TargetDecoyType::DECOY, file_identifier + scan_identifier,
                               scan_number, exp_mass, retention_time, enz, min_charge, max_charge,
                               exists_lookup, get_lookup, put);
        }
        catch (...)
        {
          // the keys added before the error are cleanup candidates as well
          if (!added_keys.empty()) added[pid_index].emplace_back(hit_index, std::move(added_keys));
          throw;
        }
        if (!added_keys.empty())
        {
          added[pid_index].emplace_back(hit_index, std::move(added_keys));
        }
        stamped = true;
      }
      return stamped;
    };

    // Up to the first stamped hit serially: it registers the new meta value names in the order of a serial loop, and
    // MetaInfo keeps meta values in the order of their name indices. The other spectra in parallel.
    size_t first_parallel = 0;
    try
    {
      while (first_parallel < peptide_ids.size() && !stamp_spectrum(first_parallel++)) {}
    }
    catch (...)
    {
      errors[first_parallel - 1] = std::current_exception();
      first_parallel = peptide_ids.size();
    }
    std::vector<std::string> names{"SpecId", "ScanNr", "Label", "CalcMass", "IsotopeError", Constants::UserParam::ISOTOPE_ERROR,
                                   "ExpMass", "deltamass", "retentiontime", "mass", "score", "peplen", "enzN", "enzC",
                                   "enzInt", "dm", "absdm", "Peptide", "Proteins"};
    for (int i = min_charge; i <= max_charge; ++i)
    {
      names.push_back("charge" + StringUtils::toStr(i));
    }
    for (const std::string& name : names)
    {
      const UInt index = MetaInfoInterface::metaRegistry().getIndex(name);
      if (index != UInt(-1)) name_index.emplace(name, index);
    }
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize pid_index = static_cast<SignedSize>(first_parallel); pid_index < static_cast<SignedSize>(peptide_ids.size()); ++pid_index)
    {
      try
      {
        stamp_spectrum(static_cast<size_t>(pid_index));
      }
      catch (...)
      {
        errors[pid_index] = std::current_exception();
      }
    }

    for (size_t pid_index = 0; pid_index < peptide_ids.size(); ++pid_index)
    {
      for (auto& [hit_index, keys] : added[pid_index])
      {
        added_meta_values.emplace_hint(added_meta_values.end(), std::make_pair(pid_index, hit_index), std::move(keys));
      }
      for (const auto& [hit_index, reason] : skipped_hits[pid_index])
      {
        if (reason == NO_EVIDENCE)
        {
          OPENMS_LOG_WARN << "PSM (PeptideHit) without protein reference found. "
                  << "This may indicate incomplete mapping during PeptideIndexing "
                     "(e.g., wrong enzyme settings). Will skip this PSM." << endl;
        }
        else
        {
          OPENMS_LOG_WARN << "PSM without target/decoy information found. "
                  << "This may indicate incomplete mapping during PeptideIndexing "
                     "(e.g., wrong decoy prefix settings). Will skip this PSM." << endl;
        }
        skipped.emplace_hint(skipped.end(), pid_index, hit_index);
      }
      if (errors[pid_index])
      {
        // restore the spectra stamped after the failing one
        for (size_t later = pid_index + 1; later < peptide_ids.size(); ++later)
        {
          auto& hits = peptide_ids[later].getHits();
          for (const auto& [hit_index, keys] : added[later])
          {
            for (const std::string& key : keys) hits[hit_index].removeMetaValue(key);
          }
          for (const auto& [hit_index, key, value] : overwritten[later])
          {
            hits[hit_index].setMetaValue(key, value);
          }
        }
        std::rethrow_exception(errors[pid_index]);
      }
    }
    return skipped;
  }

  double PercolatorInfile::getFeatureValue(const DataValue& value, const std::string& feature)
  {
    double number = std::numeric_limits<double>::quiet_NaN();
    const DataValue::DataType type = value.valueType();
    if (type == DataValue::INT_VALUE || type == DataValue::DOUBLE_VALUE)
    {
      number = static_cast<double>(value);
    }
    else if (type == DataValue::STRING_VALUE)
    {
      // the .pin writer prints the string as it is, and the executable parses it as a number
      try
      {
        number = StringUtils::toDouble(value.toString());
      }
      catch (const Exception::ConversionError&)
      {
        // not a number: reported below
      }
    }
    // the executable also stops at a feature that is not finite ("Reached strange feature")
    if (!std::isfinite(number))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Percolator feature '" + feature + "' does not have a finite numeric value", value.toString());
    }
    return number;
  }

  Size PercolatorInfile::buildFeatureMatrix(const PeptideIdentificationList& peptide_ids,
    const StringList& features,
    const std::string& enz,
    int min_charge,
    int max_charge,
    RescoreInput& input,
    std::vector<std::pair<size_t, size_t>>& rows)
  {
    input = RescoreInput();
    rows.clear();
    input.feature_names = features;

    // The columns: a PIN column (its PinColumn_) or a meta value of the hit (ENGINE, read by name index; an
    // unregistered name, index -1, exists on no hit), as the stamped feature rows were read.
    constexpr int ENGINE = -1;
    struct Column
    {
      int pin;
      UInt index;
    };
    const std::vector<std::string> column_names = pinColumnNames(min_charge, max_charge);
    std::vector<Column> columns;
    for (const std::string& feature : features)
    {
      const auto it = std::find(column_names.begin(), column_names.end(), feature);
      if (it == column_names.end())
      {
        columns.push_back({ENGINE, MetaInfoInterface::metaRegistry().getIndex(feature)});
        continue;
      }
      const int pin = static_cast<int>(it - column_names.begin());
      if (pin == PIN_SPEC_ID || pin == PIN_PEPTIDE || pin == PIN_PROTEINS)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "PIN column '" + feature + "' is not a numeric feature", feature);
      }
      columns.push_back({pin, UInt(-1)});
    }
    const size_t n_features = columns.size();
    if (peptide_ids.empty()) return 0;

    const MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const UInt td_index = registry.getIndex("target_decoy");
    const UInt lookup_index[] = {registry.getIndex("CalcMass"), registry.getIndex("IsotopeError"),
                                 registry.getIndex(Constants::UserParam::ISOTOPE_ERROR)};

    // Scan-number regex is derived from the first pid's scan identifier.
    const std::string sid = getScanIdentifier(peptide_ids[0], 0);
    RegularExpression scan_regex = RegularExpression(SpectrumLookup::getRegExFromNativeID(sid));

    // The hits of spectrum (PeptideIdentification) i are hit_offset[i] onwards. Their state (first pass): skipped
    // (0) or stamped (STAMPED, i.e. evidence and target/decoy type), a decoy (DECOY), a row (ROW: every feature).
    constexpr unsigned char STAMPED = 1, DECOY = 2, ROW = 4;
    const size_t n_spectra = peptide_ids.size();
    std::vector<size_t> hit_offset(n_spectra + 1, 0);
    for (size_t i = 0; i < n_spectra; ++i)
    {
      hit_offset[i + 1] = hit_offset[i] + peptide_ids[i].getHits().size();
    }
    std::vector<unsigned char> state(hit_offset[n_spectra], 0);
    // Per spectrum: the first row (row_offset), the file identifier (spectra with rows), the skipped hits, and the
    // first error of the stamping (with the hit it stopped at; SIZE_MAX: none) and of the rows, as stamping and
    // then the reading of the stamped feature rows raised them.
    enum SkipReason { NO_EVIDENCE, UNKNOWN_TARGET_DECOY };
    constexpr size_t NONE = std::numeric_limits<size_t>::max();
    std::vector<size_t> row_offset(n_spectra + 1, 0);
    std::vector<std::string> file_identifier(n_spectra);
    std::vector<std::vector<std::pair<size_t, SkipReason>>> skipped_hits(n_spectra);
    std::vector<size_t> stamp_error_hit(n_spectra, NONE);
    std::vector<std::exception_ptr> stamp_errors(n_spectra), row_errors(n_spectra);

    // First pass: the hits that are stamped and rows, the number of rows of each spectrum.
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize si = 0; si < static_cast<SignedSize>(n_spectra); ++si)
    {
      const size_t i = static_cast<size_t>(si);
      const std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      size_t n_rows = 0;
      for (size_t j = 0; j < hits.size(); ++j)
      {
        try
        {
          const PeptideHit& hit = hits[j];
          if (hit.getPeptideEvidences().empty())
          {
            skipped_hits[i].emplace_back(j, NO_EVIDENCE);
            continue;
          }
          const PeptideHit::TargetDecoyType target_decoy = targetDecoyType(hit, td_index);
          if (target_decoy == PeptideHit::TargetDecoyType::UNKNOWN)
          {
            skipped_hits[i].emplace_back(j, UNKNOWN_TARGET_DECOY);
            continue;
          }
          unsigned char s = STAMPED;
          if (target_decoy == PeptideHit::TargetDecoyType::DECOY) s |= DECOY;
          if (std::all_of(columns.begin(), columns.end(),
                          [&](const Column& c) { return c.pin != ENGINE || hit.metaValueExists(c.index); }))
          {
            s |= ROW;
            ++n_rows;
          }
          state[hit_offset[i] + j] = s;
        }
        catch (...)
        {
          stamp_error_hit[i] = j;
          stamp_errors[i] = std::current_exception();
          break;
        }
      }
      row_offset[i + 1] = n_rows;
      if (n_rows > 0)
      {
        try
        {
          file_identifier[i] = getFileIdentifier(peptide_ids[i]);
        }
        catch (...)
        {
          row_errors[i] = std::current_exception();
        }
      }
    }

    // The first row of each spectrum, and the spectrum files in order of appearance, as Percolator numbers the
    // FileName column of the pin file: it identifies a spectrum by file, scan number and precursor m/z (e.g. for
    // target-decoy competition), so spectra of different inputs with the same scan number must not share a file.
    std::vector<int> spec_file(n_spectra, 0);
    {
      std::map<std::string, int> spec_files;
      for (size_t i = 0; i < n_spectra; ++i)
      {
        if (row_offset[i + 1] > 0)
        {
          spec_file[i] = spec_files.try_emplace(file_identifier[i], static_cast<int>(spec_files.size())).first->second;
        }
        row_offset[i + 1] += row_offset[i];
      }
    }
    file_identifier = std::vector<std::string>();
    const size_t n_rows = row_offset[n_spectra];
    // The rows are allocated only without an error: an allocation failure must not replace the error rethrown below.
    auto has_error = [](const std::vector<std::exception_ptr>& errors)
    {
      return std::any_of(errors.begin(), errors.end(), [](const std::exception_ptr& e) { return bool(e); });
    };
    const bool fill = !has_error(stamp_errors) && !has_error(row_errors);
    std::vector<char> is_decoy;
    if (fill)
    {
      input.feature_matrix.resize(n_rows * n_features);
      is_decoy.resize(n_rows);
      input.scan_numbers.resize(n_rows);
      input.spec_file_numbers.resize(n_rows);
      input.exp_masses.resize(n_rows);
      input.calc_masses.resize(n_rows);
      rows.resize(n_rows);
    }

    // Second pass: the PIN columns of every stamped hit (the stamping's errors) and the rows (the rows' errors: a
    // feature without a finite numeric value, read in the order of the features up to a missing one). Per thread:
    // the PIN columns of a hit, and the features of a hit that is not a row.
#ifdef _OPENMP
    const size_t n_buffers = static_cast<size_t>(omp_get_max_threads());
#else
    const size_t n_buffers = 1;
#endif
    const size_t pin_width = PIN_CHARGE + static_cast<size_t>(std::max(0, max_charge - min_charge + 1));
    std::vector<double> pin_buffer(n_buffers * pin_width), scratch_buffer(n_buffers * n_features);
    const std::string no_spec_id;
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize si = 0; si < static_cast<SignedSize>(n_spectra); ++si)
    {
      const size_t i = static_cast<size_t>(si);
#ifdef _OPENMP
      const size_t thread = static_cast<size_t>(omp_get_thread_num());
#else
      const size_t thread = 0;
#endif
      double* pin_values = pin_buffer.data() + thread * pin_width;
      double* scratch_row = scratch_buffer.data() + thread * n_features;
      const PeptideIdentification& pid = peptide_ids[i];
      const std::vector<PeptideHit>& hits = pid.getHits();
      size_t row = row_offset[i];
      try
      {
        const std::string scan_identifier = getScanIdentifier(pid, i + 1);
        const Int scan_number = SpectrumLookup::extractScanNumber(scan_identifier, scan_regex, /*no_error=*/true);
        const double exp_mass = pid.getMZ();
        const double retention_time = pid.getRT();
        // the hits up to an error of the first pass (a target/decoy value)
        const size_t n_hits = std::min(hits.size(), stamp_error_hit[i]);
        for (size_t j = 0; j < n_hits; ++j)
        {
          const unsigned char s = state[hit_offset[i] + j];
          if (!(s & STAMPED)) continue;
          const PeptideHit& hit = hits[j];
          double calc_mass;
          try
          {
            auto exists = [&](int lookup)
            {
              return lookup_index[lookup] != UInt(-1) && hit.metaValueExists(lookup_index[lookup]);
            };
            auto get = [&](int lookup) -> const DataValue& { return hit.getMetaValue(lookup_index[lookup]); };
            auto put = [&](int column, const auto& value)
            {
              if constexpr (std::is_arithmetic_v<std::decay_t<decltype(value)>>) pin_values[column] = static_cast<double>(value);
            };
            calc_mass = pinHitColumns_<false>(hit, (s & DECOY) != 0, no_spec_id, scan_number, exp_mass, retention_time,
                                              enz, min_charge, max_charge, exists, get, put);
            pin_values[PIN_CALC_MASS] = calc_mass; // (also when read from the hit)
          }
          catch (...)
          {
            stamp_error_hit[i] = j;
            stamp_errors[i] = std::current_exception();
            break;
          }
          if (row_errors[i]) continue; // the stamping goes on: its errors come first
          try
          {
            // the features in the order of the feature set, up to the first missing one
            const bool is_row = fill && (s & ROW);
            double* out = is_row ? input.feature_matrix.data() + row * n_features : scratch_row;
            bool complete = true;
            for (size_t f = 0; f < n_features; ++f)
            {
              const Column& c = columns[f];
              double value;
              if (c.pin == ENGINE)
              {
                if (!hit.metaValueExists(c.index))
                {
                  complete = false;
                  break;
                }
                value = getFeatureValue(hit.getMetaValue(c.index), features[f]);
              }
              else
              {
                value = pin_values[c.pin];
                // (a value that is not finite is reported as for its stamped meta value)
                if (!std::isfinite(value)) getFeatureValue(DataValue(value), features[f]);
              }
              out[f] = value;
            }
            if (!complete) continue;
            if (!std::isfinite(calc_mass)) getFeatureValue(DataValue(calc_mass), "CalcMass");
            if (!is_row) continue;
            is_decoy[row] = (s & DECOY) != 0;
            input.scan_numbers[row] = scan_number;
            input.spec_file_numbers[row] = spec_file[i];
            input.exp_masses[row] = exp_mass;
            input.calc_masses[row] = calc_mass;
            rows[row] = {i, j};
            ++row;
          }
          catch (...)
          {
            row_errors[i] = std::current_exception();
          }
        }
      }
      catch (...)
      {
        // before the first hit, as the stamping of the spectrum
        stamp_error_hit[i] = 0;
        stamp_errors[i] = std::current_exception();
      }
    }

    // As stamping its first hit: register the names of the PIN columns in that order (MetaInfo orders meta values
    // by name index, also those that later code adds under these names), unless the stamping failed before.
    for (size_t i = 0; i < n_spectra; ++i)
    {
      size_t first = 0;
      while (first < hit_offset[i + 1] - hit_offset[i] && !(state[hit_offset[i] + first] & STAMPED)) ++first;
      if (stamp_error_hit[i] <= first) break; // (also no stamped hit before the error)
      if (first < hit_offset[i + 1] - hit_offset[i])
      {
        for (const std::string& name : pinColumnNamesInStampingOrder(min_charge, max_charge))
        {
          MetaInfoInterface::metaRegistry().registerName(name);
        }
        break;
      }
    }

    // The warnings about skipped hits and the stamping's first error in spectrum order, then the rows' first error.
    Size skipped = 0;
    for (size_t i = 0; i < n_spectra; ++i)
    {
      for (const auto& [hit_index, reason] : skipped_hits[i])
      {
        if (hit_index >= stamp_error_hit[i]) break;
        if (reason == NO_EVIDENCE)
        {
          OPENMS_LOG_WARN << "PSM (PeptideHit) without protein reference found. "
                  << "This may indicate incomplete mapping during PeptideIndexing "
                     "(e.g., wrong enzyme settings). Will skip this PSM." << endl;
        }
        else
        {
          OPENMS_LOG_WARN << "PSM without target/decoy information found. "
                  << "This may indicate incomplete mapping during PeptideIndexing "
                     "(e.g., wrong decoy prefix settings). Will skip this PSM." << endl;
        }
        ++skipped;
      }
      if (stamp_errors[i]) std::rethrow_exception(stamp_errors[i]);
    }
    for (const std::exception_ptr& error : row_errors)
    {
      if (error) std::rethrow_exception(error);
    }
    input.is_decoy.assign(is_decoy.begin(), is_decoy.end());
    return skipped;
  }

  void PercolatorInfile::refreshPinMetaValues(PeptideIdentificationList& peptide_ids,
    const std::string& enz,
    int min_charge,
    int max_charge)
  {
    // The PIN columns that stamping overwrites (all but CalcMass, which it reads): their name indices, if registered.
    // A hit can only carry those that are.
    const MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const std::vector<std::string> column_names = pinColumnNames(min_charge, max_charge);
    std::vector<UInt> index(column_names.size(), UInt(-1));
    for (size_t c = 0; c < column_names.size(); ++c)
    {
      if (c != PIN_CALC_MASS) index[c] = registry.getIndex(column_names[c]);
    }
    if (peptide_ids.empty() || std::all_of(index.begin(), index.end(), [](UInt i) { return i == UInt(-1); })) return;

    const UInt td_index = registry.getIndex("target_decoy");
    const UInt lookup_index[] = {registry.getIndex("CalcMass"), registry.getIndex("IsotopeError"),
                                 registry.getIndex(Constants::UserParam::ISOTOPE_ERROR)};
    const std::string sid = getScanIdentifier(peptide_ids[0], 0);
    RegularExpression scan_regex = RegularExpression(SpectrumLookup::getRegExFromNativeID(sid));

    // A hit whose computation fails is left as it is: the same computation fails in buildFeatureMatrix(), which
    // reports the error in order.
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize si = 0; si < static_cast<SignedSize>(peptide_ids.size()); ++si)
    {
      PeptideIdentification& pid = peptide_ids[si];
      std::optional<std::tuple<std::string, Int>> spectrum; // SpecId and ScanNr, when a hit needs them
      for (PeptideHit& hit : pid.getHits())
      {
        try
        {
          if (std::none_of(index.begin(), index.end(), [&](UInt i) { return i != UInt(-1) && hit.metaValueExists(i); }))
          {
            continue;
          }
          if (hit.getPeptideEvidences().empty()) continue;
          const PeptideHit::TargetDecoyType target_decoy = targetDecoyType(hit, td_index);
          if (target_decoy == PeptideHit::TargetDecoyType::UNKNOWN) continue;
          if (!spectrum)
          {
            const std::string scan_identifier = getScanIdentifier(pid, static_cast<size_t>(si) + 1);
            spectrum.emplace(getFileIdentifier(pid) + scan_identifier,
                             SpectrumLookup::extractScanNumber(scan_identifier, scan_regex, /*no_error=*/true));
          }
          auto exists = [&](int lookup)
          {
            return lookup_index[lookup] != UInt(-1) && hit.metaValueExists(lookup_index[lookup]);
          };
          auto get = [&](int lookup) -> const DataValue& { return hit.getMetaValue(lookup_index[lookup]); };
          auto put = [&](int column, const auto& value)
          {
            if (index[column] != UInt(-1) && hit.metaValueExists(index[column])) hit.setMetaValue(index[column], value);
          };
          pinHitColumns_<true>(hit, target_decoy == PeptideHit::TargetDecoyType::DECOY, std::get<0>(*spectrum),
                               std::get<1>(*spectrum), pid.getMZ(), pid.getRT(), enz, min_charge, max_charge,
                               exists, get, put);
        }
        catch (...)
        {
        }
      }
    }
  }

  TextFile PercolatorInfile::preparePin_(
    const PeptideIdentificationList& peptide_ids,
    const StringList& feature_set,
    const std::string& enz,
    int min_charge,
    int max_charge)
  {
    TextFile txt;
    if (peptide_ids.empty())
    {
      txt.addLine(ListUtils::concatenate(feature_set, "\t"));
      OPENMS_LOG_WARN << "No identifications provided. Creating empty percolator input." << endl;
      return txt;
    }

    // Work on a local copy so we don't mutate caller state.
    PeptideIdentificationList stamped(peptide_ids);
    const auto skipped = stampPinFeaturesOnHits(stamped, enz, min_charge, max_charge);

    // Rows: the kept hits that carry every column of the feature set.
    std::vector<std::pair<size_t, size_t>> rows;
    std::set<std::string> files;
    size_t missing_meta_value_count = 0;
    set<std::string> missing_meta_values;
    for (size_t pid_idx = 0; pid_idx < stamped.size(); ++pid_idx)
    {
      const auto& hits = stamped[pid_idx].getHits();
      for (size_t hit_idx = 0; hit_idx < hits.size(); ++hit_idx)
      {
        if (skipped.contains({pid_idx, hit_idx})) continue;
        const PeptideHit& hit = hits[hit_idx];
        bool complete = true;
        for (const std::string& feat : feature_set)
        {
          if (!hit.metaValueExists(feat))
          {
            complete = false;
            missing_meta_values.insert(feat);
          }
        }
        if (!complete)
        {
          ++missing_meta_value_count;
          continue;
        }
        rows.emplace_back(pid_idx, hit_idx);
        files.insert(getFileIdentifier(stamped[pid_idx]));
      }
    }

    // Percolator numbers the spectrum files from an optional 'FileName' column, without which all PSMs count as
    // one file: spectra of different files with the same ScanNr and ExpMass would be one spectrum. Optional columns
    // precede the features, so it follows the leading optional columns of the feature set.
    StringList columns = feature_set;
    std::optional<size_t> file_column;
    if (files.size() > 1)
    {
      auto lower = [](std::string name)
      {
        std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return name;
      };
      const std::set<std::string> optional_columns {"scannr", "expmass", "calcmass", "rt", "retentiontime", "filename", "spectrafile"};
      const bool named = std::any_of(columns.begin(), columns.end(), [&](const std::string& c) { const auto l = lower(c); return l == "filename" || l == "spectrafile"; });
      if (!named)
      {
        size_t position = std::min<size_t>(2, columns.size()); // after SpecId and Label
        while (position < columns.size() && optional_columns.contains(lower(columns[position]))) ++position;
        columns.insert(columns.begin() + position, "FileName");
        file_column = position;
      }
    }
    txt.addLine(ListUtils::concatenate(columns, "\t"));

    for (const auto& [pid_idx, hit_idx] : rows)
    {
      const PeptideHit& hit = stamped[pid_idx].getHits()[hit_idx];
      StringList feats;
      feats.reserve(columns.size());
      for (const std::string& feat : feature_set)
      {
        feats.push_back(hit.getMetaValue(feat).toString());
      }
      if (file_column) feats.insert(feats.begin() + *file_column, getFileIdentifier(stamped[pid_idx]));
      txt.addLine(ListUtils::concatenate(feats, "\t"));
    }

    if (missing_meta_value_count != 0)
    {
      OPENMS_LOG_WARN << "There were peptide hits with missing features/meta values. Skipped peptide hits: "
                      << missing_meta_value_count << endl;
      OPENMS_LOG_WARN << "Names of missing meta values: " << endl;
      for (const auto& f : missing_meta_values) OPENMS_LOG_WARN << f << endl;
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
  Size PercolatorInfile::countEnzymatic_(const std::string& peptide, const string& enz)
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
