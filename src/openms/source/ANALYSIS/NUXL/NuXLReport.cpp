// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLReport.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <boost/range/adaptor/reversed.hpp>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>

using namespace std;

namespace OpenMS
{
  using Internal::IDBoostGraph;

  String NuXLReportRow::getString(const String& separator) const
  {
    StringList sl;

    // rt mz
    sl << String::number(rt, 3) 
       << String::number(original_mz, 4);

    // id if available
    if (no_id)
    {
      for (Size i = 0; i != 12; ++i) sl << "";
    }
    else
    {
      sl << accessions 
         << peptide 
         << NA 
         << String(charge) 
         << String(score)
         << String(rank)
         << best_localization_score 
         << localization_scores 
         << best_localization
         << String::number(peptide_weight, 4) << String::number(NA_weight, 4)
         << String::number(peptide_weight + NA_weight, 4);
    }

    // write out meta value columns
    for (const String& v : meta_values)
    {
      sl << v;
    }

    // marker ions
    for (auto it = marker_ions.cbegin(); it != marker_ions.cend(); ++it)
    {
      for (Size i = 0; i != it->second.size(); ++i)
      {
        sl << String::number(it->second[i].second * 100.0, 2);
      }
    }

    // id error and multiple charged mass
    if (no_id)
    {
      for (Size i = 0; i != 7; ++i) sl << "";
    }
    else
    {
      // error
      sl << String::number(abs_prec_error, 4)
         << String::number(rel_prec_error, 1);

      // weight
      sl << String::number(m_H, 4)
         << String::number(m_2H, 4)
         << String::number(m_3H, 4)
         << String::number(m_4H, 4);

      sl << fragment_annotation;
    }

    return ListUtils::concatenate(sl, separator);
  }

  String NuXLReportRowHeader::getString(const String& separator, const StringList& meta_values_to_export)
  {
    StringList sl;
    sl << "#RT" 
       << "m/z" 
       << "proteins" 
       << "peptide"
       << "NA" 
       << "charge" 
       << "score"
       << "rank"
       << "best localization score" 
       << "localization scores" 
       << "best localization(s)"
       << "peptide weight" 
       << "NA weight" 
       << "cross-link weight";

    for (const String& s : meta_values_to_export) sl << s;

    // marker ion fields
    NuXLMarkerIonExtractor::MarkerIonsType marker_ions = NuXLMarkerIonExtractor::extractMarkerIons(PeakSpectrum(), 0.0); // call only to generate header entries
    for (auto const & ma : marker_ions)
    {
      for (Size i = 0; i != ma.second.size(); ++i)
      {
        sl << String(ma.first + "_" + ma.second[i].first);
      }
    }
    sl << "abs prec. error Da" 
       << "rel. prec. error ppm" 
       << "M+H" 
       << "M+2H" 
       << "M+3H" 
       << "M+4H"
       << Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM; 
    return ListUtils::concatenate(sl, separator);
  }

  std::vector<NuXLReportRow> NuXLReport::annotate(const PeakMap& spectra, std::vector<PeptideIdentification>& peptide_ids, const StringList& meta_values_to_export, double marker_ions_tolerance)
  {
    std::map<Size, Size> map_spectra_to_id;
    for (Size i = 0; i != peptide_ids.size(); ++i)
    {
      OPENMS_PRECONDITION(!peptide_ids[i].getHits().empty(), "Error: no empty peptide ids allowed.");
      Size scan_index = (unsigned int)peptide_ids[i].getMetaValue("scan_index");
      map_spectra_to_id[scan_index] = i;
    }

    std::vector<NuXLReportRow> csv_rows;

    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      std::vector<Precursor> precursor = s_it->getPrecursors();

      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
      if (s_it->getMSLevel() == 2 && precursor.size() == 1)
      {
        Size charge = precursor[0].getCharge();
        double mz = precursor[0].getMZ();
        NuXLMarkerIonExtractor::MarkerIonsType marker_ions = NuXLMarkerIonExtractor::extractMarkerIons(*s_it, marker_ions_tolerance);

        double rt = s_it->getRT();

        NuXLReportRow row;

        // case 1: no peptide identification: store rt, mz, charge and marker ion intensities
        if (map_spectra_to_id.find(scan_index) == map_spectra_to_id.end())
        {
          row.no_id = true;
          row.rt = rt;
          row.original_mz = mz;
          row.charge = charge;
          row.marker_ions = marker_ions;
          csv_rows.push_back(row);
          continue;
        }

        PeptideIdentification& pi = peptide_ids[map_spectra_to_id[scan_index]];
        std::vector<PeptideHit>& phs = pi.getHits();

        // case 2: identification data present for spectrum
        Size rank(0);
        for (PeptideHit& ph : phs)
        {
          ++rank;
          
          for (const String& meta_key : meta_values_to_export)
          {
            row.meta_values.emplace_back(ph.getMetaValue(meta_key).toString());
          }

          PeptideHit::PeakAnnotation::writePeakAnnotationsString_(row.fragment_annotation, ph.getPeakAnnotations()); 

          // total weight = precursor NA weight + peptide weight
          // this ensures that sequences with additional reported partial loss match the total weight
          // Note that the partial loss is only relevent on the MS2 and would otherwise be added to the totalweight
          String sequence_string = ph.getSequence().toString();

          const AASequence sequence = AASequence::fromString(sequence_string);

          double peptide_weight = sequence.getMonoWeight();
          String rna_name = ph.getMetaValue("NuXL:NA");
          double rna_weight = ph.getMetaValue("NuXL:NA_MASS_z0");
          int isotope_error = ph.getMetaValue("isotope_error");
          // crosslink weight for different charge states
          double weight_z1 = (peptide_weight + rna_weight + 1.0 * Constants::PROTON_MASS_U);
          double weight_z2 = (peptide_weight + rna_weight + 2.0 * Constants::PROTON_MASS_U) / 2.0;
          double weight_z3 = (peptide_weight + rna_weight + 3.0 * Constants::PROTON_MASS_U) / 3.0;
          double weight_z4 = (peptide_weight + rna_weight + 4.0 * Constants::PROTON_MASS_U) / 4.0;

          double xl_weight = peptide_weight + rna_weight;
          double theo_mz = (xl_weight + static_cast<double>(charge) * Constants::PROTON_MASS_U) / (double)charge;

          double corr_mz = mz - (double)isotope_error * Constants::PROTON_MASS_U / (double)charge;
          double absolute_difference = theo_mz - corr_mz;
          double ppm_difference =  Math::getPPM(corr_mz, theo_mz);

          String protein_accessions;
          std::set<String> accs = ph.extractProteinAccessionsSet();

          // concatenate set into String
          for (std::set<String>::const_iterator a_it = accs.begin(); a_it != accs.end(); ++a_it)
          {
            if (a_it != accs.begin())
            {
              protein_accessions += ",";
            }
            protein_accessions += *a_it;
          }

          row.no_id = false;
          row.rt = rt;
          row.original_mz = mz;
          row.accessions = protein_accessions;
          row.NA = rna_name;
          row.peptide = ph.getSequence().toString();
          row.charge = charge;
          row.score = ph.getScore();
          row.peptide_weight = peptide_weight;
          row.NA_weight = rna_weight;
          row.xl_weight = peptide_weight + rna_weight;
          row.rank = rank;

          ph.setMetaValue("NuXL:peptide_mass_z0", DataValue(peptide_weight));
          ph.setMetaValue("NuXL:xl_mass_z0", xl_weight);

          for (NuXLMarkerIonExtractor::MarkerIonsType::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
          {
            for (Size i = 0; i != it->second.size(); ++i)
            {
              ph.setMetaValue(it->first + "_" + it->second[i].first, static_cast<double>(it->second[i].second * 100.0));
            }
          }

          row.marker_ions = marker_ions;
          row.abs_prec_error = absolute_difference;
          row.rel_prec_error = ppm_difference;
          row.m_H = weight_z1;
          row.m_2H = weight_z2;
          row.m_3H = weight_z3;
          row.m_4H = weight_z4;

          if (ph.metaValueExists("NuXL:best_localization_score") && 
              ph.metaValueExists("NuXL:localization_scores") && 
              ph.metaValueExists("NuXL:best_localization"))
          {
            row.best_localization_score = ph.getMetaValue("NuXL:best_localization_score");
            row.localization_scores = ph.getMetaValue("NuXL:localization_scores");
            row.best_localization = ph.getMetaValue("NuXL:best_localization");;
          }

          ph.setMetaValue("NuXL:Da difference", (double)absolute_difference);
          ph.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM, (double)ppm_difference);
          ph.setMetaValue("NuXL:z1 mass", (double)weight_z1);
          ph.setMetaValue("NuXL:z2 mass", (double)weight_z2);
          ph.setMetaValue("NuXL:z3 mass", (double)weight_z3);
          ph.setMetaValue("NuXL:z4 mass", (double)weight_z4);
          csv_rows.push_back(row);
      }
    }
  }
  return csv_rows;
}

  // crosslink efficiency = frequency of the crosslinked amino acid / frequency of the amino acid in all crosslink spectrum matches
  map<char, double> NuXLProteinReport::getCrossLinkEfficiency(const vector<PeptideIdentification>& peps)
  {
    map<char, double> aa_xl_freq;
    map<char, double> aa_freq;

    // store modification statistic for every protein    
    for (const PeptideIdentification& pep : peps)
    {
      auto& hits = pep.getHits();
      if (hits.empty()) continue;
      const PeptideHit& ph = hits[0]; // only consider top hit
      if (ph.getMetaValue("target_decoy") == "decoy" || ph.getMetaValue("NuXL:isXL") == "false") continue;
      const int best_localization = ph.getMetaValue("NuXL:best_localization_position");

      if (best_localization >= 0)
      {
        const AASequence& aas = ph.getSequence();
        char c = aas.toUnmodifiedString()[best_localization];
        aa_xl_freq[c] += 1;
        for (char c : aas.toUnmodifiedString())
        {
          aa_freq[c] += 1;
        }
      }
    }

    double xl_sum{};
    for (const auto& m : aa_xl_freq) { xl_sum += m.second; }
    for (auto& m : aa_xl_freq) { m.second /= xl_sum; }

    double aa_sum{};
    for (const auto& m : aa_freq) { aa_sum += m.second; }
    for (auto& m : aa_freq) { m.second /= aa_sum; }
    
    for (auto& m : aa_xl_freq) {m.second /= aa_freq[m.first]; }    

    return aa_xl_freq;
  }

  // returns map of adduct to counts
  map<String, size_t> NuXLProteinReport::countAdducts(const vector<PeptideIdentification>& peps)
  {
    map<String, size_t> adduct2count;
    for (const PeptideIdentification& pep : peps)
    {
      auto& hits = pep.getHits();
      if (hits.empty()) continue;
      const PeptideHit& ph = hits[0]; // only consider top hit
      const String NA = ph.getMetaValue("NuXL:NA", String("none"));
      adduct2count[NA] += 1;
    }
    return adduct2count;
  }
/*
Output format:
+----------------------+-----+--------+---------+------+--------------------------+---------------------+---------------------------+-----------------------+------------------------+------------------------------+------------------------------+-------------------------------+------------------------------+---------------------------+---------------------------+----------------------------------+----------------------------------+---------------+-------------------------+
|      accession       |  AA |   pos. |   start |  end |  adducts (loc. + unique) |  NT (loc. + unique) |   charges (loc. + unique) |  CSMs (loc. + unique) |   CSMs (loc. + shared) |   precursors (loc. + unique) |   precursors (loc. + shared) |   adducts (\wo loc. + unique) |  charges (\wo loc. + unique) |  CSMs (\wo loc. + unique) |  CSMs (\wo loc. + shared) |   precursors (\wo loc. + unique) |   precursors (\wo loc. + shared) |   ambiguities |         peptide         | peptide-XL q-value
+----------------------+-----+--------+---------+------+--------------------------+---------------------+---------------------------+-----------------------+------------------------+------------------------------+------------------------------+-------------------------------+------------------------------+---------------------------+---------------------------+----------------------------------+----------------------------------+---------------+-------------------------+
| sp|P19338|NUCL_HUMAN |   P |    302 |     297 |  317 |  U                       |  U                  |                         3 |                     1 |                      0 |                            1 |                            0 |                               |                              |                         0 |                         0 |                                0 |                                0 |               |   VEGTEPTTAFNLFVGNLNFNK |  0.0
| sp|P19338|NUCL_HUMAN |   F |    309 |     297 |  317 |  U,U-H2O1                |   U                 |                       2,3 |                     3 |                      0 |                            3 |                            0 |                               |                              |                         0 |                         0 |                                0 |                                0 |               |   VEGTEPTTAFNLFVGNLNFNK |  0.0
+----------------------+-----+--------+---------+------+--------------------------+---------------------+---------------------------+-----------------------+------------------------+------------------------------+------------------------------+-------------------------------+------------------------------+---------------------------+---------------------------+----------------------------------+----------------------------------+---------------+-------------------------+
*/
  
  struct AALevelLocalization
  {
    struct LocalizedXL
    {
      std::string adduct;
      std::string NT;
      int charge = 0;
    };
    // Note protein accession and position are stored in ProteinReportEntry
    string AA;

    // note: we use a vector as entries might be duplicated and we need to obtain proper counts
    std::map<std::string, vector<LocalizedXL>> peptide2XL; // observed peptide -> adduct,NA,charge tuples    
  };

  struct RegionLevelLocalization
  {
    struct UnlocalizedXL
    {
      std::string adduct;
      int charge = 0;
    };

    // note: we use a vector as entries might be duplicated and we need to obtain proper counts
    std::map<std::string, vector<UnlocalizedXL>> peptide2unlocalizedXL; 
  };

  // all localization information for protein accession
  struct ProteinReport
  {
    String sequence; //< the protein sequence
    size_t CSMs_of_shared_peptides = 0; // XL spectral count of shared peptides
    size_t CSMs_of_unique_peptides = 0; // XL spectral count of unique peptides
    map<size_t, AALevelLocalization> aa_level_localization; // position in protein to loc info
    map<pair<size_t, size_t>, RegionLevelLocalization> region_level_localization;
  };

  // all proteins
  using ProteinsReport = map<std::string, ProteinReport>; //< protein accession to details
  std::unordered_map<String, double> peptide_seq2XLFDR;

  ProteinsReport getProteinReportEntries(
//    vector<ProteinIdentification>& prot_ids, 
    const vector<PeptideIdentification>& peps,
    const map<String, ProteinHit*>& acc2protein_targets,
    const std::map<string, set<string>>& peptide2proteins
    )
  {
    ProteinsReport report; // map accession to reporting data

    // go through all CSMs and collect information for individual XL sites
    for (const PeptideIdentification& pep : peps)
    {
      auto& hits = pep.getHits();

      if (hits.empty()) continue;

      const PeptideHit& ph = hits[0]; // only consider top hit
      const int best_localization = ph.getMetaValue("NuXL:best_localization_position");        
      const String& NA = ph.getMetaValue("NuXL:NA"); // adduct
      const String& NT = ph.getMetaValue("NuXL:NT"); // XLed nucleotide
      const int charge = ph.getCharge();
      const AASequence& peptide_sequence = ph.getSequence();
      
      // get mapping of peptide sequence to protein(s)
      const std::vector<PeptideEvidence>& ph_evidences = ph.getPeptideEvidences();
      const std::string peptide_sequence_string = peptide_sequence.toUnmodifiedString();

      // the peptide-level FDR in the group of cross-linked peptides
      double peptide_XL_level_qvalue = (double)ph.getMetaValue(Constants::UserParam::PEPTIDE_Q_VALUE, 0.0);
      peptide_seq2XLFDR[peptide_sequence_string] = peptide_XL_level_qvalue;

      // loop over all target proteins the peptide maps to
      const std::set<std::string> proteins = peptide2proteins.at(peptide_sequence_string);
      const bool is_unique = proteins.size() == 1;

      for (const String& acc : proteins)
      {
        // add basic protein information first time we encounter a protein accession
        int protein_length{};
        std::string protein_sequence;
        if (auto it = report.find(acc); it == report.end())
        {          
          auto& protein_sequence = acc2protein_targets.at(acc)->getSequence();
          ProteinReport pe;
          pe.sequence = protein_sequence;
          report.emplace(acc, pe); // add to report
          protein_length = protein_sequence.size();
        }
        else
        {
          protein_length = (int)it->second.sequence.size();
        }
        
        // TODO: potentially could save that loop. Doesn't handle peptides mapping twice into same protein (but different position)
        const PeptideEvidence& ph_evidence = *find_if(ph_evidences.begin(), ph_evidences.end(), [&acc] (const PeptideEvidence& e) 
          { return e.getProteinAccession() == acc; } );

        const int peptide_start_in_protein = ph_evidence.getStart();

        if (peptide_start_in_protein < 0) continue; // TODO: can this happen?
 
        if (best_localization >= 0) 
        { // XL was localized
          // calculate position in protein
          int xl_pos_in_protein = peptide_start_in_protein + best_localization;

          if (xl_pos_in_protein >= protein_length) continue; // TODO: can this happen?

          // create basic site information for this protein
          if (auto it = report[acc].aa_level_localization.find(xl_pos_in_protein); 
            it == report[acc].aa_level_localization.end())
          {
            auto& protein_sequence = acc2protein_targets.at(acc)->getSequence();            
            report[acc].aa_level_localization[xl_pos_in_protein].AA = protein_sequence[xl_pos_in_protein];
          }

          AALevelLocalization::LocalizedXL xl;
          xl.adduct = NA;
          xl.NT = NT;
          xl.charge = charge;          
          report[acc].aa_level_localization[xl_pos_in_protein].peptide2XL[peptide_sequence_string].push_back(xl);
        }
        else
        { // not localized? annotate region
          RegionLevelLocalization::UnlocalizedXL xl;
          xl.adduct = NA;
          xl.charge = charge;          
          int start = ph_evidence.getStart();
          int end = ph_evidence.getEnd();
          report[acc].region_level_localization[{start, end}].peptide2unlocalizedXL[peptide_sequence_string].push_back(xl);
        }
        
        if (is_unique)
        {
          report[acc].CSMs_of_unique_peptides++; // count CSM (localized and unlocalized) of unique peptides
        }
        else
        {
          report[acc].CSMs_of_shared_peptides++; // count CSM (localized and unlocalized)
        }
      }
    }    

    return report;
  }

  set<string> printXLSiteDetails(
    TextFile& tsv_file, 
    const std::string& accession, 
    size_t position, 
    const AALevelLocalization& aa_loc,
    map<string, vector<RegionLevelLocalization::UnlocalizedXL>>& peptides2unlocalizedXL,    
    std::map<std::string, std::set<std::string>>& peptide2proteins,
    map<string, map<string, set<pair<size_t, size_t>>>>& peptides2proteins2regions,
    map<string, set<string>>& protein2proteins)
  {
    set<string> printed_peptides;

    // one row per localized peptide
    const string line_start = accession + "\t" + aa_loc.AA + "\t" + String(position) + "\t";
    for (const auto& [peptide, localizedXLs] : aa_loc.peptide2XL)
    {
      printed_peptides.insert(peptide);      
      // protein, AA, position
      String l = line_start;

      // TODO: handle all set entries (e.g., if peptide maps multiple times in same protein)
      // start and end of peptide in protein
      string region = String(peptides2proteins2regions[peptide][accession].begin()->first) + "\t" + String(peptides2proteins2regions[peptide][accession].begin()->second) +"\t";
      l += region;

      bool is_unique = peptide2proteins[peptide].size() == 1;
      vector<RegionLevelLocalization::UnlocalizedXL>* unlocalized = nullptr;
      if (auto it = peptides2unlocalizedXL.find(peptide); it != peptides2unlocalizedXL.end())
      {
        unlocalized = &(it->second);
      }      

      // condense information down
      set<string> adduct_set, nt_set, charge_set;      
      set<string> unique_peptidoforms, shared_peptidoforms; // XLs that differ in either adduct, nt or charge

      size_t unique_localized_CSM_count{};
      size_t shared_localized_CSM_count{}; // shared
      for (const auto& xls : localizedXLs)
      {
        adduct_set.insert(xls.adduct);
        nt_set.insert(xls.NT);
        charge_set.insert(String(xls.charge));        
        if (is_unique)
        {
          unique_localized_CSM_count++;
          unique_peptidoforms.insert(xls.adduct + xls.NT + String(xls.charge));
        }
        else
        {
          shared_localized_CSM_count++;
          shared_peptidoforms.insert(xls.adduct + xls.NT + String(xls.charge));
        } 
      }

      // peptide was found but not localized
      set<string> unlocalized_adduct_set, unlocalized_charge_set;      
      set<string> unique_unlocalized_peptidoforms, shared_unlocalized_peptidoforms; // XLs that differ in either adduct, nt or charge
      size_t unique_unlocalized_CSM_count{};
      size_t shared_unlocalized_CSM_count{};
      if (unlocalized != nullptr)
      {
        for (const auto& xls : *unlocalized)
        {
          unlocalized_adduct_set.insert(xls.adduct);
          unlocalized_charge_set.insert(String(xls.charge));        
          if (is_unique)
          {
            unique_unlocalized_CSM_count++;
            unique_unlocalized_peptidoforms.insert(xls.adduct + String(xls.charge));
          }
          else
          {
            shared_unlocalized_CSM_count++;
            shared_unlocalized_peptidoforms.insert(xls.adduct + String(xls.charge));
          }     
        }
      }

      // print adducts, nucleotides and charge sets
      l += ListUtils::concatenate(adduct_set, ",") + "\t"; 
      l += ListUtils::concatenate(nt_set, ",") + "\t";
      l += ListUtils::concatenate(charge_set, ",") + "\t";
      l += String(unique_localized_CSM_count) + "\t";
      l += String(shared_localized_CSM_count) + "\t";
      l += String(unique_peptidoforms.size()) + "\t"; // peptide counts
      l += String(shared_peptidoforms.size()) + "\t";

      // print adducts, nucleotides and charge sets of unlocalized peptides
      l += ListUtils::concatenate(unlocalized_adduct_set, ",") + "\t"; 
      l += ListUtils::concatenate(unlocalized_charge_set, ",") + "\t";
      l += String(unique_unlocalized_CSM_count) + "\t";
      l += String(shared_unlocalized_CSM_count) + "\t";
      l += String(unique_unlocalized_peptidoforms.size()) + "\t"; // peptide counts
      l += String(shared_unlocalized_peptidoforms.size()) + "\t";

      // create string with other proteins
      auto ambiguities = peptide2proteins[peptide];
      protein2proteins[accession].insert(ambiguities.begin(), ambiguities.end()); // note: add same protein to group as well
      ambiguities.erase(accession);
      l += ListUtils::concatenate(ambiguities, ",") + "\t";

      // add peptide sequence and sequence level FDR (in the group of XLs)
      l += peptide + "\t";
      l += peptide_seq2XLFDR[peptide];
      tsv_file.addLine(l);
    }
    return printed_peptides;
  }

  void printXLRegionDetails(
    TextFile& tsv_file, 
    const std::string& accession,   
    const RegionLevelLocalization& region_loc,    
    const set<string>& remaining_peptides,    
    std::map<std::string, std::set<std::string>>& peptide2proteins,
    map<string, map<string, set<pair<size_t, size_t>>>>& peptides2proteins2regions,
    map<string, set<string>>& protein2proteins
    )
  {
    // one row per unlocalized peptide    
    for (const auto& [peptide, unlocalizedXLs] : region_loc.peptide2unlocalizedXL)
    {
      if (remaining_peptides.find(peptide) == remaining_peptides.end()) continue;

      // protein, AA, position
      String l = accession + "\t-\t-\t";

      // TODO: handle all set entries (e.g., if peptide maps multiple times in same protein)
      string region = String(peptides2proteins2regions[peptide][accession].begin()->first) + "\t" + String(peptides2proteins2regions[peptide][accession].begin()->second);
      l += region + "\t";

      bool is_unique = peptide2proteins[peptide].size() == 1;

      // peptide was found but not localized
      set<string> unlocalized_adduct_set, unlocalized_charge_set;      
      set<string> unique_unlocalized_peptidoforms, shared_unlocalized_peptidoforms; // XLs that differ in either adduct, nt or charge
      size_t unique_unlocalized_CSM_count{};
      size_t shared_unlocalized_CSM_count{};
      for (const auto& xls : unlocalizedXLs)
      {
        unlocalized_adduct_set.insert(xls.adduct);
        unlocalized_charge_set.insert(String(xls.charge));        
        if (is_unique)
        {
          unique_unlocalized_CSM_count++;
          unique_unlocalized_peptidoforms.insert(xls.adduct + String(xls.charge));
        }
        else
        {
          shared_unlocalized_CSM_count++;
          shared_unlocalized_peptidoforms.insert(xls.adduct + String(xls.charge));
        }  
      }

      // print adducts, nucleotides and charge sets
      l += "-\t-\t-\t0\t0\t0\t0\t"; 

      // print adducts, nucleotides and charge sets of unlocalized peptides
      l += ListUtils::concatenate(unlocalized_adduct_set, ",") + "\t"; 
      l += ListUtils::concatenate(unlocalized_charge_set, ",") + "\t";
      l += String(unique_unlocalized_CSM_count) + "\t";
      l += String(shared_unlocalized_CSM_count)+ "\t";
      l += String(unique_unlocalized_peptidoforms.size()) + "\t"; // peptide counts
      l += String(shared_unlocalized_peptidoforms.size()) + "\t";

      // create string with other proteins
      auto ambiguities = peptide2proteins[peptide];
      protein2proteins[accession].insert(ambiguities.begin(), ambiguities.end()); // note: add same protein to group as well
      ambiguities.erase(accession);
      l += ListUtils::concatenate(ambiguities, ",") + "\t";      
      l += peptide + "\t"; // add peptide sequence
      l += peptide_seq2XLFDR[peptide];
      tsv_file.addLine(l);
    }
  }

  // static 
  void  NuXLProteinReport::mapAccessionToTDProteins(ProteinIdentification& prot_id, std::map<String, ProteinHit*>& acc2protein_targets, std::map<String, ProteinHit*>& acc2protein_decoys)
  {
    std::vector<ProteinHit>& proteins = prot_id.getHits();
    for (ProteinHit& protein : proteins)
    {
      if (protein.getMetaValue("target_decoy").toString().hasPrefix("target"))
      {
        acc2protein_targets[protein.getAccession()] = &protein;
      }
      else
      {
        acc2protein_decoys[protein.getAccession()] = &protein;
      }
    }
  }

/*
  map<string, set<string>> acc2grouptype = calculateXLGroupType(const ProteinsReport& report, const ProteinIdentification& prot_id)
  {
    map<string, string> a2g;
    set<string> single_protein;
    
    map<set<string>, size_t> indist;
    size_t ind_index{};
    for (const auto& pg : prot_id.getIndistinguishableProteins())
    {
      set<string> s(pg.accessions.begin(), pg.accessions.end());
      if (auto it = indist.find(s); it == indist.end())
      {
        indist[s] = ind_index++;
      }      
    }

    for (const auto& [accession, pr] : report)
    {
      if (pr.CSMs_of_unique_peptides > 0)  // at least one unique peptide -> protein has unique evidence for XL
      {
        a2g.insert("single protein");
      }
      else if (indist.find(accession) != indist.end())
      {
        a2g.insert("ind. protein group");
      }
      else
      {
        a2g.insert("gen. protein group");
      }
    }

    return a2g;
  }
  */

  void NuXLProteinReport::annotateProteinModificationForTopHits(
    vector<ProteinIdentification>& prot_ids, 
    const vector<PeptideIdentification>& peps, 
    TextFile& tsv_file)
  {
    assert(prot_ids.size() == 1); // support for one run only

    // create lookup accession -> protein    
    ProteinIdentification& prot_id = prot_ids[0];

    // create lookup accession -> protein
    map<String, ProteinHit*> acc2protein_targets, acc2protein_decoys;
    NuXLProteinReport::mapAccessionToTDProteins(prot_id, acc2protein_targets, acc2protein_decoys);

    size_t CSMs_sum{}; // total number of XLed spectra

    // map peptide sequence to protein(s)
    std::map<std::string, std::set<std::string>> peptide2proteins;
    std::map<std::string, std::set<std::string>> protein2peptides;
    // lookup from peptide to its proteins and the region it maps to
    map<string, map<string, set<pair<size_t, size_t>>>> peptides2proteins2regions;
    for (const PeptideIdentification& pep : peps)
    {
      auto& hits = pep.getHits();
      if (hits.empty()) continue;
      const PeptideHit& ph = hits[0]; // only consider top hit
      auto peptide_sequence = ph.getSequence().toUnmodifiedString();
      const std::vector<PeptideEvidence>& ph_evidences = ph.getPeptideEvidences();
      ++CSMs_sum;
      for (auto& ph_evidence : ph_evidences)
      {
        const String& acc = ph_evidence.getProteinAccession();
        bool is_target = acc2protein_targets.find(acc) != acc2protein_targets.end();
        if (!is_target) continue; // skip decoys            
        peptide2proteins[peptide_sequence].insert(acc);
        protein2peptides[acc].insert(peptide_sequence);
        peptides2proteins2regions[peptide_sequence][acc].insert({ph_evidence.getStart(), ph_evidence.getEnd()});
      }
    }

    ProteinsReport r = getProteinReportEntries(peps, acc2protein_targets, peptide2proteins);

    // copy to vector so we can sort
    vector<pair<string, ProteinReport>> report;
    copy(r.begin(), r.end(), back_inserter(report));
    r.clear();

    // sort report entries (largest number of XL PSM count first) 
    cout << "Sorting entries... " << endl;
    std::sort(report.begin(), report.end(), 
      [](const pair<string, ProteinReport> & a, const pair<string, ProteinReport> & b) -> bool
      { 
         return std::tie(a.second.CSMs_of_unique_peptides, a.second.CSMs_of_shared_peptides, a.first) 
          > std::tie(b.second.CSMs_of_unique_peptides, b.second.CSMs_of_shared_peptides, b.first);
      }); 


    // write to file
    cout << "Writing " << report.size() << " proteins to tsv file... " << endl;

    tsv_file.addLine(String("accession\tAA\tpos.\tstart\tend\t") + 
                     "adducts (loc. + unique)\tNT (loc. + unique)\tcharges (loc. + unique)\t" + 
                     "CSMs (loc. + unique)\tCSMs (loc. + shared)\tprecursors (loc. + unique)\tprecursors (loc. + shared)\t" +
                     "adducts (\\wo loc. + unique)\tcharges (\\wo loc. + unique)\t" + 
                     "CSMs (\\wo loc. + unique)\tCSMs (\\wo loc. + shared)\tprecursors (\\wo loc. + unique)\tprecursors (\\wo loc. + shared)\t" +
                     "ambiguities\tpeptide\tq-value (peptide seq. level)"
      );

    map<string, set<string>> protein2proteins;

    for (const auto& [accession, pr] : report)
    {
      // lookup to determine if and where the given peptide 
      // was identified without localization and where it maps to in the protein.
      map<string, vector<RegionLevelLocalization::UnlocalizedXL>> peptides2unlocalizedXL;
      for (const auto& [start_end, region_detail] : pr.region_level_localization)
      {
        for (const auto& [peptide, xls] : region_detail.peptide2unlocalizedXL)
        {
          peptides2unlocalizedXL[peptide] = xls;
        }
      }  

      // first write lines with site localizations
      set<string> printed_peptides;
      for (const auto& [position, aa_loc] : pr.aa_level_localization)
      {
        set<string> p = printXLSiteDetails(
          tsv_file,
          accession,
          position,
          aa_loc,
          peptides2unlocalizedXL,
          peptide2proteins,
          peptides2proteins2regions,
          protein2proteins
          );
        printed_peptides.insert(p.begin(), p.end());
      }

      // determine peptides/regions not yet printed (e.g., no site localization exists for those)
      set<string> all_peptides = protein2peptides.at(accession);

      set<string> remaining_peptides;
      std::set_difference(all_peptides.begin(), all_peptides.end(), 
        printed_peptides.begin(), printed_peptides.end(),
        std::inserter(remaining_peptides, remaining_peptides.end()));

      // write remaining unlocalized peptides (=regions)
      for (const auto& [region, region_loc] : pr.region_level_localization)
      {
        printXLRegionDetails(
          tsv_file,
          accession,
          region_loc,
          remaining_peptides,
          peptide2proteins,
          peptides2proteins2regions,
          protein2proteins);        
      }
    }

    tsv_file.addLine("\n=============================================================");
    tsv_file.addLine("Run summary:");
    tsv_file.addLine("CSMs:\t" + String(CSMs_sum));
    tsv_file.addLine("Proteins:\t" + String(report.size()));

    tsv_file.addLine("\n=============================================================");
    tsv_file.addLine("Protein summary:");
    tsv_file.addLine("accession\tCSMs (unique pep.)\tCSMs (shared pep.)\tgroup type");

    // calculate indistinguishable groups
    vector<ProteinIdentification::ProteinGroup> ipg;
    if (!prot_id.getHits().empty())
    {
      vector<PeptideIdentification> pep_copy{peps}; // TODO: why copy needed?
      IDBoostGraph ibg{prot_id, pep_copy, true, false, false}; // only consider top hit
      ibg.calculateAndAnnotateIndistProteins(false); // only indistinguishable protein groups
      ipg = prot_id.getIndistinguishableProteins();
      std::sort(std::begin(ipg), std::end(ipg));
    }

    // calculate proteins with unique XL peptide
    set<string> accessionToUniquePeptides;
    for (const PeptideIdentification& pep : peps)
    {
      auto& hits = pep.getHits();
      if (hits.empty()) continue;
      const PeptideHit& ph = hits[0]; // only consider top hit      
      const AASequence& peptide_sequence = ph.getSequence();
      const std::string peptide_sequence_string = peptide_sequence.toUnmodifiedString();
      const std::set<std::string>& proteins = peptide2proteins.at(peptide_sequence_string);
      const bool is_unique = proteins.size() == 1;
      if (is_unique) 
      {
        accessionToUniquePeptides.insert(*proteins.begin());
      }
    }

    map<string, ProteinIdentification::ProteinGroup*> accessionToIndistinguishableGroup;
    for (auto& pg : ipg)
    {
      for (const auto& a : pg.accessions)
      {
        accessionToIndistinguishableGroup[a] = &pg;
      }
    }

    // print single proteins and ind. protein groups
    set<string> printed_ind_group;
    for (const auto& [accession, pr] : report)
    {
      if (accessionToUniquePeptides.count(accession) == 1)
      { // protein with unique peptide
        String group_type = "protein";
        tsv_file.addLine(accession + "\t" + String(pr.CSMs_of_unique_peptides) 
          + "\t" + String(pr.CSMs_of_shared_peptides) 
          + "\t" + group_type);        
      }
      else if (auto it = accessionToIndistinguishableGroup.find(accession);
        it != accessionToIndistinguishableGroup.end())
      { // ind. protein group
        if (printed_ind_group.count(accession) == 0)
        {
          const ProteinIdentification::ProteinGroup* pg = it->second;
          String a = ListUtils::concatenate(pg->accessions, ",");
          String group_type = "ind. protein group";
          tsv_file.addLine(a + "\t" + String(pr.CSMs_of_unique_peptides) 
            + "\t" + String(pr.CSMs_of_shared_peptides) 
            + "\t" + group_type);

          for (auto a : pg->accessions) 
          { // mark all ind. protein members as already printed
            printed_ind_group.insert(std::move(a));
          }
        }          
      }
      else
      { // general protein groups
        set<string> group_neighbors = protein2proteins[accession];
        group_neighbors.erase(accession);
        tsv_file.addLine(accession + "\t" + String(pr.CSMs_of_unique_peptides) 
          + "\t" + String(pr.CSMs_of_shared_peptides) 
          + "\tgen. protein group (shares XL peptides with: " + ListUtils::concatenate(group_neighbors, ",") + ")");
      }
    }

    tsv_file.addLine("\n=============================================================");
    tsv_file.addLine("Crosslink efficiency (AA freq. / AA freq. in all CSMs):");
    auto aa_xl_freq = getCrossLinkEfficiency(peps);
    for (auto& m : aa_xl_freq) 
    { 
      tsv_file.addLine(String(m.first) + "\t" + String(m.second));
    }

    tsv_file.addLine("\n=============================================================");
    tsv_file.addLine("Precursor adduct summary:");
    tsv_file.addLine("Precursor adduct:\tPSMs:\tPSMs(%)");

    map<String, size_t> adduct2count = countAdducts(peps);
    vector<pair<size_t, String>> count2adduct;
    size_t total_psms{};
    for (const auto& ac : adduct2count)
    {
      count2adduct.push_back({ac.second, ac.first});
      total_psms += ac.second;
    }

    std::sort(count2adduct.begin(), count2adduct.end(), 
      [](const pair<size_t, String> & a, const pair<size_t, String> & b) -> bool
      { 
         return std::tie(a.first, a.second) > std::tie(b.first, b.second);
      }); 

    for (const auto& ca : count2adduct)
    {
      tsv_file.addLine(ca.second + "\t" + String(ca.first) + "\t" + String(100.0 * (double)ca.first / (double)total_psms));
    }

  }

}
