// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzTabFile.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <QtCore/QString>

#include <boost/regex.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

// TODO fix all the shadowed "String s"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

namespace OpenMS
{

  MzTabFile::MzTabFile():
  store_protein_reliability_(false),
  store_peptide_reliability_(false),
  store_psm_reliability_(false),
  store_smallmolecule_reliability_(false),
  store_protein_uri_(false),
  store_peptide_uri_(false),
  store_psm_uri_(false),
  store_smallmolecule_uri_(false),
    store_protein_goterms_(false),
    store_nucleic_acid_reliability_(false),
    store_oligonucleotide_reliability_(false),
    store_osm_reliability_(false),
    store_nucleic_acid_uri_(false),
    store_oligonucleotide_uri_(false),
    store_osm_uri_(false),
    store_nucleic_acid_goterms_(false)
  {
  }

  MzTabFile::~MzTabFile() = default;

  std::pair<int, int> MzTabFile::extractIndexPairsFromBrackets_(const String & s)
  {
  std::pair<Int, Int> pair(0,0);

  // ^        # Match the start of the line
  // .*?      # Non-greedy match anything
  // \[       # Upto the first opening bracket (escaped)
  // (\d+)    # Match a digit string (one or more)
  // \]       # Match closing bracket
  // .*       # Match the rest of the line
  // $        # Match the end of the line
  boost::sregex_token_iterator end;

  boost::regex rx_first_number(R"(^.*?\[(\d+)\].*$)");
  boost::sregex_token_iterator it1(s.begin(), s.end(), rx_first_number, 1);
  if (it1 != end)
  {
    pair.first = String(*it1++).toInt();
  }

  boost::regex  rx_second_number(R"(^.*?\[\d+\].*?\[(\d+)\].*$)");
  boost::sregex_token_iterator it2(s.begin(), s.end(), rx_second_number, 1);
  if (it2 != end)
  {
    pair.second = String(*it2++).toInt();
  }

  return pair;
  }

  void MzTabFile::load(const String& filename, MzTab& mz_tab)
  {
  TextFile tf(filename, true);

  MzTabMetaData mz_tab_metadata;
  MzTabProteinSectionRows mz_tab_protein_section_data;
  MzTabPeptideSectionRows mz_tab_peptide_section_data;
  MzTabPSMSectionRows mz_tab_psm_section_data;
  MzTabSmallMoleculeSectionRows mz_tab_small_molecule_section_data;
  map<Size, String> comment_rows;
  vector<Size> empty_rows;

  map<String, Size> protein_custom_opt_columns;  // map column name to original column index
  map<String, Size> peptide_custom_opt_columns;
  map<String, Size> psm_custom_opt_columns;
  map<String, Size> smallmolecule_custom_opt_columns;

  Size count_study_variable_description = 0;
  Size count_ms_run_location = 0;

  //  protein section column information
  Size protein_accession_index = 0;
  Size protein_description_index = 0;
  Size protein_taxid_index = 0;
  Size protein_species_index = 0;
  Size protein_database_index = 0;
  Size protein_database_version_index = 0;
  Size protein_search_engine_index = 0;
  map<Size, Size> protein_best_search_engine_score_to_column_index;  // score number to column index
  map<Size, std::pair<Size, Size> > protein_column_index_to_score_runs_pair;   // map column of the protein section to a search_engine_score[index1]_ms_run[index] index pair
  Size protein_reliability_index = 0;
  map<Size, Size> protein_num_psms_ms_run_indices;
  map<Size, Size> protein_num_peptides_distinct_ms_run_indices;
  map<Size, Size> protein_num_peptides_unique_ms_run_indices;
  Size protein_ambiguity_members_index = 0;
  Size protein_modifications_index = 0;
  Size protein_uri_index = 0;
  Size protein_go_terms_index = 0;
  Size protein_coverage_index = 0;
  map<Size, Size> protein_abundance_assay_indices;
  map<Size, Size> protein_abundance_study_variable_to_column_indices;
  map<Size, Size> protein_abundance_stdev_study_variable_to_column_indices;
  map<Size, Size> protein_abundance_std_error_study_variable_to_column_indices;

  //  peptide section column information
  Size peptide_sequence_index = 0;
  Size peptide_accession_index = 0;
  Size peptide_unique_index = 0;
  Size peptide_database_index = 0;
  Size peptide_database_version_index = 0;
  Size peptide_search_engine_index = 0;
  map<Size, Size> peptide_best_search_engine_score_to_column_index;
  map<Size, std::pair<Size, Size> > peptide_column_index_to_score_runs_pair; // map column of the peptide section to a search_engine_score[index1]_ms_run[index] index pair
  Size peptide_reliability_index = 0;
  Size peptide_modifications_index = 0;
  Size peptide_retention_time_index = 0;
  Size peptide_retention_time_window_index = 0;
  Size peptide_charge_index = 0;
  Size peptide_mass_to_charge_index = 0;
  // Size peptide_uri_index = 0;
  Size peptide_spectra_ref_index = 0;
  map<Size, Size> peptide_abundance_assay_indices;
  map<Size, Size> peptide_abundance_study_variable_to_column_indices;
  map<Size, Size> peptide_abundance_study_variable_stdev_to_column_indices;
  map<Size, Size> peptide_abundance_study_variable_std_error_to_column_indices;

  // psm section column information
  Size psm_sequence_index = 0;
  Size psm_psm_id_index = 0;
  Size psm_accession_index = 0;
  Size psm_unique_index = 0;
  Size psm_database_index = 0;
  Size psm_database_version_index = 0;
  Size psm_search_engine_index = 0;
  map<Size, Size> psm_search_engine_score_to_column_index;
  Size psm_reliability_index = 0;
  Size psm_modifications_index = 0;
  Size psm_retention_time_index = 0;
  Size psm_charge_index = 0;
  Size psm_exp_mass_to_charge_index = 0;
  Size psm_calc_mass_to_charge_index = 0;
  // Size psm_uri_index = 0;
  Size psm_spectra_ref_index = 0;
  Size psm_pre_index = 0;
  Size psm_post_index = 0;
  Size psm_start_index = 0;
  Size psm_end_index = 0;

  // small molecule column information
  Size smallmolecule_identifier_index = 0;
  Size smallmolecule_chemical_formula_index = 0;
  Size smallmolecule_smiles_index = 0;
  Size smallmolecule_inchi_key_index = 0;
  Size smallmolecule_description_index = 0;
  Size smallmolecule_exp_mass_to_charge_index = 0;
  Size smallmolecule_calc_mass_to_charge_index = 0;
  Size smallmolecule_charge_index = 0;
  Size smallmolecule_retention_time_index = 0;
  Size smallmolecule_taxid_index = 0;
  Size smallmolecule_species_index = 0;
  Size smallmolecule_database_index = 0;
  Size smallmolecule_database_version_index = 0;
  Size smallmolecule_reliability_index = 0;
  Size smallmolecule_uri_index = 0;
  Size smallmolecule_spectra_ref_index = 0;
  Size smallmolecule_search_engine_index = 0;
  map<Size, Size> smallmolecule_best_search_engine_score_to_column_index;
  map<Size, std::pair<Size, Size> > smallmolecule_column_index_to_score_runs_pair; // map column of the small molecule section to a search_engine_score[index1]_ms_run[index] index pair
  Size smallmolecule_modifications_index = 0;
  map<Size, Size> smallmolecule_abundance_assay_indices;
  map<Size, Size> smallmolecule_abundance_study_variable_indices;
  map<Size, Size> smallmolecule_abundance_stdev_study_variable_indices;
  map<Size, Size> smallmolecule_abundance_std_error_study_variable_indices;

  // potentially mandatory meta values (depending on mzTab type, mode and sections that are present)
  set<String> mandatory_meta_values;

  // mzTab sections present in the file. Influences compulsoriness of meta-values.
  set<String> sections_present;

  Size count_protein_search_engine_score = 0;
  Size count_peptide_search_engine_score = 0;
  Size count_psm_search_engine_score = 0;
  Size count_smallmolecule_search_engine_score = 0;

  Size line_number = 0;
  for (TextFile::ConstIterator sit = tf.begin(); sit != tf.end(); ++sit, ++line_number)
  {
    //  std::cout << *sit << std::endl;
    String s = *sit;

    // skip empty lines or lines that are too short
    if (s.trim().size() < 3)
    {
      empty_rows.push_back(line_number); // preserve empty lines to map comments to correct position
      continue;
    }

    const String section = s.prefix(3);

    // discard comments
    if (section == "COM")
    {
      comment_rows[line_number] = s;
      continue;
    }

    StringList cells;
    s.split("\t", cells);

    if (cells.size() < 3)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "Error parsing MzTab line: " + String(s) + ". Did you forget to use tabulator as separator?");
    }

    // parse metadata section
    if (section == "MTD")
    {
      sections_present.insert("MTD");
      StringList meta_key_fields; // the "-" separated fields of the metavalue key
      cells[1].split("-", meta_key_fields);
      String meta_key = meta_key_fields[0];

      if (cells[1].hasPrefix("mzTab-version"))
      {
        mz_tab_metadata.mz_tab_version.fromCellString(cells[2]);
        mandatory_meta_values.insert("mzTab-version");
      }
      else if (cells[1].hasPrefix("mzTab-mode"))
      {
        mz_tab_metadata.mz_tab_mode.fromCellString(cells[2]);
        mandatory_meta_values.insert("mzTab-mode");
      }
      else if (cells[1].hasPrefix("mzTab-type"))
      {
        mz_tab_metadata.mz_tab_type.fromCellString(cells[2]);
        mandatory_meta_values.insert("mzTab-type");
      }
      else if (cells[1].hasPrefix("mzTab-ID"))
      {
        mz_tab_metadata.mz_tab_id.fromCellString(cells[2]);
      }
      else if (meta_key == "title")
      {
        mz_tab_metadata.title.set(cells[2]);
      }
      else if (meta_key == "description")
      {
        mz_tab_metadata.description.set(cells[2]);
        mandatory_meta_values.insert("description");
      }
      else if (meta_key.hasPrefix("sample_processing["))
      {
        Int n = meta_key.substitute("sample_processing[", "").substitute("]","").trim().toInt();
        MzTabParameterList pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.sample_processing[n] = pl;
      }
      else if (meta_key.hasPrefix("instrument[") && meta_key_fields[1] == "name")
      {
        Int n = meta_key_fields[0].substitute("instrument[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.instrument[n].name = p;
      }
      else if (meta_key.hasPrefix("instrument[") && meta_key_fields[1] == "source")
      {
        Int n = meta_key_fields[0].substitute("instrument[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.instrument[n].source = p;
      }
      else if (meta_key.hasPrefix("instrument[") && meta_key_fields.size() == 2 && meta_key_fields[1].hasPrefix("analyzer["))
      {
        Int n = meta_key_fields[0].substitute("instrument[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("analyzer[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.instrument[n].analyzer[m] = p;
      }
      else if (meta_key.hasPrefix("instrument[") && meta_key_fields[1] == "detector")
      {
        Int n = meta_key_fields[0].substitute("instrument[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.instrument[n].detector = p;
      }
      else if (meta_key.hasPrefix("software[") && meta_key_fields.size() == 1)
      {
        Int n = meta_key.substitute("software[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.software[n].software = p;
      }
      else if (meta_key.hasPrefix("software[") && meta_key_fields.size() == 2 && meta_key_fields[1].hasPrefix("setting["))
      {
        Int n = meta_key_fields[0].substitute("software[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("setting[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.software[n].setting[m] = p;
      }
      else if (meta_key.hasPrefix("protein_search_engine_score["))
      {
        Size n = (Size)meta_key_fields[0].substitute("protein_search_engine_score[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.protein_search_engine_score[n] = p;
        count_protein_search_engine_score = std::max(n, count_protein_search_engine_score); // will be checked to match number of entries in map to detect skipped entries or wrong numbering
      }
      else if (meta_key.hasPrefix("peptide_search_engine_score["))
      {
        Size n = (Size)meta_key_fields[0].substitute("peptide_search_engine_score[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.peptide_search_engine_score[n] = p;
        count_peptide_search_engine_score = std::max(n, count_peptide_search_engine_score); // will be checked to match number of entries in map to detect skipped entries or wrong numbering
      }
      else if (meta_key.hasPrefix("psm_search_engine_score["))
      {
        Size n = (Size)meta_key_fields[0].substitute("psm_search_engine_score[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.psm_search_engine_score[n] = p;
        count_psm_search_engine_score = std::max(n, count_psm_search_engine_score); // will be checked to match number of entries in map to detect skipped entries or wrong numbering
      }
      else if (meta_key.hasPrefix("smallmolecule_search_engine_score["))
      {
        Size n = (Size)meta_key_fields[0].substitute("smallmolecule_search_engine_score[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.smallmolecule_search_engine_score[n] = p;
        count_smallmolecule_search_engine_score = std::max(n, count_smallmolecule_search_engine_score); // will be checked to match number of entries in map to detect skipped entries or wrong numbering
      }
      else if (meta_key == "false_discovery_rate")
      {
        MzTabParameterList pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.false_discovery_rate = pl;
      }
      else if (meta_key.hasPrefix("publication["))
      {
        Int n = meta_key.substitute("publication[", "").substitute("]","").trim().toInt();
        MzTabString sl;
        sl.fromCellString(cells[2]);
        mz_tab_metadata.publication[n] = sl;
      }
      else if (meta_key.hasPrefix("contact") && meta_key_fields[1] == "name")
      {
        Int n = meta_key_fields[0].substitute("contact[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.contact[n].name = s;
      }
      else if (meta_key.hasPrefix("contact") && meta_key_fields[1] == "affiliation")
      {
        Int n = meta_key_fields[0].substitute("contact[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.contact[n].affiliation = s;
      }
      else if (meta_key.hasPrefix("contact") && meta_key_fields[1] == "email")
      {
        Int n = meta_key_fields[0].substitute("contact[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.contact[n].email = s;
      }
      else if (meta_key.hasPrefix("uri["))
      {
        Int n = meta_key_fields[0].substitute("uri[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.uri[n] = s;
      }
      //TODO: add mandatory check
      else if (meta_key.hasPrefix("variable_mod[") &&  meta_key_fields.size() == 1)
      {
        Int n = meta_key.substitute("variable_mod[", "").substitute("]","").trim().toInt();
        MzTabParameter pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.variable_mod[n].modification = pl;
      }
      else if (meta_key.hasPrefix("variable_mod[") &&  meta_key_fields[1] == "site") // variable_mod[1-n]-site
      {
        Int n = meta_key.substitute("variable_mod[", "").substitute("]","").trim().toInt();
        MzTabString pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.variable_mod[n].site = pl;
      }
      else if (meta_key.hasPrefix("variable_mod[") &&  meta_key_fields[1] == "position") // variable_mod[1-n]-position
      {
        Int n = meta_key.substitute("variable_mod[", "").substitute("]","").trim().toInt();
        MzTabString pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.variable_mod[n].position = pl;
      }
      //TODO: add mandatory check
      else if (meta_key.hasPrefix("fixed_mod[") &&  meta_key_fields.size() == 1) // fixed_mod[1-n]
      {
        Int n = meta_key.substitute("fixed_mod[", "").substitute("]","").trim().toInt();
        MzTabParameter pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.fixed_mod[n].modification = pl;
      }
      else if (meta_key.hasPrefix("fixed_mod[") &&  meta_key_fields[1] == "site") // fixed_mod[1-n]-site
      {
        Int n = meta_key.substitute("fixed_mod[", "").substitute("]","").trim().toInt();
        MzTabString pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.fixed_mod[n].site = pl;
      }
      else if (meta_key.hasPrefix("fixed_mod[") &&  meta_key_fields[1] == "position") // fixed_mod[1-n]-position
      {
        Int n = meta_key.substitute("fixed_mod[", "").substitute("]","").trim().toInt();
        MzTabString pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.fixed_mod[n].position = pl;
      }
      else if (meta_key == "quantification_method")
      {
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.quantification_method = p;
      }
      else if (meta_key == "protein" && meta_key_fields[1] == "quantification_unit")
      {
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.protein_quantification_unit = p;
        mandatory_meta_values.insert("protein-quantification_unit");
      }
      else if (meta_key == "peptide" && meta_key_fields[1] == "quantification_unit")
      {
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.peptide_quantification_unit = p;
        mandatory_meta_values.insert("peptide-quantification_unit");
      }
      else if (meta_key == "small_molecule" && meta_key_fields[1] == "quantification_unit")
      {
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.small_molecule_quantification_unit = p;
        mandatory_meta_values.insert("small_molecule-quantification_unit");
      }
      else if (meta_key.hasPrefix("ms_run[") && meta_key_fields[1] == "format")
      {
        Int n = meta_key_fields[0].substitute("ms_run[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.ms_run[n].format = p;
      }
      else if (meta_key.hasPrefix("ms_run[") && meta_key_fields[1] == "location")
      {
        Int n = meta_key_fields[0].substitute("ms_run[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.ms_run[n].location = p;
        count_ms_run_location = std::max((Size)n, (Size)count_ms_run_location); // will be checked to match number of entries in map to detect skipped entries or wrong numbering
      }
      else if (meta_key.hasPrefix("ms_run[") && meta_key_fields[1] == "id_format")
      {
        Int n = meta_key_fields[0].substitute("ms_run[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.ms_run[n].id_format = p;
      }
      else if (meta_key.hasPrefix("ms_run[") && meta_key_fields[1] == "fragmentation_method")
      {
        Int n = meta_key_fields[0].substitute("ms_run[", "").substitute("]","").trim().toInt();
        MzTabParameterList p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.ms_run[n].fragmentation_method = p;
      }
      else if (meta_key.hasPrefix("custom["))
      {
        Int n = meta_key.substitute("custom[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.custom[n] = p;
      }
      else if (meta_key.hasPrefix("sample[") && meta_key_fields[1].hasPrefix("species["))
      {
        Int n = meta_key_fields[0].substitute("sample[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("species[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.sample[n].species[m] = p;
      }
      else if (meta_key.hasPrefix("sample[") && meta_key_fields[1].hasPrefix("tissue["))
      {
        Int n = meta_key_fields[0].substitute("sample[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("tissue[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.sample[n].tissue[m] = p;
      }
      else if (meta_key.hasPrefix("sample[") && meta_key_fields[1].hasPrefix("cell_type["))
      {
        Int n = meta_key_fields[0].substitute("sample[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("cell_type[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.sample[n].cell_type[m] = p;
      }
      else if (meta_key.hasPrefix("sample[") && meta_key_fields[1].hasPrefix("disease["))
      {
        Int n = meta_key_fields[0].substitute("sample[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("disease[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.sample[n].disease[m] = p;
      }
      else if (meta_key.hasPrefix("sample[") && meta_key_fields[1] == "description")
      {
        Int n = meta_key_fields[0].substitute("sample[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.sample[n].description = p;
      }
      else if (meta_key.hasPrefix("sample[") && meta_key_fields[1].hasPrefix("custom["))
      {
        Int n = meta_key_fields[0].substitute("sample[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("custom[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.sample[n].custom[m] = p;
      }
      else if (meta_key.hasPrefix("assay[") && meta_key_fields[1] == "quantification_reagent")
      {
        Int n = meta_key_fields[0].substitute("assay[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.assay[n].quantification_reagent = p;
      }
      else if (meta_key.hasPrefix("assay[") && meta_key_fields[1].hasPrefix("quantification_mod[") && meta_key_fields.size() == 2) // assay[]-quantification_mod[]
      {
        Int n = meta_key_fields[0].substitute("assay[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("quantification_mod[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.assay[n].quantification_mod[m].modification = p;
      }
      else if (meta_key.hasPrefix("assay[") && meta_key_fields[1].hasPrefix("quantification_mod[") && meta_key_fields[2] == "site")
      {
        Int n = meta_key_fields[0].substitute("assay[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("quantification_mod[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.assay[n].quantification_mod[m].site = s;
      }
      else if (meta_key.hasPrefix("assay[") && meta_key_fields[1].hasPrefix("quantification_mod[") && meta_key_fields[2] == "position")
      {
        Int n = meta_key_fields[0].substitute("assay[", "").substitute("]","").trim().toInt();
        Int m = meta_key_fields[1].substitute("quantification_mod[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.assay[n].quantification_mod[m].position = s;
      }
      else if (meta_key.hasPrefix("assay[") && meta_key_fields[1] == "sample_ref")
      {
        Int n = meta_key_fields[0].substitute("assay[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.assay[n].sample_ref = p;
      }
      else if (meta_key.hasPrefix("cv[") && meta_key_fields[1] == "label")
      {
        Int n = meta_key_fields[0].substitute("cv[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.cv[n].label = p;
      }
      else if (meta_key.hasPrefix("cv[") && meta_key_fields[1] == "full_name")
      {
        Int n = meta_key_fields[0].substitute("cv[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.cv[n].full_name = p;
      }
      else if (meta_key.hasPrefix("cv[") && meta_key_fields[1] == "version")
      {
        Int n = meta_key_fields[0].substitute("cv[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.cv[n].version = p;
      }
      else if (meta_key.hasPrefix("cv[") && meta_key_fields[1] == "url")
      {
        Int n = meta_key_fields[0].substitute("cv[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.cv[n].url = p;
      }
      else if (meta_key.hasPrefix("assay[") && meta_key_fields[1] == "ms_run_ref")
      {
        Int n = meta_key_fields[0].substitute("assay[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        s.substitute("ms_run[","").substitute("]","");
        vector<String> ms_run;
        s.split(',', ms_run);
        for (auto& a : ms_run)
        {
          a.trim();
          mz_tab_metadata.assay[n].ms_run_ref.push_back(a.toInt());
        }
      }
      else if (meta_key.hasPrefix("study_variable[") && meta_key_fields[1] == "assay_refs")
      {
        Int n = meta_key_fields[0].substitute("study_variable[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        s.substitute("assay[","").substitute("]","");
        vector<String> assays;
        s.split(',', assays);
        for (auto& a : assays)
        {
          a.trim();
          mz_tab_metadata.study_variable[n].assay_refs.push_back(a.toInt());
        }
      }
      else if (meta_key.hasPrefix("study_variable[") && meta_key_fields[1] == "sample_refs")
      {
        Int n = meta_key_fields[0].substitute("study_variable[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        s.substitute("sample[","").substitute("]","");

        vector<String> assays;
        s.split(',', assays);
        for (auto& a : assays)
        {
          a.trim();
          mz_tab_metadata.study_variable[n].sample_refs.push_back(a.toInt());
      }
      }
      else if (meta_key.hasPrefix("study_variable[") && meta_key_fields[1] == "description")
      {
        Int n = meta_key_fields[0].substitute("study_variable[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.study_variable[n].description = p;
        count_study_variable_description = std::max((Size)n, count_study_variable_description);
      }
      else if (meta_key.hasPrefix("colunit") && meta_key_fields[1] == "protein")
      {
        Int n = meta_key_fields[0].substitute("colunit[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        mz_tab_metadata.colunit_protein[n] = s;
      }
      else if (meta_key.hasPrefix("colunit") && meta_key_fields[1] == "peptide")
      {
        Int n = meta_key_fields[0].substitute("colunit[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        mz_tab_metadata.colunit_peptide[n] = s;
      }
      else if (meta_key.hasPrefix("colunit") && meta_key_fields[1] == "psm")
      {
        Int n = meta_key_fields[0].substitute("colunit[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        mz_tab_metadata.colunit_psm[n] = s;
      }
      else if (meta_key.hasPrefix("colunit") && meta_key_fields[1] == "small_molecule")
      {
        Int n = meta_key_fields[0].substitute("colunit[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        mz_tab_metadata.colunit_small_molecule[n] = s;
      }
    }

    // parse protein header section
    if (section == "PRH")
    {
      for (Size i = 0; i != cells.size(); ++i)
      {
        if (cells[i] == "accession")
        {
          protein_accession_index = i;
        }
        else if (cells[i] == "description")
        {
          protein_description_index = i;
        }
        else if (cells[i] == "taxid")
        {
          protein_taxid_index = i;
        }
        else if (cells[i] == "species")
        {
          protein_species_index = i;
        }
        else if (cells[i] == "database")
        {
          protein_database_index = i;
        }
        else if (cells[i] == "database_version")
        {
          protein_database_version_index = i;
        }
        else if (cells[i] =="search_engine")
        {
          protein_search_engine_index = i;
        }
        else if (cells[i].hasPrefix("best_search_engine_score["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("best_search_engine_score[", "").substitute("]","").trim().toInt();
          protein_best_search_engine_score_to_column_index[n] = i;
        }
        else if (cells[i].hasPrefix("search_engine_score["))
        {
          std::pair<Size, Size> pair = extractIndexPairsFromBrackets_(cells[i]);
          protein_column_index_to_score_runs_pair[i] = pair;
        }
        else if (cells[i] == "reliability")
        {
          protein_reliability_index = i;
        }
        else if (cells[i].hasPrefix("num_psms_ms_run["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("num_psms_ms_run[", "").substitute("]","").trim().toInt();
          protein_num_psms_ms_run_indices[n] = i;
        }
        else if (cells[i].hasPrefix("num_peptides_distinct_ms_run["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("num_peptides_distinct_ms_run[", "").substitute("]","").trim().toInt();
          protein_num_peptides_distinct_ms_run_indices[n] = i;
        }
        else if (cells[i].hasPrefix("num_peptides_unique_ms_run["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("num_peptides_unique_ms_run[", "").substitute("]","").trim().toInt();
          protein_num_peptides_unique_ms_run_indices[n] = i;
        }
        else if (cells[i] == "ambiguity_members")
        {
          protein_ambiguity_members_index = i;
        }
        else if (cells[i] == "modifications")
        {
          protein_modifications_index = i;
        }
        else if (cells[i] == "uri")
        {
          protein_uri_index = i;
        }
        else if (cells[i] == "go_terms")
        {
          protein_go_terms_index = i;
        }
        else if (cells[i] == "protein_coverage")
        {
          protein_coverage_index = i;
        }
        else if (cells[i].hasPrefix("protein_abundance_assay["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("protein_abundance_assay[", "").substitute("]","").trim().toInt();
          protein_abundance_assay_indices[n] = i;
        }
        else if (cells[i].hasPrefix("protein_abundance_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("protein_abundance_study_variable[", "").substitute("]","").trim().toInt();
          protein_abundance_study_variable_to_column_indices[n] = i;
        }
        else if (cells[i].hasPrefix("protein_abundance_stdev_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("protein_abundance_stdev_study_variable[", "").substitute("]","").trim().toInt();
          protein_abundance_stdev_study_variable_to_column_indices[n] = i;
        }
        else if (cells[i].hasPrefix("protein_abundance_std_error_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("protein_abundance_std_error_study_variable[", "").substitute("]","").trim().toInt();
          protein_abundance_std_error_study_variable_to_column_indices[n] = i;
        }
        else if (cells[i].hasPrefix("opt_"))
        {
          protein_custom_opt_columns[cells[i]] = i;
        }
      }
      continue;
    }

    // parse protein section
    if (section == "PRT")
    {
      sections_present.insert("PRT");
      // check if all mandatory columns are present
      if (protein_accession_index == 0)
      {
        cout << "Error: mandatory protein accession column missing" << endl;
      }

      if (protein_description_index == 0)
      {
        cout << "Error: mandatory protein description column missing" << endl;
      }

      if (protein_taxid_index == 0)
      {
        cout << "Error: mandatory protein taxid column missing" << endl;
      }

      if (protein_species_index == 0)
      {
        cout << "Error: mandatory protein species column missing" << endl;
      }

      if (protein_database_index == 0)
      {
        cout << "Error: mandatory protein database column missing" << endl;
      }

      if (protein_database_version_index == 0)
      {
        cout << "Error: mandatory protein database_version column missing" << endl;
      }

      if (protein_search_engine_index == 0)
      {
        cout << "Error: mandatory protein search_engine column missing" << endl;
      }

      if (protein_best_search_engine_score_to_column_index.empty())
      {
        cout << "Error: mandatory protein best_search_engine_score[1-n] column missing" << endl;
      }

      if (mz_tab_metadata.mz_tab_mode.toCellString() == "Complete")
      {
        if (protein_column_index_to_score_runs_pair.size() != mz_tab_metadata.ms_run.size())
        {
          cout << "Error: mandatory protein search_engine_score_ms_run column(s) missing. Expected "
               << mz_tab_metadata.ms_run.size() << " ms_runs (meta data section) but only " << protein_column_index_to_score_runs_pair.size()
               << " ms_runs provide score columns (protein section)." << endl;
        }
      }

      if (protein_num_psms_ms_run_indices.empty() && mz_tab_metadata.mz_tab_mode.toCellString() == "Complete"
          && mz_tab_metadata.mz_tab_type.toCellString() == "Identification")
      {
        cout << "Error: mandatory protein num_psms_ms_run column(s) missing" << endl;
      }

      if (protein_num_peptides_distinct_ms_run_indices.empty() && mz_tab_metadata.mz_tab_mode.toCellString() == "Complete"
          && mz_tab_metadata.mz_tab_type.toCellString() == "Identification")
      {
        cout << "Error: mandatory protein num_peptides_distinct_ms_run column(s) missing" << endl;
      }

      if (protein_num_peptides_unique_ms_run_indices.empty() && mz_tab_metadata.mz_tab_mode.toCellString() == "Complete"
          && mz_tab_metadata.mz_tab_type.toCellString() == "Identification")
      {
        cout << "Error: mandatory protein num_peptides_unique_ms_run column(s) missing" << endl;
      }

      if (protein_ambiguity_members_index == 0)
      {
        cout << "Error: mandatory protein ambiguity_members column missing" << endl;
      }

      if (protein_modifications_index == 0)
      {
        cout << "Error: mandatory protein modifications column missing" << endl;
      }

      if (protein_coverage_index == 0)
      {
        cout << "Error: mandatory protein coverage column missing" << endl;
      }

      if (protein_abundance_assay_indices.empty()&& mz_tab_metadata.mz_tab_mode.toCellString() == "Complete"
          && mz_tab_metadata.mz_tab_type.toCellString() == "Quantification")
      {
        cout << "Error: mandatory protein protein_abundance_assay column(s) missing" << endl;
      }

      if (protein_abundance_study_variable_to_column_indices.empty() && mz_tab_metadata.mz_tab_type.toCellString() == "Quantification")
      {
        cout << "Error: mandatory protein_abundance_study_variable column(s) missing" << endl;
      }

      if (protein_abundance_stdev_study_variable_to_column_indices.empty() && mz_tab_metadata.mz_tab_type.toCellString() == "Quantification")
      {
        cout << "Error: mandatory protein_abundance_stdev_study_variable column(s) missing" << endl;
      }

      if (protein_abundance_std_error_study_variable_to_column_indices.empty() && mz_tab_metadata.mz_tab_type.toCellString() == "Quantification")
      {
        cout << "Error: mandatory protein_abundance_stderr_study_variable column(s) missing" << endl;
      }

      MzTabProteinSectionRow row;
      row.accession.fromCellString(cells[protein_accession_index]);
      row.description.fromCellString(cells[protein_description_index]);
      row.taxid.fromCellString(cells[protein_taxid_index]);
      row.species.fromCellString(cells[protein_species_index]);
      row.database.fromCellString(cells[protein_database_index]);
      row.database_version.fromCellString(cells[protein_database_version_index]);
      row.search_engine.fromCellString(cells[protein_search_engine_index]);

      for (map<Size, Size>::const_iterator it = protein_best_search_engine_score_to_column_index.begin(); it != protein_best_search_engine_score_to_column_index.end(); ++it)
      {
        row.best_search_engine_score[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, std::pair<Size, Size> >::const_iterator it = protein_column_index_to_score_runs_pair.begin(); it != protein_column_index_to_score_runs_pair.end(); ++it)
      {
        row.search_engine_score_ms_run[it->second.first][it->second.second].fromCellString(cells[it->first]);
      }

      if (protein_reliability_index != 0)
      {
        row.reliability.fromCellString(cells[protein_reliability_index]);
      }

      for (map<Size, Size>::const_iterator it = protein_num_psms_ms_run_indices.begin(); it != protein_num_psms_ms_run_indices.end(); ++it)
      {
        row.num_psms_ms_run[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = protein_num_peptides_distinct_ms_run_indices.begin(); it != protein_num_peptides_distinct_ms_run_indices.end(); ++it)
      {
        row.num_peptides_distinct_ms_run[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = protein_num_peptides_unique_ms_run_indices.begin(); it != protein_num_peptides_unique_ms_run_indices.end(); ++it)
      {
        row.num_peptides_unique_ms_run[it->first].fromCellString(cells[it->second]);
      }

      row.ambiguity_members.fromCellString(cells[protein_ambiguity_members_index]);
      row.modifications.fromCellString(cells[protein_modifications_index]);

      if (protein_uri_index != 0)
      {
        row.uri.fromCellString(cells[protein_uri_index]);
      }

      if (protein_go_terms_index != 0)
      {
        row.go_terms.fromCellString(cells[protein_go_terms_index]);
      }

      row.coverage.fromCellString(cells[protein_coverage_index]);

      // quantification data
      for (map<Size, Size>::const_iterator it = protein_abundance_assay_indices.begin(); it != protein_abundance_assay_indices.end(); ++it)
      {
        row.protein_abundance_assay[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = protein_abundance_study_variable_to_column_indices.begin(); it != protein_abundance_study_variable_to_column_indices.end(); ++it)
      {
        row.protein_abundance_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = protein_abundance_stdev_study_variable_to_column_indices.begin(); it != protein_abundance_stdev_study_variable_to_column_indices.end(); ++it)
      {
        row.protein_abundance_stdev_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = protein_abundance_std_error_study_variable_to_column_indices.begin(); it != protein_abundance_std_error_study_variable_to_column_indices.end(); ++it)
      {
        row.protein_abundance_std_error_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<String, Size>::const_iterator it = protein_custom_opt_columns.begin(); it != protein_custom_opt_columns.end(); ++it)
      {
        MzTabString s;
        s.fromCellString(cells[it->second]);
        MzTabOptionalColumnEntry e(it->first, s);
        row.opt_.push_back(e);
      }

      mz_tab_protein_section_data.push_back(row);
      continue;
    }

    // parse peptide header section
    if (section == "PEH")
    {
      for (Size i = 0; i != cells.size(); ++i)
      {
        if (cells[i] == "sequence")
        {
          peptide_sequence_index = i;
        }
        else if (cells[i] == "accession")
        {
          peptide_accession_index = i;
        }
        else if (cells[i] == "unique")
        {
          peptide_unique_index = i;
        }
        else if (cells[i] == "database")
        {
          peptide_database_index = i;
        }
        else if (cells[i] == "database_version")
        {
          peptide_database_version_index = i;
        }
        else if (cells[i] == "search_engine")
        {
          peptide_search_engine_index = i;
        }
        else if (cells[i].hasPrefix("best_search_engine_score["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("best_search_engine_score[", "").substitute("]","").trim().toInt();
          peptide_best_search_engine_score_to_column_index[n] = i;
        }
        else if (cells[i].hasPrefix("search_engine_score["))
        {
          std::pair<Size, Size> pair = extractIndexPairsFromBrackets_(cells[i].toQString());
          peptide_column_index_to_score_runs_pair[i] = pair;
        }
        else if (cells[i] == "reliability")
        {
          peptide_reliability_index = i;
        }
        else if (cells[i] == "modifications")
        {
          peptide_modifications_index = i;
        }
        else if (cells[i] == "retention_time")
        {
          peptide_retention_time_index = i;
        }
        else if (cells[i] == "retention_time_window")
        {
          peptide_retention_time_window_index = i;
        }
        else if (cells[i] == "charge")
        {
          peptide_charge_index = i;
        }
        else if (cells[i] == "mass_to_charge")
        {
          peptide_mass_to_charge_index = i;
        }
        else if (cells[i] == "uri")
        {
          protein_uri_index = i;
        }
        else if (cells[i] == "spectra_ref")
        {
          peptide_spectra_ref_index = i;
        }
        else if (cells[i].hasPrefix("peptide_abundance_assay["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("peptide_abundance_assay[", "").substitute("]","").trim().toInt();
          peptide_abundance_assay_indices[n] = i;
        }
        else if (cells[i].hasPrefix("peptide_abundance_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("peptide_abundance_study_variable[", "").substitute("]","").trim().toInt();
          peptide_abundance_study_variable_to_column_indices[n] = i;
        }
        else if (cells[i].hasPrefix("peptide_abundance_stdev_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("peptide_abundance_stdev_study_variable[", "").substitute("]","").trim().toInt();
          peptide_abundance_study_variable_stdev_to_column_indices[n] = i;
        }
        else if (cells[i].hasPrefix("peptide_abundance_std_error_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("peptide_abundance_std_error_study_variable[", "").substitute("]","").trim().toInt();
          peptide_abundance_study_variable_std_error_to_column_indices[n] = i;
        }
        else if (cells[i].hasPrefix("opt_"))
        {
          peptide_custom_opt_columns[cells[i]] = i;
        }
      }
      continue;
    }

    // parse peptide section
    if (section == "PEP")
    {
      sections_present.insert("PRT");
      MzTabPeptideSectionRow row;
      row.sequence.fromCellString(cells[peptide_sequence_index]);
      row.accession.fromCellString(cells[peptide_accession_index]);
      row.unique.fromCellString(cells[peptide_unique_index]);
      row.database.fromCellString(cells[peptide_database_index]);
      row.database_version.fromCellString(cells[peptide_database_version_index]);
      row.search_engine.fromCellString(cells[peptide_search_engine_index]);

      for (map<Size, Size>::const_iterator it = peptide_best_search_engine_score_to_column_index.begin(); it != peptide_best_search_engine_score_to_column_index.end(); ++it)
      {
        row.best_search_engine_score[it->first].fromCellString(cells[it->second]);
      }

      for (auto it = peptide_column_index_to_score_runs_pair.begin(); it != peptide_column_index_to_score_runs_pair.end(); ++it)
      {
        row.search_engine_score_ms_run[it->second.first][it->second.second].fromCellString(cells[it->first]);
      }

      if (peptide_reliability_index != 0)
      {
        row.reliability.fromCellString(cells[peptide_reliability_index]);
      }

      row.modifications.fromCellString(cells[peptide_modifications_index]);
      row.retention_time.fromCellString(cells[peptide_retention_time_index]);
      row.retention_time_window.fromCellString(cells[peptide_retention_time_window_index]);
      row.charge.fromCellString(cells[peptide_charge_index]);
      row.mass_to_charge.fromCellString(cells[peptide_mass_to_charge_index]);

      // if (peptide_uri_index != 0) // always false
      // {
      //   row.uri.fromCellString(cells[peptide_uri_index]);
      // }

      row.spectra_ref.fromCellString(cells[peptide_spectra_ref_index]);

      // quantification data
      for (map<Size, Size>::const_iterator it = peptide_abundance_assay_indices.begin(); it != peptide_abundance_assay_indices.end(); ++it)
      {
        row.peptide_abundance_assay[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = peptide_abundance_study_variable_to_column_indices.begin(); it != peptide_abundance_study_variable_to_column_indices.end(); ++it)
      {
        row.peptide_abundance_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = peptide_abundance_study_variable_stdev_to_column_indices.begin(); it != peptide_abundance_study_variable_stdev_to_column_indices.end(); ++it)
      {
        row.peptide_abundance_stdev_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = peptide_abundance_study_variable_std_error_to_column_indices.begin(); it != peptide_abundance_study_variable_std_error_to_column_indices.end(); ++it)
      {
        row.peptide_abundance_std_error_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<String, Size>::const_iterator it = peptide_custom_opt_columns.begin(); it != peptide_custom_opt_columns.end(); ++it)
      {
        MzTabString s;
        s.fromCellString(cells[it->second]);
        MzTabOptionalColumnEntry e(it->first, s);
        row.opt_.push_back(e);
      }

      mz_tab_peptide_section_data.push_back(row);
      continue;
    }

    // parse PSM header section
    if (section == "PSH")
    {
      for (Size i = 0; i != cells.size(); ++i)
      {
        if (cells[i] == "sequence")
        {
          psm_sequence_index = i;
        }
        else if (cells[i] == "PSM_ID")
        {
          psm_psm_id_index = i;
        }
        else if (cells[i] == "accession")
        {
          psm_accession_index = i;
        }
        else if (cells[i] == "unique")
        {
          psm_unique_index = i;
        }
        else if (cells[i] == "database")
        {
          psm_database_index = i;
        }
        else if (cells[i] == "database_version")
        {
          psm_database_version_index = i;
        }
        else if (cells[i] == "search_engine")
        {
          psm_search_engine_index = i;
        }
        else if (cells[i].hasPrefix("search_engine_score["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("search_engine_score[", "").substitute("]","").trim().toInt();
          psm_search_engine_score_to_column_index[n] = i;
        }
        else if (cells[i].hasPrefix("reliability"))
        {
          psm_reliability_index = i;
        }
        else if (cells[i] == "modifications")
        {
          psm_modifications_index = i;
        }
        else if (cells[i] == "retention_time")
        {
          psm_retention_time_index = i;
        }
        else if (cells[i] == "charge")
        {
          psm_charge_index = i;
        }
        else if (cells[i] == "exp_mass_to_charge")
        {
          psm_exp_mass_to_charge_index = i;
        }
        else if (cells[i] == "calc_mass_to_charge")
        {
          psm_calc_mass_to_charge_index = i;
        }
        else if (cells[i] == "uri")
        {
          protein_uri_index = i;
        }
        else if (cells[i] == "spectra_ref")
        {
          psm_spectra_ref_index = i;
        }
        else if (cells[i] == "pre")
        {
          psm_pre_index = i;
        }
        else if (cells[i] == "post")
        {
          psm_post_index = i;
        }
        else if (cells[i] == "start")
        {
          psm_start_index = i;
        }
        else if (cells[i] == "end")
        {
          psm_end_index = i;
        }
        else if (cells[i] == "opt_")
        {
          psm_custom_opt_columns[cells[i]] = i;
        }
      }
      continue;
    }

    // parse peptide section
    if (section == "PSM")
    {
      MzTabPSMSectionRow row;
      row.sequence.fromCellString(cells[psm_sequence_index]);
      row.PSM_ID.fromCellString(cells[psm_psm_id_index]);
      row.accession.fromCellString(cells[psm_accession_index]);
      row.unique.fromCellString(cells[psm_unique_index]);
      row.database.fromCellString(cells[psm_database_index]);
      row.database_version.fromCellString(cells[psm_database_version_index]);
      row.search_engine.fromCellString(cells[psm_search_engine_index]);

      for (map<Size, Size>::const_iterator it = psm_search_engine_score_to_column_index.begin(); it != psm_search_engine_score_to_column_index.end(); ++it)
      {
        row.search_engine_score[it->first].fromCellString(cells[it->second]);
      }

      if (psm_reliability_index != 0)
      {
        row.reliability.fromCellString(cells[psm_reliability_index]);
      }

      row.modifications.fromCellString(cells[psm_modifications_index]);
      row.retention_time.fromCellString(cells[psm_retention_time_index]);
      row.charge.fromCellString(cells[psm_charge_index]);
      row.exp_mass_to_charge.fromCellString(cells[psm_exp_mass_to_charge_index]);
      row.calc_mass_to_charge.fromCellString(cells[psm_calc_mass_to_charge_index]);

      // always false
      // if (psm_uri_index != 0)
      // {
      //   row.uri.fromCellString(cells[psm_uri_index]);
      // }

      row.spectra_ref.fromCellString(cells[psm_spectra_ref_index]);
      row.pre.fromCellString(cells[psm_pre_index]);
      row.post.fromCellString(cells[psm_post_index]);
      row.start.fromCellString(cells[psm_start_index]);
      row.end.fromCellString(cells[psm_end_index]);

      for (map<String, Size>::const_iterator it = psm_custom_opt_columns.begin(); it != psm_custom_opt_columns.end(); ++it)
      {
        MzTabString s;
        s.fromCellString(cells[it->second]);
        MzTabOptionalColumnEntry e(it->first, s);
        row.opt_.push_back(e);
      }

      mz_tab_psm_section_data.push_back(row);
      continue;
    }

    // parse small molecule header section
    if (section == "SMH")
    {
      for (Size i = 0; i != cells.size(); ++i)
      {
        if (cells[i] == "identifier")
        {
          smallmolecule_identifier_index = i;
        }
        else if (cells[i] == "chemical_formula")
        {
          smallmolecule_chemical_formula_index = i;
        }
        else if (cells[i] == "smiles")
        {
          smallmolecule_smiles_index = i;
        }
        else if (cells[i] == "inchi_key")
        {
          smallmolecule_inchi_key_index = i;
        }
        else if (cells[i] == "description")
        {
          smallmolecule_description_index = i;
        }
        else if (cells[i] == "exp_mass_to_charge")
        {
          smallmolecule_exp_mass_to_charge_index = i;
        }
        else if (cells[i] == "calc_mass_to_charge")
        {
          smallmolecule_calc_mass_to_charge_index = i;
        }
        else if (cells[i] == "charge")
        {
          smallmolecule_charge_index = i;
        }
        else if (cells[i] == "retention_time")
        {
          smallmolecule_retention_time_index = i;
        }
        else if (cells[i] == "taxid")
        {
          smallmolecule_taxid_index = i;
        }
        else if (cells[i] == "species")
        {
          smallmolecule_species_index = i;
        }
        else if (cells[i] == "database")
        {
          smallmolecule_database_index = i;
        }
        else if (cells[i] == "database_version")
        {
          smallmolecule_database_version_index = i;
        }
        else if (cells[i] == "reliability")
        {
          smallmolecule_reliability_index = i;
        }
        else if (cells[i] == "uri")
        {
          smallmolecule_uri_index = i;
        }
        else if (cells[i] == "spectra_ref")
        {
          smallmolecule_spectra_ref_index = i;
        }
        else if (cells[i] == "search_engine")
        {
          smallmolecule_search_engine_index = i;
        }
        else if (cells[i].hasPrefix("best_search_engine_score["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("best_search_engine_score[", "").substitute("]","").trim().toInt();
          smallmolecule_best_search_engine_score_to_column_index[n] = i;
        }
        else if (cells[i].hasPrefix("search_engine_score["))
        {
          std::pair<Size, Size> pair = extractIndexPairsFromBrackets_(cells[i].toQString());
          smallmolecule_column_index_to_score_runs_pair[i] = pair;
        }
        else if (cells[i] == "modifications")
        {
          smallmolecule_modifications_index = i;
        }
        else if (cells[i].hasPrefix("smallmolecule_abundance_assay["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("smallmolecule_abundance_assay[", "").substitute("]","").trim().toInt();
          smallmolecule_abundance_assay_indices[n] = i;
        }
        else if (cells[i].hasPrefix("smallmolecule_abundance_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("smallmolecule_abundance_study_variable[", "").substitute("]","").trim().toInt();
          smallmolecule_abundance_study_variable_indices[n] = i;
        }
        else if (cells[i].hasPrefix("smallmolecule_abundance_stdev_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("smallmolecule_abundance_stdev_study_variable[", "").substitute("]","").trim().toInt();
          smallmolecule_abundance_stdev_study_variable_indices[n] = i;
        }
        else if (cells[i].hasPrefix("smallmolecule_abundance_std_error_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("smallmolecule_abundance_std_error_study_variable[", "").substitute("]","").trim().toInt();
          smallmolecule_abundance_std_error_study_variable_indices[n] = i;
        }
        else if (cells[i].hasPrefix("opt_"))
        {
          smallmolecule_custom_opt_columns[cells[i]] = i;
        }
      }
      continue;
    }

    // parse small molecule section
    if (section == "SML")
    {
      sections_present.insert("SML");
      MzTabSmallMoleculeSectionRow row;
      row.identifier.fromCellString(cells[smallmolecule_identifier_index]);
      row.chemical_formula.fromCellString(cells[smallmolecule_chemical_formula_index]);
      row.smiles.fromCellString(cells[smallmolecule_smiles_index]);
      row.inchi_key.fromCellString(cells[smallmolecule_inchi_key_index]);
      row.description.fromCellString(cells[smallmolecule_description_index]);
      row.exp_mass_to_charge.fromCellString(cells[smallmolecule_exp_mass_to_charge_index]);
      row.calc_mass_to_charge.fromCellString(cells[smallmolecule_calc_mass_to_charge_index]);
      row.charge.fromCellString(cells[smallmolecule_charge_index]);
      row.retention_time.fromCellString(cells[smallmolecule_retention_time_index]);
      row.taxid.fromCellString(cells[smallmolecule_taxid_index]);
      row.species.fromCellString(cells[smallmolecule_species_index]);
      row.database.fromCellString(cells[smallmolecule_database_index]);
      row.database_version.fromCellString(cells[smallmolecule_database_version_index]);

      if (smallmolecule_reliability_index != 0)
      {
        row.reliability.fromCellString(cells[smallmolecule_reliability_index]);
      }

      if (smallmolecule_uri_index != 0)
      {
        row.uri.fromCellString(cells[smallmolecule_uri_index]);
      }

      row.spectra_ref.fromCellString(cells[smallmolecule_spectra_ref_index]);

      row.search_engine.fromCellString(cells[smallmolecule_search_engine_index]);

      for (map<Size, Size>::const_iterator it = smallmolecule_best_search_engine_score_to_column_index.begin(); it != smallmolecule_best_search_engine_score_to_column_index.end(); ++it)
      {
        row.best_search_engine_score[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, std::pair<Size, Size> >::const_iterator it = smallmolecule_column_index_to_score_runs_pair.begin(); it != smallmolecule_column_index_to_score_runs_pair.end(); ++it)
      {
        row.search_engine_score_ms_run[it->second.first][it->second.second].fromCellString(cells[it->first]);
      }

      row.modifications.fromCellString(cells[smallmolecule_modifications_index]);

      // quantification data
      for (map<Size, Size>::const_iterator it = smallmolecule_abundance_assay_indices.begin(); it != smallmolecule_abundance_assay_indices.end(); ++it)
      {
        row.smallmolecule_abundance_assay[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = smallmolecule_abundance_study_variable_indices.begin(); it != smallmolecule_abundance_study_variable_indices.end(); ++it)
      {
        row.smallmolecule_abundance_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = smallmolecule_abundance_stdev_study_variable_indices.begin(); it != smallmolecule_abundance_stdev_study_variable_indices.end(); ++it)
      {
        row.smallmolecule_abundance_stdev_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = smallmolecule_abundance_std_error_study_variable_indices.begin(); it != smallmolecule_abundance_std_error_study_variable_indices.end(); ++it)
      {
        row.smallmolecule_abundance_std_error_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<String, Size>::const_iterator it = smallmolecule_custom_opt_columns.begin(); it != smallmolecule_custom_opt_columns.end(); ++it)
      {
        MzTabString s;
        s.fromCellString(cells[it->second]);
        MzTabOptionalColumnEntry e(it->first, s);
        row.opt_.push_back(e);
      }

      mz_tab_small_molecule_section_data.push_back(row);
      continue;
    }
  }

  // TODO: check compulsoriness
  //hasMandatoryMetaDataKeys_(mandatory_meta_values, sections_present, mz_tab_metadata);

  mz_tab.setMetaData(mz_tab_metadata);
  mz_tab.setProteinSectionRows(mz_tab_protein_section_data);
  mz_tab.setPeptideSectionRows(mz_tab_peptide_section_data);
  mz_tab.setPSMSectionRows(mz_tab_psm_section_data);
  mz_tab.setSmallMoleculeSectionRows(mz_tab_small_molecule_section_data);
  mz_tab.setEmptyRows(empty_rows);
  mz_tab.setCommentRows(comment_rows);
  }

  void MzTabFile::generateMzTabMetaDataSection_(const MzTabMetaData& md, StringList& sl) const
  {
  sl.push_back(String("MTD\tmzTab-version\t") + md.mz_tab_version.toCellString());
  sl.push_back(String("MTD\tmzTab-mode\t") + md.mz_tab_mode.toCellString());
  sl.push_back(String("MTD\tmzTab-type\t") + md.mz_tab_type.toCellString());

  if (!md.title.isNull())
  {
    String s = String("MTD\ttitle\t") + md.title.toCellString();
    sl.push_back(s);
  }

  if (!md.mz_tab_id.isNull())
  {
    String s = String("MTD\tmzTab-ID\t") + md.mz_tab_id.toCellString();
    sl.push_back(s);
  }

  sl.push_back(String("MTD\tdescription\t") + md.description.toCellString());

  for (map<Size, MzTabParameterList>::const_iterator it = md.sample_processing.begin(); it != md.sample_processing.end(); ++it)
  {
    String s = "MTD\tsample_processing[" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.protein_search_engine_score.begin(); it != md.protein_search_engine_score.end(); ++it)
  {
    String s = "MTD\tprotein_search_engine_score[" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.peptide_search_engine_score.begin(); it != md.peptide_search_engine_score.end(); ++it)
  {
    String s = "MTD\tpeptide_search_engine_score[" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.psm_search_engine_score.begin(); it != md.psm_search_engine_score.end(); ++it)
  {
    String s = "MTD\tpsm_search_engine_score[" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.smallmolecule_search_engine_score.begin(); it != md.smallmolecule_search_engine_score.end(); ++it)
  {
    String s = "MTD\tsmallmolecule_search_engine_score[" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
  }

    for (map<Size, MzTabParameter>::const_iterator it = md.nucleic_acid_search_engine_score.begin(); it != md.nucleic_acid_search_engine_score.end(); ++it)
    {
      String s = "MTD\tnucleic_acid_search_engine_score[" + String(it->first) + "]\t" + it->second.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabParameter>::const_iterator it = md.oligonucleotide_search_engine_score.begin(); it != md.oligonucleotide_search_engine_score.end(); ++it)
    {
      String s = "MTD\toligonucleotide_search_engine_score[" + String(it->first) + "]\t" + it->second.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabParameter>::const_iterator it = md.osm_search_engine_score.begin(); it != md.osm_search_engine_score.end(); ++it)
    {
      String s = "MTD\tosm_search_engine_score[" + String(it->first) + "]\t" + it->second.toCellString();
      sl.push_back(s);
    }

  for (map<Size, MzTabInstrumentMetaData>::const_iterator it = md.instrument.begin(); it != md.instrument.end(); ++it)
  {
    const MzTabInstrumentMetaData & imd = it->second;

    if (!imd.name.isNull())
    {
      String s = "MTD\tinstrument[" + String(it->first) + "]-name\t" +  imd.name.toCellString();
      sl.push_back(s);
    }

    if (!imd.source.isNull())
    {
      String s = "MTD\tinstrument[" + String(it->first) + "]-source\t" + imd.source.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabParameter>::const_iterator mit = imd.analyzer.begin(); mit != imd.analyzer.end(); ++mit)
    {
      if (!mit->second.isNull())
      {
        String s = "MTD\tinstrument[" + String(it->first) + "]-analyzer[" + String(mit->first) + "]\t" + mit->second.toCellString();
        sl.push_back(s);
      }
    }

    if (!imd.detector.isNull())
    {
      String s = "MTD\tinstrument[" + String(it->first) + "]-detector\t" + imd.detector.toCellString();
      sl.push_back(s);
    }
  }

  for (map<Size, MzTabSoftwareMetaData>::const_iterator it = md.software.begin(); it != md.software.end(); ++it)
  {
    String s = "MTD\tsoftware[" + String(it->first) + "]\t" + it->second.software.toCellString();
    sl.push_back(s);

    for (map<Size, MzTabString>::const_iterator jt = it->second.setting.begin(); jt != it->second.setting.end(); ++jt)
    {
      String s = "MTD\tsoftware[" + String(it->first) + "]-setting[" + String(jt->first) + String("]\t") + jt->second.toCellString();
      sl.push_back(s);
    }
  }

  if (!md.false_discovery_rate.isNull())
  {
    String s = "MTD\tfalse_discovery_rate\t" + md.false_discovery_rate.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.publication.begin(); it != md.publication.end(); ++it)
  {
    String s = "MTD\tpublication[" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabContactMetaData>::const_iterator it = md.contact.begin(); it != md.contact.end(); ++it)
  {
    const MzTabContactMetaData & mdc = it->second;
    if (!mdc.name.isNull())
    {
      String s = "MTD\tcontact[" + String(it->first) + "]-name\t" + mdc.name.toCellString();
      sl.push_back(s);
    }

    if (!mdc.affiliation.isNull())
    {
      String s = "MTD\tcontact[" + String(it->first) + "]-affiliation\t" + mdc.affiliation.toCellString();
      sl.push_back(s);
    }

    if (!mdc.email.isNull())
    {
      String s = "MTD\tcontact[" + String(it->first) + "]-email\t" + mdc.email.toCellString();
      sl.push_back(s);
    }
  }


  for (map<Size, MzTabString>::const_iterator it = md.uri.begin(); it != md.uri.end(); ++it)
  {
    String s = "MTD\turi[" + String(it->first) + String("]\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabModificationMetaData>::const_iterator it = md.fixed_mod.begin(); it != md.fixed_mod.end(); ++it)
  {
    const MzTabModificationMetaData & mod_md = it->second;
    if (!mod_md.modification.isNull())
    {
      String s = "MTD\tfixed_mod[" + String(it->first) + String("]\t")+ mod_md.modification.toCellString();
      sl.push_back(s);
    }
    else
    {
      //TODO: add CV for no fixed modification searched when it is available
    }

    if (!mod_md.site.isNull())
    {
      String s = "MTD\tfixed_mod[" + String(it->first) + String("]-site\t") + mod_md.site.toCellString();
      sl.push_back(s);
    }

    if (!mod_md.position.isNull())
    {
      String s = "MTD\tfixed_mod[" + String(it->first) + String("]-position\t") + mod_md.position.toCellString();
      sl.push_back(s);
    }
  }

  for (map<Size, MzTabModificationMetaData>::const_iterator it = md.variable_mod.begin(); it != md.variable_mod.end(); ++it)
  {
    const MzTabModificationMetaData & mod_md = it->second;
    if (!mod_md.modification.isNull())
    {
      String s = "MTD\tvariable_mod[" + String(it->first) + String("]\t") + mod_md.modification.toCellString();
      sl.push_back(s);
    }
    else
    {
      //TODO: add CV for no variable modification searched when it is available
    }

    if (!mod_md.site.isNull())
    {
      String s = "MTD\tvariable_mod[" + String(it->first) + String("]-site\t") + mod_md.site.toCellString();
      sl.push_back(s);
    }

    if (!mod_md.position.isNull())
    {
      String s = "MTD\tvariable_mod[" + String(it->first) + String("]-position\t")+ mod_md.position.toCellString();
      sl.push_back(s);
    }
  }

  // quantification_method
  if (!md.quantification_method.isNull())
  {
    String s = "MTD\tquantification_method\t" + md.quantification_method.toCellString();
    sl.push_back(s);
  }

  // protein-quantification_unit
  if (!md.protein_quantification_unit.isNull())
  {
    String s = "MTD\tprotein-quantification_unit\t" + md.protein_quantification_unit.toCellString();
    sl.push_back(s);
  }

  // peptide-quantification_unit
  if (!md.peptide_quantification_unit.isNull())
  {
    String s = "MTD\tpeptide-quantification_unit\t" + md.peptide_quantification_unit.toCellString();
    sl.push_back(s);
  }

  // small_molecule-quantification_unit
  if (!md.small_molecule_quantification_unit.isNull())
  {
    String s = "MTD\tsmall_molecule-quantification_unit\t" + md.small_molecule_quantification_unit.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabMSRunMetaData>::const_iterator it = md.ms_run.begin(); it != md.ms_run.end(); ++it)
  {
    const MzTabMSRunMetaData & msmd = it->second;

    if (!msmd.format.isNull())
    {
      String s = "MTD\tms_run[" + String(it->first) + "]-format\t" + msmd.format.toCellString();
      sl.push_back(s);
    }

    if (!msmd.location.isNull())
    {
      String s = "MTD\tms_run[" + String(it->first) + "]-location\t" + msmd.location.toCellString();
      sl.push_back(s);
    }

    if (!msmd.id_format.isNull())
    {
      String s = "MTD\tms_run[" + String(it->first) + "]-id_format\t" + msmd.id_format.toCellString();
      sl.push_back(s);
    }

    if (!msmd.fragmentation_method.isNull())
    {
      String s = "MTD\tms_run[" + String(it->first) + "]-fragmentation_method\t" + msmd.fragmentation_method.toCellString();
      sl.push_back(s);
    }
  }

  // custom
  for (map<Size, MzTabParameter>::const_iterator it = md.custom.begin(); it != md.custom.end(); ++it)
  {
    String s = "MTD\tcustom[" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabSampleMetaData>::const_iterator it = md.sample.begin(); it != md.sample.end(); ++it)
  {
    for (map<Size, MzTabParameter>::const_iterator sit = it->second.species.begin(); sit != it->second.species.end(); ++sit)
    {
      String s = "MTD\tsample[" + String(it->first) + "]-species[" + String(sit->first) + "]\t" + sit->second.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabParameter>::const_iterator sit = it->second.tissue.begin(); sit != it->second.tissue.end(); ++sit)
    {
      String s = "MTD\tsample[" + String(it->first) + "]-tissue[" + String(sit->first) + "]\t" + sit->second.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabParameter>::const_iterator sit = it->second.cell_type.begin(); sit != it->second.cell_type.end(); ++sit)
    {
      String s = "MTD\tsample[" + String(it->first) + "]-cell_type[" + String(sit->first) + "]\t" + sit->second.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabParameter>::const_iterator sit = it->second.disease.begin(); sit != it->second.disease.end(); ++sit)
    {
      String s = "MTD\tsample[" + String(it->first) + "]-disease[" + String(sit->first) + "]\t" + sit->second.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabParameter>::const_iterator sit = it->second.custom.begin(); sit != it->second.custom.end(); ++sit)
    {
      String s = "MTD\tsample[" + String(it->first) + "]-custom[" + String(sit->first) + "]\t" + sit->second.toCellString();
      sl.push_back(s);
    }

    if (!it->second.description.isNull())
    {
      String s = "MTD\tsample[" + String(it->first) + String("]-description\t") + it->second.description.toCellString();
      sl.push_back(s);
    }
  }

  for (map<Size, MzTabAssayMetaData>::const_iterator it = md.assay.begin(); it != md.assay.end(); ++it)
  {
    const MzTabAssayMetaData & amd = it->second;
    if (!amd.quantification_reagent.isNull())
    {
      String s = "MTD\tassay[" + String(it->first) + "]-quantification_reagent\t" + amd.quantification_reagent.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabModificationMetaData>::const_iterator mit = amd.quantification_mod.begin(); mit != amd.quantification_mod.end(); ++mit)
    {
      const MzTabModificationMetaData & mod = mit->second;
      if (!mod.modification.isNull())
      {
        String s = "MTD\tassay[" + String(it->first) + String("]-quantification_mod[") + String(mit->first) + String("]\t") + mod.modification.toCellString();
        sl.push_back(s);
      }

      if (!mod.site.isNull())
      {
        String s = "MTD\tassay[" + String(it->first) + String("]-quantification_mod[") + String(mit->first) + String("]-site\t") + mod.site.toCellString();
        sl.push_back(s);
      }

      if (!mod.position.isNull())
      {
        String s = "MTD\tassay[" + String(it->first) + String("]-quantification_mod[") + String(mit->first) + String("]-position\t") + mod.position.toCellString();
        sl.push_back(s);
      }
    }

    if (!amd.sample_ref.isNull())
    {
      String s = "MTD\tassay[" + String(it->first) + String("]-sample_ref\t") + amd.sample_ref.toCellString();
      sl.push_back(s);
    }

    if (!amd.ms_run_ref.empty())
    {
      String s = "MTD\tassay[" + String(it->first) + "]-ms_run_ref\t";
      bool first(true);
      for (auto const & a : amd.ms_run_ref)
      {
        if (!first) { s += ","; } else { first = false; }
        s += "ms_run[" + String(a) + "]";
      }
    sl.push_back(s);
  }
  }


  for (map<Size, MzTabStudyVariableMetaData>::const_iterator it = md.study_variable.begin(); it != md.study_variable.end(); ++it)
  {
    const MzTabStudyVariableMetaData & smd = it->second;

    if (!smd.assay_refs.empty())
    {
      String s = "MTD\tstudy_variable[" + String(it->first) + "]-assay_refs\t";
      bool first(true);
      for (auto const & a : smd.assay_refs)
      {
        if (!first) { s += ","; } else { first = false; }
        s += "assay[" + String(a) + "]";
        }
      sl.push_back(s);
    }

    if (!smd.sample_refs.empty())
    {
      String s = "MTD\tstudy_variable[" + String(it->first) + String("]-sample_refs\t");
      bool first(true);
      for (auto const & a : smd.sample_refs)
      {
        if (!first) { s += ","; } else { first = false; }
        s += "sample[" + String(a) + "]";
        }
      sl.push_back(s);
    }

    if (!smd.description.isNull())
    {
      String s = "MTD\tstudy_variable[" + String(it->first) + String("]-description\t") + smd.description.toCellString();
      sl.push_back(s);
    }
  }

  for (map<Size, MzTabCVMetaData>::const_iterator it = md.cv.begin(); it != md.cv.end(); ++it)
  {
    const MzTabCVMetaData & mdcv = it->second;

    if (!mdcv.label.isNull())
    {
      String s = "MTD\tcv[" + String(it->first) + String("]-label\t") + mdcv.label.toCellString();
      sl.push_back(s);
    }

    if (!mdcv.full_name.isNull())
    {
      String s = "MTD\tcv[" + String(it->first) + String("]-full_name\t") + mdcv.full_name.toCellString();
      sl.push_back(s);
    }

    if (!mdcv.version.isNull())
    {
      String s = "MTD\tcv[" + String(it->first) + String("]-version\t") + mdcv.version.toCellString();
      sl.push_back(s);
    }

    if (!mdcv.url.isNull())
    {
      String s = "MTD\tcv[" + String(it->first) + String("]-url\t") + mdcv.url.toCellString();
      sl.push_back(s);
    }
  }

  // colunit-protein
  for (Size i = 0; i != md.colunit_protein.size(); ++i)
  {
    String s = String("MTD\tcolunit-protein") + md.colunit_protein[i];
    sl.push_back(s);
  }

  // colunit-peptide
  for (Size i = 0; i != md.colunit_peptide.size(); ++i)
  {
    String s = String("MTD\tcolunit-peptide") + md.colunit_peptide[i];
    sl.push_back(s);
  }

  // colunit-PSM
  for (Size i = 0; i != md.colunit_psm.size(); ++i)
  {
    String s = String("MTD\tcolunit-PSM") + md.colunit_psm[i];
    sl.push_back(s);
  }

  // colunit-small_molecule
  for (Size i = 0; i != md.colunit_small_molecule.size(); ++i)
  {
    String s = String("MTD\tcolunit-small_molecule") + md.colunit_small_molecule[i];
    sl.push_back(s);
  }
  }

  String MzTabFile::generateMzTabProteinHeader_(
      const MzTabProteinSectionRow& reference_row,
      const Size n_best_search_engine_scores,
      const std::vector<String>& optional_columns,
      const MzTabMetaData& meta, 
      size_t& n_columns) const
  {
    Size n_search_engine_scores = reference_row.search_engine_score_ms_run.size();

    StringList header;
    header.push_back("PRH");
    header.push_back("accession");
    header.push_back("description");
    header.push_back("taxid");
    header.push_back("species");
    header.push_back("database");
    header.push_back("database_version");
    header.push_back("search_engine");

    for (Size i = 0; i != n_best_search_engine_scores; ++i)
    {
      header.push_back(String("best_search_engine_score[") + String(i + 1) + String("]"));
    }

    if (n_search_engine_scores != 0)
    {
      // get number of runs for the first search score type (should be the same for every score)
      for (Size i = 0; i != reference_row.search_engine_score_ms_run.begin()->second.size(); ++i)
      {
        for (std::map<Size, std::map<Size, MzTabDouble> >::const_iterator search_it = reference_row.search_engine_score_ms_run.begin(); search_it != reference_row.search_engine_score_ms_run.end(); ++search_it)
        {
          header.push_back(String("search_engine_score[" + String(search_it->first) + "]_ms_run[") + String(i + 1) + String("]"));
        }
      }
    }

    if (store_protein_reliability_)
    {
      header.push_back("reliability");
    }

    for (std::map<Size, MzTabInteger>::const_iterator it = reference_row.num_psms_ms_run.begin(); it != reference_row.num_psms_ms_run.end(); ++it)
    {
      header.push_back(String("num_psms_ms_run[") + String(it->first) + String("]"));
    }

    for (std::map<Size, MzTabInteger>::const_iterator it = reference_row.num_peptides_distinct_ms_run.begin(); it != reference_row.num_peptides_distinct_ms_run.end(); ++it)
    {
      header.push_back(String("num_peptides_distinct_ms_run[") + String(it->first) + String("]"));
    }

    for (std::map<Size, MzTabInteger>::const_iterator it = reference_row.num_peptides_unique_ms_run.begin(); it != reference_row.num_peptides_unique_ms_run.end(); ++it)
    {
      header.push_back(String("num_peptides_unique_ms_run[") + String(it->first) + String("]"));
    }

    header.push_back("ambiguity_members");
    header.push_back("modifications");

    if (store_protein_uri_)
    {
      header.push_back("uri");
    }

    if (store_protein_goterms_)
    {
      header.push_back("go_terms");
    }

    header.push_back("protein_coverage");

    for (const auto& a : meta.assay)
    {
      header.push_back(String("protein_abundance_assay[") + String(a.first) + String("]"));
    }

    for (const auto& s : meta.study_variable)
    {
      header.push_back(String("protein_abundance_study_variable[") + String(s.first) + String("]"));
      header.push_back(String("protein_abundance_stdev_study_variable[") + String(s.first) + String("]"));
      header.push_back(String("protein_abundance_std_error_study_variable[") + String(s.first) + String("]"));
    }

    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  String MzTabFile::generateMzTabSectionRow_(
    const MzTabProteinSectionRow& row, 
    const vector<String>& optional_columns,
    const MzTabMetaData& meta, 
    size_t& n_columns) const
  {
  StringList s;
  s.push_back("PRT");
  s.push_back(row.accession.toCellString());
  s.push_back(row.description.toCellString());
  s.push_back(row.taxid.toCellString());
  s.push_back(row.species.toCellString());
  s.push_back(row.database.toCellString());
  s.push_back(row.database_version.toCellString());
  s.push_back(row.search_engine.toCellString());

  for (map<Size, MzTabDouble>::const_iterator it = row.best_search_engine_score.begin(); it != row.best_search_engine_score.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  for (std::map<Size, std::map<Size, MzTabDouble> >::const_iterator it = row.search_engine_score_ms_run.begin(); it != row.search_engine_score_ms_run.end(); ++it)
  {
    for (std::map<Size, MzTabDouble>::const_iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
    {
      s.push_back(sit->second.toCellString());
    }
  }

  if (store_protein_reliability_)
  {
    s.push_back(row.reliability.toCellString());
  }

  for (std::map<Size, MzTabInteger>::const_iterator it = row.num_psms_ms_run.begin(); it != row.num_psms_ms_run.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  for (std::map<Size, MzTabInteger>::const_iterator it = row.num_peptides_distinct_ms_run.begin(); it != row.num_peptides_distinct_ms_run.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  for (std::map<Size, MzTabInteger>::const_iterator it = row.num_peptides_unique_ms_run.begin(); it != row.num_peptides_unique_ms_run.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  s.push_back(row.ambiguity_members.toCellString());
  s.push_back(row.modifications.toCellString());

  if (store_protein_uri_)
  {
    s.push_back(row.uri.toCellString());
  }

  if (store_protein_goterms_)
  {
    s.push_back(row.go_terms.toCellString());
  }

  s.push_back(row.coverage.toCellString());

  // Quantification columns:
  // Assays
  for (const auto& kv : meta.assay)
  {
    const auto& k = kv.first;
    const auto& ity = row.protein_abundance_assay.find(k);

    // make sure that everything is set for this SV
    if (ity != row.protein_abundance_assay.end())
    {
      s.emplace_back(ity->second.toCellString());
    }
    else
    {
      // putting "null" directly would be a little faster, but with this we can keep consistency
      s.emplace_back(MzTabString().toCellString());
    }
  }

  // Study variables
  // go over all study variables that should be present and fill with either values
  // or uninitialized MzTabDouble()
  for (const auto& kv : meta.study_variable)
  {
    const auto& k = kv.first;
    const auto& ity = row.protein_abundance_study_variable.find(k);
    const auto& sd = row.protein_abundance_stdev_study_variable.find(k);
    const auto& err = row.protein_abundance_std_error_study_variable.find(k);

    // make sure that everything is set for this SV
    if (ity != row.protein_abundance_study_variable.end() &&
        sd != row.protein_abundance_stdev_study_variable.end() &&
        err != row.protein_abundance_std_error_study_variable.end())
    {
      s.emplace_back(ity->second.toCellString());
      s.emplace_back(sd->second.toCellString());
      s.emplace_back(err->second.toCellString());
    }
    else
    {
      // putting "null" directly would be a little faster, but with this we can keep consistency
      s.emplace_back(MzTabString().toCellString());
      s.emplace_back(MzTabString().toCellString());
      s.emplace_back(MzTabString().toCellString());
    }
  }

    addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();
    return ListUtils::concatenate(s, "\t");
  }

  String MzTabFile::generateMzTabPeptideHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_scores, Size assays, Size study_variables, const vector<String>& optional_columns,
    size_t& n_columns) const
  {
    StringList header;
    header.push_back("PEH");
    header.push_back("sequence");
    header.push_back("accession");
    header.push_back("unique");
    header.push_back("database");
    header.push_back("database_version");
    header.push_back("search_engine");

    for (Size i = 0; i != n_best_search_engine_scores; ++i)
    {
      header.push_back(String("best_search_engine_score[") + String(i + 1) + String("]"));
    }

    for (Size i = 0; i != search_ms_runs; ++i)
    {
      for (Size j = 0; j != n_search_engine_scores; ++j)
      {
        header.push_back(String("search_engine_score[" + String(j + 1) + "]_ms_run[") + String(i + 1) + String("]"));
      }
    }

    if (store_peptide_reliability_)
    {
      header.push_back("reliability");
    }

    header.push_back("modifications");
    header.push_back("retention_time");
    header.push_back("retention_time_window");
    header.push_back("charge");
    header.push_back("mass_to_charge");

    if (store_peptide_uri_)
    {
      header.push_back("uri");
    }

    header.push_back("spectra_ref");

    for (Size i = 0; i != assays; ++i)
    {
      header.push_back(String("peptide_abundance_assay[") + String(i + 1) + String("]"));
    }

    for (Size i = 0; i != study_variables; ++i)
    {
      header.push_back(String("peptide_abundance_study_variable[") + String(i + 1) + String("]"));
      header.push_back(String("peptide_abundance_stdev_study_variable[") + String(i + 1) + String("]"));
      header.push_back(String("peptide_abundance_std_error_study_variable[") + String(i + 1) + String("]"));
    }

    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  String MzTabFile::generateMzTabPSMHeader_(Size n_search_engine_scores, const vector<String>& optional_columns, size_t& n_columns) const
  {
  StringList header;
  header.push_back("PSH");
  header.push_back("sequence");
  header.push_back("PSM_ID");
  header.push_back("accession");
  header.push_back("unique");
  header.push_back("database");
  header.push_back("database_version");
  header.push_back("search_engine");

  for (Size i = 0; i != n_search_engine_scores; ++i)
  {
    header.push_back("search_engine_score[" + String(i + 1) + "]");
  }

  if (store_psm_reliability_)
  {
    header.push_back("reliability");
  }

  header.push_back("modifications");
  header.push_back("retention_time");
  header.push_back("charge");
  header.push_back("exp_mass_to_charge");
  header.push_back("calc_mass_to_charge");

  if (store_psm_uri_)
  {
    header.push_back("uri");
  }

  header.push_back("spectra_ref");
  header.push_back("pre");
  header.push_back("post");
  header.push_back("start");
  header.push_back("end");

  std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
  n_columns = header.size();
  return ListUtils::concatenate(header, "\t");
  }

  String MzTabFile::generateMzTabSectionRow_(const MzTabPeptideSectionRow& row, const vector<String>& optional_columns,
                                             const MzTabMetaData& /*meta*/, size_t& n_columns) const
  {
  StringList s;
  s.push_back("PEP");
  s.push_back(row.sequence.toCellString());
  s.push_back(row.accession.toCellString());
  s.push_back(row.unique.toCellString());
  s.push_back(row.database.toCellString());
  s.push_back(row.database_version.toCellString());
  s.push_back(row.search_engine.toCellString());

  // best score(s)
  for (auto const & bs : row.best_search_engine_score) { s.push_back(bs.second.toCellString()); }

  // run-level best score(s)
  for (auto const & bsrun : row.search_engine_score_ms_run)
  {
    for (auto const & bs : bsrun.second) { s.push_back(bs.second.toCellString()); }
  }

  if (store_peptide_reliability_)
  {
    s.push_back(row.reliability.toCellString());
  }

  s.push_back(row.modifications.toCellString());
  s.push_back(row.retention_time.toCellString());
  s.push_back(row.retention_time_window.toCellString());
  s.push_back(row.charge.toCellString());
  s.push_back(row.mass_to_charge.toCellString());

  if (store_peptide_uri_)
  {
    s.push_back(row.uri.toCellString());
  }

  s.push_back(row.spectra_ref.toCellString());

  // quantification columns
  for (std::map<Size, MzTabDouble>::const_iterator it = row.peptide_abundance_assay.begin(); it != row.peptide_abundance_assay.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  std::map<Size, MzTabDouble>::const_iterator sv_it = row.peptide_abundance_study_variable.begin();
  std::map<Size, MzTabDouble>::const_iterator sv_stdev_it = row.peptide_abundance_stdev_study_variable.begin();
  std::map<Size, MzTabDouble>::const_iterator sv_error_it = row.peptide_abundance_std_error_study_variable.begin();

  for (;
       sv_it != row.peptide_abundance_study_variable.end()
       && sv_stdev_it != row.peptide_abundance_stdev_study_variable.end()
       && sv_error_it != row.peptide_abundance_std_error_study_variable.end();
       ++sv_it, ++sv_stdev_it, ++sv_error_it)
  {
    s.push_back(sv_it->second.toCellString());
    s.push_back(sv_stdev_it->second.toCellString());
    s.push_back(sv_error_it->second.toCellString());
  }

    addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();

    return ListUtils::concatenate(s, "\t");
  }

  String MzTabFile::generateMzTabSectionRow_(const MzTabPSMSectionRow& row, const vector<String>& optional_columns,
                                             const MzTabMetaData& /*meta*/, size_t& n_columns) const
  {
    StringList s;
    s.push_back("PSM");
    s.push_back(row.sequence.toCellString());
    s.push_back(row.PSM_ID.toCellString());
    s.push_back(row.accession.toCellString());
    s.push_back(row.unique.toCellString());
    s.push_back(row.database.toCellString());
    s.push_back(row.database_version.toCellString());
    s.push_back(row.search_engine.toCellString());

    if (row.search_engine_score.empty())
    { // workaround for PeptideIDs without hits (QC export)
      s.push_back("null");
    }
    else
    {
      for (map<Size, MzTabDouble>::const_iterator it = row.search_engine_score.begin(); it != row.search_engine_score.end(); ++it)
      {
        s.push_back(it->second.toCellString());
      }
    }

    if (store_psm_reliability_)
    {
      s.push_back(row.reliability.toCellString());
    }

    s.push_back(row.modifications.toCellString());
    s.push_back(row.retention_time.toCellString());
    s.push_back(row.charge.toCellString());
    s.push_back(row.exp_mass_to_charge.toCellString());
    s.push_back(row.calc_mass_to_charge.toCellString());

    if (store_psm_uri_)
    {
      s.push_back(row.uri.toCellString());
    }

    s.push_back(row.spectra_ref.toCellString());
    s.push_back(row.pre.toCellString());
    s.push_back(row.post.toCellString());
    s.push_back(row.start.toCellString());
    s.push_back(row.end.toCellString());

    addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();
    return ListUtils::concatenate(s, "\t");
  }

  String MzTabFile::generateMzTabSmallMoleculeHeader_(Size ms_runs, Size n_best_search_engine_scores, Size n_search_engine_scores, Size assays, Size study_variables, const vector<String>& optional_smallmolecule_columns, size_t& n_columns) const
  {
    StringList header;
    header.push_back("SMH");
    header.push_back("identifier");
    header.push_back("chemical_formula");
    header.push_back("smiles");
    header.push_back("inchi_key");
    header.push_back("description");
    header.push_back("exp_mass_to_charge");
    header.push_back("calc_mass_to_charge");
    header.push_back("charge");
    header.push_back("retention_time");
    header.push_back("taxid");
    header.push_back("species");
    header.push_back("database");
    header.push_back("database_version");

    if (store_smallmolecule_reliability_)
    {
      header.push_back("reliability");
    }

    if (store_smallmolecule_uri_)
    {
      header.push_back("uri");
    }

    header.push_back("spectra_ref");
    header.push_back("search_engine");

    for (Size i = 0; i != n_best_search_engine_scores; ++i)
    {
      header.push_back(String("best_search_engine_score[") + String(i + 1) + String("]"));
    }

    for (Size i = 0; i != ms_runs; ++i)
    {
      for (Size j = 0; j != n_search_engine_scores; ++j)
      {
        header.push_back(String("search_engine_score[" + String(j + 1) + "]_ms_run[") + String(i + 1) + String("]"));
      }
    }

    header.push_back("modifications");

    for (Size i = 0; i != assays; ++i)
    {
      header.push_back(String("smallmolecule_abundance_assay[") + String(i + 1) + String("]"));
    }

    for (Size i = 0; i != study_variables; ++i)
    {
      header.push_back(String("smallmolecule_abundance_study_variable[") + String(i + 1) + String("]"));
      header.push_back(String("smallmolecule_abundance_stdev_study_variable[") + String(i + 1) + String("]"));
      header.push_back(String("smallmolecule_abundance_std_error_study_variable[") + String(i + 1) + String("]"));
    }

    // copy optional column names to header
    std::copy(optional_smallmolecule_columns.begin(), optional_smallmolecule_columns.end(), std::back_inserter(header));

    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  String MzTabFile::generateMzTabSectionRow_(
      const MzTabSmallMoleculeSectionRow& row,
      const std::vector<String>& optional_columns,
      const MzTabMetaData& /*meta*/, size_t& n_columns) const
  {
  StringList s;
  s.push_back("SML");
  s.push_back(row.identifier.toCellString());
  s.push_back(row.chemical_formula.toCellString());
  s.push_back(row.smiles.toCellString());
  s.push_back(row.inchi_key.toCellString());
  s.push_back(row.description.toCellString());
  s.push_back(row.exp_mass_to_charge.toCellString());
  s.push_back(row.calc_mass_to_charge.toCellString());
  s.push_back(row.charge.toCellString());
  s.push_back(row.retention_time.toCellString());
  s.push_back(row.taxid.toCellString());
  s.push_back(row.species.toCellString());
  s.push_back(row.database.toCellString());
  s.push_back(row.database_version.toCellString());

  if (store_smallmolecule_reliability_)
  {
    s.push_back(row.reliability.toCellString());
  }

  if (store_smallmolecule_uri_)
  {
    s.push_back(row.uri.toCellString());
  }

  s.push_back(row.spectra_ref.toCellString());
  s.push_back(row.search_engine.toCellString());

  for (map<Size, MzTabDouble>::const_iterator it = row.best_search_engine_score.begin(); it != row.best_search_engine_score.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  for (auto it = row.search_engine_score_ms_run.begin(); it != row.search_engine_score_ms_run.end(); ++it)
  {
    for (auto sit = it->second.begin(); sit != it->second.end(); ++sit)
    {
      s.push_back(sit->second.toCellString());
    }
  }

  s.push_back(row.modifications.toCellString());

  // quantification columns
  std::map<Size, MzTabDouble>::const_iterator sv_it = row.smallmolecule_abundance_study_variable.begin();
  std::map<Size, MzTabDouble>::const_iterator sv_stdev_it = row.smallmolecule_abundance_stdev_study_variable.begin();
  std::map<Size, MzTabDouble>::const_iterator sv_error_it = row.smallmolecule_abundance_std_error_study_variable.begin();

  for (;
       sv_it != row.smallmolecule_abundance_study_variable.end()
       && sv_stdev_it != row.smallmolecule_abundance_stdev_study_variable.end()
       && sv_error_it != row.smallmolecule_abundance_std_error_study_variable.end();
       ++sv_it, ++sv_stdev_it, ++sv_error_it)
  {
    s.push_back(sv_it->second.toCellString());
    s.push_back(sv_stdev_it->second.toCellString());
    s.push_back(sv_error_it->second.toCellString());
  }

    addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();
    return ListUtils::concatenate(s, "\t");
  }

  String MzTabFile::generateMzTabNucleicAcidHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_scores, const std::vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList header;
    header.push_back("NUH");
    header.push_back("accession");
    header.push_back("description");
    header.push_back("taxid");
    header.push_back("species");
    header.push_back("database");
    header.push_back("database_version");
    header.push_back("search_engine");

    for (Size i = 0; i != n_best_search_engine_scores; ++i)
    {
      header.push_back(String("best_search_engine_score[") + String(i + 1) + String("]"));
    }

    for (Size i = 0; i != search_ms_runs; ++i)
    {
      for (Size j = 0; j != n_search_engine_scores; ++j)
      {
        header.push_back(String("search_engine_score[" + String(j + 1) + "]_ms_run[") + String(i + 1) + String("]"));
      }
    }

    if (store_nucleic_acid_reliability_)
    {
      header.push_back("reliability");
    }

    for (Size i = 0; i != search_ms_runs; ++i)
    {
      header.push_back(String("num_osms_ms_run[") + String(i) + String("]"));
    }

    for (Size i = 0; i != search_ms_runs; ++i)
    {
      header.push_back(String("num_oligos_distinct_ms_run[") + String(i) + String("]"));
    }

    for (Size i = 0; i != search_ms_runs; ++i)
    {
      header.push_back(String("num_oligos_unique_ms_run[") + String(i) + String("]"));
    }

    header.push_back("ambiguity_members");
    header.push_back("modifications");

    if (store_nucleic_acid_uri_)
    {
      header.push_back("uri");
    }

    if (store_nucleic_acid_goterms_)
    {
      header.push_back("go_terms");
    }

    header.push_back("sequence_coverage");

    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  String MzTabFile::generateMzTabSectionRow_(
      const MzTabNucleicAcidSectionRow& row,
      const vector<String>& optional_columns,
      const MzTabMetaData& /*meta*/, size_t& n_columns) const
  {
    StringList s;
    s.push_back("NUC");
    s.push_back(row.accession.toCellString());
    s.push_back(row.description.toCellString());
    s.push_back(row.taxid.toCellString());
    s.push_back(row.species.toCellString());
    s.push_back(row.database.toCellString());
    s.push_back(row.database_version.toCellString());
    s.push_back(row.search_engine.toCellString());

    for (map<Size, MzTabDouble>::const_iterator it = row.best_search_engine_score.begin(); it != row.best_search_engine_score.end(); ++it)
    {
      s.push_back(it->second.toCellString());
    }

    for (std::map<Size, std::map<Size, MzTabDouble> >::const_iterator it = row.search_engine_score_ms_run.begin(); it != row.search_engine_score_ms_run.end(); ++it)
    {
      for (std::map<Size, MzTabDouble>::const_iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
      {
        s.push_back(sit->second.toCellString());
      }
    }

    if (store_nucleic_acid_reliability_)
    {
      s.push_back(row.reliability.toCellString());
    }

    for (std::map<Size, MzTabInteger>::const_iterator it = row.num_osms_ms_run.begin(); it != row.num_osms_ms_run.end(); ++it)
    {
      s.push_back(it->second.toCellString());
    }

    for (std::map<Size, MzTabInteger>::const_iterator it = row.num_oligos_distinct_ms_run.begin(); it != row.num_oligos_distinct_ms_run.end(); ++it)
    {
      s.push_back(it->second.toCellString());
    }

    for (std::map<Size, MzTabInteger>::const_iterator it = row.num_oligos_unique_ms_run.begin(); it != row.num_oligos_unique_ms_run.end(); ++it)
    {
      s.push_back(it->second.toCellString());
    }

    s.push_back(row.ambiguity_members.toCellString());
    s.push_back(row.modifications.toCellString());

    if (store_nucleic_acid_uri_)
    {
      s.push_back(row.uri.toCellString());
    }

    if (store_nucleic_acid_goterms_)
    {
      s.push_back(row.go_terms.toCellString());
    }

    s.push_back(row.coverage.toCellString());

    addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();
    return ListUtils::concatenate(s, "\t");
  }

  String MzTabFile::generateMzTabOligonucleotideHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_scores, const vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList header;
    header.push_back("OLH");
    header.push_back("sequence");
    header.push_back("accession");
    header.push_back("unique");
    header.push_back("search_engine");

    for (Size i = 0; i != n_best_search_engine_scores; ++i)
    {
      header.push_back(String("best_search_engine_score[") + String(i + 1) + String("]"));
    }

    for (Size i = 0; i != search_ms_runs; ++i)
    {
      for (Size j = 0; j != n_search_engine_scores; ++j)
      {
        header.push_back(String("search_engine_score[" + String(j + 1) + "]_ms_run[") + String(i + 1) + String("]"));
      }
    }

    if (store_oligonucleotide_reliability_)
    {
      header.push_back("reliability");
    }

    header.push_back("modifications");
    header.push_back("retention_time");
    header.push_back("retention_time_window");

    if (store_oligonucleotide_uri_)
    {
      header.push_back("uri");
    }

    header.push_back("pre");
    header.push_back("post");
    header.push_back("start");
    header.push_back("end");

    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  String MzTabFile::generateMzTabSectionRow_(
      const MzTabOligonucleotideSectionRow& row,
      const vector<String>& optional_columns,
      const MzTabMetaData& /*meta*/, size_t& n_columns) const
  {
    StringList s;
    s.push_back("OLI");
    s.push_back(row.sequence.toCellString());
    s.push_back(row.accession.toCellString());
    s.push_back(row.unique.toCellString());
    s.push_back(row.search_engine.toCellString());

    for (map<Size, MzTabDouble>::const_iterator it = row.best_search_engine_score.begin(); it != row.best_search_engine_score.end(); ++it)
    {
      s.push_back(it->second.toCellString());
    }

    for (map<Size, map<Size, MzTabDouble> >::const_iterator it = row.search_engine_score_ms_run.begin(); it != row.search_engine_score_ms_run.end(); ++it)
    {
      for (map<Size, MzTabDouble>::const_iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
      {
        s.push_back(sit->second.toCellString());
      }
    }

    if (store_oligonucleotide_reliability_)
    {
      s.push_back(row.reliability.toCellString());
    }

    s.push_back(row.modifications.toCellString());
    s.push_back(row.retention_time.toCellString());
    s.push_back(row.retention_time_window.toCellString());

    if (store_oligonucleotide_uri_)
    {
      s.push_back(row.uri.toCellString());
    }

    s.push_back(row.pre.toCellString());
    s.push_back(row.post.toCellString());
    s.push_back(row.start.toCellString());
    s.push_back(row.end.toCellString());

    addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();

    return ListUtils::concatenate(s, "\t");
  }

  String MzTabFile::generateMzTabOSMHeader_(Size n_search_engine_scores, const vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList header;
    header.push_back("OSH");
    header.push_back("sequence");
    header.push_back("search_engine");

    for (Size i = 0; i != n_search_engine_scores; ++i)
    {
      header.push_back("search_engine_score[" + String(i + 1) + "]");
    }

    if (store_osm_reliability_)
    {
      header.push_back("reliability");
    }

    header.push_back("modifications");
    header.push_back("retention_time");
    header.push_back("charge");
    header.push_back("exp_mass_to_charge");
    header.push_back("calc_mass_to_charge");

    if (store_osm_uri_)
    {
      header.push_back("uri");
    }

    header.push_back("spectra_ref");

    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  String MzTabFile::generateMzTabSectionRow_(
      const MzTabOSMSectionRow& row,
      const vector<String>& optional_columns,
      const MzTabMetaData& /*meta*/, size_t& n_columns) const
  {
    StringList s;
    s.push_back("OSM");
    s.push_back(row.sequence.toCellString());
    s.push_back(row.search_engine.toCellString());

    for (map<Size, MzTabDouble>::const_iterator it = row.search_engine_score.begin(); it != row.search_engine_score.end(); ++it)
    {
      s.push_back(it->second.toCellString());
    }

    if (store_osm_reliability_)
    {
      s.push_back(row.reliability.toCellString());
    }

    s.push_back(row.modifications.toCellString());
    s.push_back(row.retention_time.toCellString());
    s.push_back(row.charge.toCellString());
    s.push_back(row.exp_mass_to_charge.toCellString());
    s.push_back(row.calc_mass_to_charge.toCellString());

    if (store_osm_uri_)
    {
      s.push_back(row.uri.toCellString());
    }

    s.push_back(row.spectra_ref.toCellString());

    addOptionalColumnsToSectionRow_(optional_columns, row.opt_, s);
    n_columns = s.size();
    return ListUtils::concatenate(s, "\t");
  }

  void MzTabFile::addOptionalColumnsToSectionRow_(const vector<String>& column_names, const vector<MzTabOptionalColumnEntry>& column_entries, StringList& output)
  {
    for (vector<String>::const_iterator it = column_names.begin(); it != column_names.end(); ++it)
    {
      bool found = false;
      for (Size i = 0; i != column_entries.size(); ++i)
      {
        if (column_entries[i].first == *it)
        {
          output.push_back(column_entries[i].second.toCellString());
          found = true;
          break;
        }
      }
      if (!found)
      {
        output.push_back(MzTabString("null").toCellString());
      }
    }
  }
    // stream IDs to file
  void MzTabFile::store(
        const String& filename,
        const std::vector<ProteinIdentification>& protein_identifications,
        const std::vector<PeptideIdentification>& peptide_identifications,
        bool first_run_inference_only,
        bool export_empty_pep_ids,
        bool export_all_psms,
        const String& title)
  {
    if (!(FileHandler::hasValidExtension(filename, FileTypes::MZTAB) || FileHandler::hasValidExtension(filename, FileTypes::TSV)))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '"
      + FileTypes::typeToName(FileTypes::MZTAB) + "' or '" + FileTypes::typeToName(FileTypes::TSV) + "'");
    }

    vector<const PeptideIdentification*> pep_ids_ptr;
    pep_ids_ptr.reserve(peptide_identifications.size());
    for (const PeptideIdentification& pi : peptide_identifications) { pep_ids_ptr.push_back(&pi); }

    vector<const ProteinIdentification*> prot_ids_ptr;
    prot_ids_ptr.reserve(protein_identifications.size());
    for (const ProteinIdentification& pi : protein_identifications) { prot_ids_ptr.push_back(&pi); }

    ofstream tab_file;
    tab_file.open(filename, ios::out | ios::trunc);

    MzTab::IDMzTabStream s(
      prot_ids_ptr,
      pep_ids_ptr,
      filename,
      first_run_inference_only,
      export_empty_pep_ids,
      export_all_psms,
      title);      

    // generate full meta data section and write to file
    MzTabMetaData meta_data = s.getMetaData();

    {
      StringList out;
      generateMzTabMetaDataSection_(meta_data, out);
      for (const String & line : out) { tab_file << line << "\n"; }
    }
   
    Size n_best_search_engine_score = meta_data.protein_search_engine_score.size();

    {
      MzTabProteinSectionRow row;
      bool first = true;
      size_t n_header_columns = 0;
      while (s.nextPRTRow(row))
      {
        if (first)
        { // add header
          tab_file << "\n" << generateMzTabProteinHeader_(
            row,
            n_best_search_engine_score,
            s.getProteinOptionalColumnNames(),
            meta_data,
            n_header_columns) + "\n";
          first = false;
        }
        size_t n_section_columns = 0;
        tab_file << generateMzTabSectionRow_(row, s.getProteinOptionalColumnNames(), meta_data, n_section_columns) + "\n";
        if (n_header_columns != n_section_columns)  throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Protein header and content differs in columns. Please report this bug to the OpenMS developers.");
      }
    }

    Size n_search_engine_scores = meta_data.psm_search_engine_score.size();

    if (n_search_engine_scores == 0)
    {
      OPENMS_LOG_WARN << "No search engine scores given. Please check your input data." << endl;
    }

    {
      MzTabPSMSectionRow row;
      bool first = true;
      size_t n_header_columns = 0;
      while (s.nextPSMRow(row))
      {
        // TODO better return a State enum instead of relying on some uninitialized
        // parts of a row.. at least it is a mandatory field and therefore it would not make
        // sense writing that row anyway
        if (!row.sequence.isNull())
        {
          if (first)
          { // add header
            tab_file << "\n" << generateMzTabPSMHeader_(n_search_engine_scores, s.getPSMOptionalColumnNames(), n_header_columns) + "\n";
            first = false;
          }
          size_t n_section_columns = 0;
          tab_file << generateMzTabSectionRow_(row, s.getPSMOptionalColumnNames(), meta_data, n_section_columns) + "\n";
          if (n_header_columns != n_section_columns)  throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "PSM header and content differs in columns. Please report this bug to the OpenMS developers.");
        }
      }
    }

    tab_file.close();
    
  }

  void MzTabFile::store(
      const String& filename, 
      const ConsensusMap& cmap,
      const bool first_run_inference_only,
      const bool export_unidentified_features,
      const bool export_unassigned_ids,
      const bool export_subfeatures,
      const bool export_empty_pep_ids,
      const bool export_all_psms) const
  {
    if (!(FileHandler::hasValidExtension(filename, FileTypes::MZTAB) || FileHandler::hasValidExtension(filename, FileTypes::TSV)))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '"
      + FileTypes::typeToName(FileTypes::MZTAB) + "' or '" + FileTypes::typeToName(FileTypes::TSV) + "'");
    }

    ofstream tab_file;
    tab_file.open(filename, ios::out | ios::trunc);

    MzTab::CMMzTabStream s(
      cmap,
      filename,
      first_run_inference_only,
      export_unidentified_features,
      export_unassigned_ids,
      export_subfeatures,
      export_empty_pep_ids,
      export_all_psms,
      "ConsensusMap export from OpenMS");      

    // generate full meta data section and write to file
    MzTabMetaData meta_data = s.getMetaData();

    {
      StringList out;
      generateMzTabMetaDataSection_(meta_data, out);
      for (const String & line : out) { tab_file << line << "\n"; }
    }
   
    Size n_best_search_engine_score = meta_data.protein_search_engine_score.size();
    // TODO: we currently only store one search engine score per PSM so we need to limit the number to the main score
    n_best_search_engine_score = std::min(n_best_search_engine_score, Size(1));
    {
      MzTabProteinSectionRow row;
      bool first = true;
      size_t n_header_columns = 0;
      while (s.nextPRTRow(row))
      {
        if (first)
        { // add header
          tab_file << "\n" << generateMzTabProteinHeader_(
            row,
            n_best_search_engine_score,
            s.getProteinOptionalColumnNames(),
            meta_data, 
            n_header_columns) + "\n";
          first = false;
        }
        size_t n_section_columns = 0;
        tab_file << generateMzTabSectionRow_(row, s.getProteinOptionalColumnNames(), meta_data, n_section_columns) + "\n";
        if (n_header_columns != n_section_columns)  throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Protein header and content differs in columns. Please report this bug to the OpenMS developers.");
      }
    }

    Size assays(0);
    Size study_variables(0);
    {
      MzTabPeptideSectionRow row;
      bool first = true;
      size_t n_header_columns = 0;
      while (s.nextPEPRow(row))
      {
        if (first)
        {
          assays = row.peptide_abundance_assay.size();
          study_variables = row.peptide_abundance_study_variable.size();
          Size n_search_engine_score = row.search_engine_score_ms_run.size(); // scores to runs          
          Size search_ms_runs = n_search_engine_score != 0 ? row.search_engine_score_ms_run.at(1).size() : 0; // take number of searched MS runs from first score. TODO: handle this more generic
          OPENMS_LOG_DEBUG << "Exporting assays: " << assays << endl;
          OPENMS_LOG_DEBUG << "Exporting study variables: " << study_variables << endl;
          OPENMS_LOG_DEBUG << "Exporting search engines scores: " << n_search_engine_score << endl;
          Size n_best_search_engine_score = row.best_search_engine_score.size();
          tab_file << "\n" << generateMzTabPeptideHeader_(search_ms_runs, n_best_search_engine_score, n_search_engine_score, assays, study_variables, s.getPeptideOptionalColumnNames(), n_header_columns) + "\n";
          first = false;
        }
        size_t n_section_columns = 0;
        tab_file << generateMzTabSectionRow_(row, s.getPeptideOptionalColumnNames(), meta_data, n_section_columns) + "\n";
        if (n_header_columns != n_section_columns)
        {
          OPENMS_LOG_ERROR << "Number of columns in header/section: " << n_header_columns << "/" << n_section_columns << endl;
          throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide header and content differs in columns. Please report this bug to the OpenMS developers.");
        }
      }
    } 

    Size n_search_engine_scores = meta_data.psm_search_engine_score.size();

    if (n_search_engine_scores == 0)
    {
      OPENMS_LOG_WARN << "No search engine scores given. Please check your input data." << endl;
    }

    {
      MzTabPSMSectionRow row;
      bool first = true;

      // TODO: we currently only store one search engine score per PSM so we need to limit the number to the main score      
      n_search_engine_scores = 1;
      size_t n_header_columns = 0;
      while (s.nextPSMRow(row))
      {
        // TODO better return a State enum instead of relying on some uninitialized
        // parts of a row.. at least it is a mandatory field and therefore it would not make
        // sense writing that row anyway
        if (!row.sequence.isNull())
        {
          if (first)
          { // add header
            tab_file << "\n" << generateMzTabPSMHeader_(n_search_engine_scores, s.getPSMOptionalColumnNames(), n_header_columns) + "\n";           
            first = false;
          }
          size_t n_section_columns = 0;
          tab_file << generateMzTabSectionRow_(row, s.getPSMOptionalColumnNames(), meta_data, n_section_columns) + "\n";
          if (n_header_columns != n_section_columns)  throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "PSM header and content differs in columns. Please report this bug to the OpenMS developers.");
        }
      }
    }

    tab_file.close();
  }

  void MzTabFile::store(const String& filename, const MzTab& mz_tab) const
  {
    if (!(FileHandler::hasValidExtension(filename, FileTypes::MZTAB) || FileHandler::hasValidExtension(filename, FileTypes::TSV)))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '"
      + FileTypes::typeToName(FileTypes::MZTAB) + "' or '" + FileTypes::typeToName(FileTypes::TSV) + "'");
    }

    StringList out;
    generateMzTabMetaDataSection_(mz_tab.getMetaData(), out);
    bool complete = (mz_tab.getMetaData().mz_tab_mode.toCellString() == "Complete");
    Size ms_runs = mz_tab.getMetaData().ms_run.size();

    const MzTabProteinSectionRows& protein_section = mz_tab.getProteinSectionRows();
    const MzTabPeptideSectionRows& peptide_section = mz_tab.getPeptideSectionRows();
    const MzTabPSMSectionRows& psm_section = mz_tab.getPSMSectionRows();
    const MzTabSmallMoleculeSectionRows& smallmolecule_section = mz_tab.getSmallMoleculeSectionRows();

    if (!protein_section.empty())
    {
      Size n_best_search_engine_score = mz_tab.getMetaData().protein_search_engine_score.size();

      // add header
      // use the first row as a reference row to get the optional cols from meta values
      out.push_back("");
      size_t n_header_columns = 0;
      out.push_back(generateMzTabProteinHeader_(protein_section[0],
          n_best_search_engine_score,
          mz_tab.getProteinOptionalColumnNames(),
          mz_tab.getMetaData(), 
          n_header_columns));

      // add section content
      generateMzTabSection_(protein_section, mz_tab.getProteinOptionalColumnNames(), mz_tab.getMetaData(), out, n_header_columns);
    }

    if (!peptide_section.empty())
    {
      Size assays = peptide_section[0].peptide_abundance_assay.size();
      Size study_variables = peptide_section[0].peptide_abundance_study_variable.size();
      Size search_ms_runs = 0;
      if (complete)
      {
        // all ms_runs mandatory
        search_ms_runs = ms_runs;
      }
      else // only report all scores if user provided at least one
      {
        const MzTabPeptideSectionRows& psr = mz_tab.getPeptideSectionRows();
        bool has_ms_run_level_scores = false;
        for (Size i = 0; i != psr.size(); ++i)
        {
          if (!psr[i].search_engine_score_ms_run.empty())
          {
            has_ms_run_level_scores = true;
          }
        }

        if (has_ms_run_level_scores) { search_ms_runs = ms_runs; }
      }
      Size n_search_engine_score = peptide_section[0].search_engine_score_ms_run.size();
      Size n_best_search_engine_score = peptide_section[0].best_search_engine_score.size();

      out.push_back("");
      size_t n_header_columns = 0;
      out.push_back(generateMzTabPeptideHeader_(search_ms_runs, n_best_search_engine_score, n_search_engine_score, assays, study_variables, mz_tab.getPeptideOptionalColumnNames(), n_header_columns));
      generateMzTabSection_(mz_tab.getPeptideSectionRows(), mz_tab.getPeptideOptionalColumnNames(), mz_tab.getMetaData(), out, n_header_columns);
    }

    if (!psm_section.empty())
    {
      //Size n_search_engine_scores = mz_tab.getMetaData().psm_search_engine_score.size();
      //TODO since we currently only support getting the main score for psms, maximally reserve one column
      Size n_search_engine_scores = std::min(mz_tab.getMetaData().psm_search_engine_score.size(), Size(1));
      if (n_search_engine_scores == 0)
      {
        // TODO warn
      }
      out.push_back("");
      size_t n_header_columns = 0;
      out.push_back(generateMzTabPSMHeader_(n_search_engine_scores, mz_tab.getPSMOptionalColumnNames(), n_header_columns));
      generateMzTabSection_(mz_tab.getPSMSectionRows(), mz_tab.getPSMOptionalColumnNames(), mz_tab.getMetaData(), out, n_header_columns);
    }

    if (!smallmolecule_section.empty())
    {
      Size assays = smallmolecule_section[0].smallmolecule_abundance_assay.size();
      Size study_variables = smallmolecule_section[0].smallmolecule_abundance_study_variable.size();
      Size n_search_engine_score = smallmolecule_section[0].search_engine_score_ms_run.size();
      Size n_best_search_engine_score = mz_tab.getMetaData().smallmolecule_search_engine_score.size();
      out.push_back("");
      size_t n_header_columns = 0;
      out.push_back(generateMzTabSmallMoleculeHeader_(ms_runs, n_best_search_engine_score, n_search_engine_score, assays, study_variables, mz_tab.getSmallMoleculeOptionalColumnNames(), n_header_columns));
      generateMzTabSection_(smallmolecule_section, mz_tab.getSmallMoleculeOptionalColumnNames(), mz_tab.getMetaData(), out, n_header_columns);
    }

    const MzTabNucleicAcidSectionRows& nucleic_acid_section = mz_tab.getNucleicAcidSectionRows();
    const MzTabOligonucleotideSectionRows& oligonucleotide_section = mz_tab.getOligonucleotideSectionRows();
    const MzTabOSMSectionRows& osm_section = mz_tab.getOSMSectionRows();

    if (!nucleic_acid_section.empty())
    {
      Size search_ms_runs = 0;
      if (complete)
      {
        // all ms_runs mandatory
        search_ms_runs = ms_runs;
      }
      else // only report all scores if user provided at least one
      {
        bool has_ms_run_level_scores = false;
        for (Size i = 0; i != nucleic_acid_section.size(); ++i)
        {
          if (!nucleic_acid_section[i].search_engine_score_ms_run.empty())
          {
            has_ms_run_level_scores = true;
          }
        }

        if (has_ms_run_level_scores)
        {
          search_ms_runs = ms_runs;
        }
      }
      Size n_search_engine_score = nucleic_acid_section[0].search_engine_score_ms_run.size();
      Size n_best_search_engine_score = mz_tab.getMetaData().nucleic_acid_search_engine_score.size();

      // add header
      out.push_back("");
      size_t n_header_columns = 0;
      out.push_back(generateMzTabNucleicAcidHeader_(search_ms_runs, n_search_engine_score, n_best_search_engine_score, mz_tab.getNucleicAcidOptionalColumnNames(), n_header_columns));

      // add section
      generateMzTabSection_(nucleic_acid_section, mz_tab.getNucleicAcidOptionalColumnNames(), mz_tab.getMetaData(), out, n_header_columns);
    }

    if (!oligonucleotide_section.empty())
    {
      Size search_ms_runs = 0;
      if (complete)
      {
        // all ms_runs mandatory
        search_ms_runs = ms_runs;
      }
      else // only report all scores if user provided at least one
      {
        bool has_ms_run_level_scores = false;
        for (Size i = 0; i != oligonucleotide_section.size(); ++i)
        {
          if (!oligonucleotide_section[i].search_engine_score_ms_run.empty())
          {
            has_ms_run_level_scores = true;
          }
        }

        if (has_ms_run_level_scores)
        {
          search_ms_runs = ms_runs;
        }
      }
      Size n_search_engine_score = oligonucleotide_section[0].search_engine_score_ms_run.size();
      Size n_best_search_engine_score = mz_tab.getMetaData().oligonucleotide_search_engine_score.size();
      out.push_back("");
      size_t n_columns = 0;
      out.push_back(generateMzTabOligonucleotideHeader_(search_ms_runs, n_best_search_engine_score, n_search_engine_score, mz_tab.getOligonucleotideOptionalColumnNames(), n_columns));
      generateMzTabSection_(mz_tab.getOligonucleotideSectionRows(), mz_tab.getOligonucleotideOptionalColumnNames(), mz_tab.getMetaData(), out, n_columns);
    }

    if (!osm_section.empty())
    {
      Size n_search_engine_scores = mz_tab.getMetaData().osm_search_engine_score.size();

      if (n_search_engine_scores == 0)
      {
        // TODO warn
      }
      out.push_back("");
      size_t n_columns = 0;
      out.push_back(generateMzTabOSMHeader_(n_search_engine_scores, mz_tab.getOSMOptionalColumnNames(), n_columns));
      generateMzTabSection_(mz_tab.getOSMSectionRows(), mz_tab.getOSMOptionalColumnNames(), mz_tab.getMetaData(), out, n_columns);
  }

    // insert comments (might provide critical cues for human reader) and empty lines
  Size line = 0;
  vector<Size> empty_rows = mz_tab.getEmptyRows();
  map<Size, String> comment_rows = mz_tab.getCommentRows();

  if (empty_rows.empty() && comment_rows.empty())
  {
    TextFile tmp_out;
    for (TextFile::ConstIterator it = out.begin(); it != out.end(); ++it)
    {
      tmp_out.addLine(*it);
    }
    tmp_out.store(filename);
  }
  else
  {
    TextFile tmp_out;
    for (TextFile::ConstIterator it = out.begin(); it != out.end(); )
    {
      if (std::binary_search(empty_rows.begin(), empty_rows.end(), line))  // check if current line was originally an empty line
      {
        tmp_out.addLine("\n");
        ++line;
      }
      else if (comment_rows.find(line) != comment_rows.end()) // check if current line was originally a comment line
      {
        tmp_out.addLine(comment_rows[line]);
        ++line;
      }
      else   // no empty line, no comment => add row
      {
        tmp_out.addLine(*it);
        ++line;
        ++it;
      }
    }
    tmp_out.store(filename);
  }
  }

}

#pragma clang diagnostic pop
