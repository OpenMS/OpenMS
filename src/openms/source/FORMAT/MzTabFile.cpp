// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <QtCore/QRegExp>
#include <fstream>

namespace OpenMS
{

MzTabFile::MzTabFile()
{
}

MzTabFile::~MzTabFile()
{
}

std::pair<int, int> MzTabFile::extractIndexPairsFromBrackets_(const String & s)
{
  QString qs = s.toQString();
  std::pair<Int, Int> pair(0,0);

  // ^        # Match the start of the line
  // .*?      # Non-greedy match anything
  // \[       # Upto the first opening bracket (escaped)
  // (\d+)    # Match a digit string (one or more)
  // \]       # Match closing bracket
  // .*       # Match the rest of the line
  // $        # Match the end of the line
  QRegExp rx_first_number(QString("^.*?\\[(\\d+)\\].*$"));

  // as above but capture on the second bracket
  QRegExp rx_second_number(QString("^.*?\\[\\d+\\].*.*?\\[(\\d+)\\].*$"));

  Int pos_first_number = rx_first_number.indexIn(qs);
  if (pos_first_number > -1)
  {
    QString value = rx_first_number.cap(1);
    pair.first = value.toInt();
  }

  Int pos_second_number = rx_second_number.indexIn(qs);
  if (pos_second_number > -1)
  {
    QString value = rx_second_number.cap(1);
    pair.second = value.toInt();
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

  map<String, Size> protein_custom_opt_columns;  // map column name to original column index
  map<String, Size> peptide_custom_opt_columns;
  map<String, Size> psm_custom_opt_columns;
  map<String, Size> smallmolecule_custom_opt_columns;

  // mandatory meta values / columns check
  bool has_mzTab_version = false;
  bool has_mzTab_mode = false;
  bool has_mzTab_type = false;
  bool has_description = false;
  Size count_search_engine_score = 0;
  Size count_fixed_mod = 0;
  Size count_variable_mod = 0;
  bool has_software = false;
  bool has_protein_quantification_unit = false;
  bool has_peptide_quantification_unit = false;
  bool has_small_molecule_quantification_unit = false;
  bool has_assay_ms_run_ref = false;
  Size count_study_variable_assay_refs = 0;
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
  map<Size, Size> protein_abundance_study_variable_indices;
  map<Size, Size> protein_abundance_stdev_study_variable_indices;
  map<Size, Size> protein_abundance_std_err_study_variable_indices;

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
  Size peptide_uri_index = 0;
  Size peptide_spectra_ref_index = 0;
  map<Size, Size> peptide_abundance_assay_indices;
  map<Size, Size> peptide_abundance_study_variable_indices;
  map<Size, Size> peptide_abundance_stdev_study_variable_indices;
  map<Size, Size> peptide_abundance_std_err_study_variable_indices;

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
  Size psm_uri_index = 0;
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
  map<Size, Size> smallmolecule_abundance_std_err_study_variable_indices;

  for (TextFile::ConstIterator sit = tf.begin(); sit != tf.end(); ++sit)
  {
    //  std::cout << *sit << std::endl;
    const String & s = *sit;

    // skip empty lines or lines that are too short
    if (s.size() < 3)
    {
      continue;
    }

    cout << s << endl;

    const String section = s.prefix(3);

    // discard comments
    if (section == "COM")
    {
      continue;
    }

    StringList cells;
    s.split("\t", cells);

    if (cells.size() < 3)
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename, "Error parsing MzTab line: " + String(s) + ". Did you forget to use tabulator as separator?");
    }

    // parse metadata section
    if (section == "MTD")
    {
      StringList meta_key_fields; // the "-" separated fields of the metavalue key
      cells[1].split("-", meta_key_fields);
      String meta_key = meta_key_fields[0];

      if (cells[1].hasPrefix("mzTab-version"))
      {
        mz_tab_metadata.mz_tab_version.fromCellString(cells[2]);
        has_mzTab_version = true;
      }
      if (cells[1].hasPrefix("mzTab-mode"))
      {
        mz_tab_metadata.mz_tab_mode.fromCellString(cells[2]);
        has_mzTab_mode = true;
      }
      if (cells[1].hasPrefix("mzTab-type"))
      {
        mz_tab_metadata.mz_tab_type.fromCellString(cells[2]);
        has_mzTab_type = true;
      }
      if (cells[1].hasPrefix("mzTab-ID"))
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
        has_description = true;
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
        mz_tab_metadata.instrument_name[n] = p;
      }
      else if (meta_key.hasPrefix("instrument[") && meta_key_fields[1] == "source")
      {
        Int n = meta_key_fields[0].substitute("instrument[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.instrument_source[n] = p;
      }
      else if (meta_key.hasPrefix("instrument[") && meta_key_fields[1] == "analyzer")
      {
        Int n = meta_key_fields[0].substitute("instrument[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.instrument_analyzer[n] = p;
      }
      else if (meta_key.hasPrefix("instrument[") && meta_key_fields[1] == "detector")
      {
        Int n = meta_key_fields[0].substitute("instrument[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.instrument_detector[n] = p;
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
      else if (meta_key.hasPrefix("search_engine_score["))
      {
        Size n = (Size)meta_key_fields[0].substitute("search_engine_score[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.search_engine_score[n] = p;
        count_search_engine_score = std::max(n, count_search_engine_score); // will be checked to match number of entries in map to detect skipped entries or wrong numbering
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
        mz_tab_metadata.contact_name[n] = s;
      }
      else if (meta_key.hasPrefix("contact") && meta_key_fields[1] == "affiliation")
      {
        Int n = meta_key_fields[0].substitute("contact[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.contact_affiliation[n] = s;
      }
      else if (meta_key.hasPrefix("contact") && meta_key_fields[1] == "email")
      {
        Int n = meta_key_fields[0].substitute("contact[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.contact_email[n] = s;
      }
      else if (meta_key.hasPrefix("uri["))
      {
        Int n = meta_key_fields[0].substitute("uri[", "").substitute("]","").trim().toInt();
        MzTabString s;
        s.fromCellString(cells[2]);
        mz_tab_metadata.uri[n] = s;
      }
        //TODO: add mandatory check
      else if (meta_key.hasPrefix("variable_mod[" &&  meta_key_fields.size() == 1))
      {
        Int n = meta_key.substitute("variable_mod[", "").substitute("]","").trim().toInt();
        MzTabParameter pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.variable_mod[n] = pl;
      }
        //TODO: add mandatory check
      else if (meta_key.hasPrefix("fixed_mod[") &&  meta_key_fields.size() == 1) // fixed_mod[1-n]
      {
        Int n = meta_key.substitute("fixed_mod[", "").substitute("]","").trim().toInt();
        MzTabParameter pl;
        pl.fromCellString(cells[2]);
        mz_tab_metadata.fixed_mod[n] = pl;
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
        has_protein_quantification_unit = true;
      }
      else if (meta_key == "peptide" && meta_key_fields[1] == "quantification_unit")
      {
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.peptide_quantification_unit = p;
        has_peptide_quantification_unit = true;
      }
      else if (meta_key == "small_molecule" && meta_key_fields[1] == "quantification_unit")
      {
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.small_molecule_quantification_unit = p;
        has_small_molecule_quantification_unit = true;
      }
      else if (meta_key.hasPrefix("ms_run[") && meta_key_fields[1] == "format")
      {
        Int n = meta_key_fields[0].substitute("ms_run[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.ms_run_format[n] = p;
      }
      else if (meta_key.hasPrefix("ms_run[") && meta_key_fields[1] == "location")
      {
        Int n = meta_key_fields[0].substitute("ms_run[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.ms_run_location[n] = p;
        count_ms_run_location = std::max((Size)n, (Size)count_ms_run_location); // will be checked to match number of entries in map to detect skipped entries or wrong numbering
      }
      else if (meta_key.hasPrefix("ms_run[") && meta_key_fields[1] == "id_format")
      {
        Int n = meta_key_fields[0].substitute("ms_run[", "").substitute("]","").trim().toInt();
        MzTabParameter p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.ms_run_id_format[n] = p;
      }
      else if (meta_key.hasPrefix("ms_run[") && meta_key_fields[1] == "fragmentation_method")
      {
        Int n = meta_key_fields[0].substitute("ms_run[", "").substitute("]","").trim().toInt();
        MzTabParameterList p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.ms_run_fragmentation_method[n] = p;
      }
      else if (meta_key.hasPrefix("custom["))
      {
        Int n = meta_key.substitute("custom[", "").substitute("]","").trim().toInt();
        MzTabString sl;
        sl.fromCellString(cells[2]);
        mz_tab_metadata.publication[n] = sl;
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
        mz_tab_metadata.assay_quantification_reagent[n] = p;
      }
      else if (meta_key.hasPrefix("assay[") && meta_key_fields[1] == "sample_ref")
      {
        Int n = meta_key_fields[0].substitute("assay[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.assay_sample_ref[n] = p;
      }
      else if (meta_key.hasPrefix("cv[") && meta_key_fields[1] == "label")
      {
        Int n = meta_key_fields[0].substitute("cv[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.cv_label[n] = p;
      }
      else if (meta_key.hasPrefix("cv[") && meta_key_fields[1] == "full_name")
      {
        Int n = meta_key_fields[0].substitute("cv[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.cv_full_name[n] = p;
      }
      else if (meta_key.hasPrefix("cv[") && meta_key_fields[1] == "version")
      {
        Int n = meta_key_fields[0].substitute("cv[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.cv_version[n] = p;
      }
      else if (meta_key.hasPrefix("cv[") && meta_key_fields[1] == "url")
      {
        Int n = meta_key_fields[0].substitute("cv[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.cv_url[n] = p;
      }
      else if (meta_key.hasPrefix("assay[") && meta_key_fields[1] == "ms_run_ref")
      {
        Int n = meta_key_fields[0].substitute("assay[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.assay_ms_run_ref[n] = p;
        has_assay_ms_run_ref = true;
      }
      else if (meta_key.hasPrefix("study_variable[") && meta_key_fields[1] == "assay_refs")
      {
        Int n = meta_key_fields[0].substitute("study_variable[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        s.substitute("assay[","").substitute("]","");
        MzTabIntegerList p;
        p.fromCellString(s);
        mz_tab_metadata.study_variable_sample_refs[n] = p;
        count_study_variable_assay_refs++;
      }
      else if (meta_key.hasPrefix("study_variable[") && meta_key_fields[1] == "sample_refs")
      {
        Int n = meta_key_fields[0].substitute("study_variable[", "").substitute("]","").trim().toInt();
        String s = cells[2];
        s.substitute("sample[","").substitute("]","");
        MzTabIntegerList p;
        p.fromCellString(s);
        mz_tab_metadata.study_variable_sample_refs[n] = p;
      }
      else if (meta_key.hasPrefix("study_variable[") && meta_key_fields[1] == "description")
      {
        Int n = meta_key_fields[0].substitute("study_variable[", "").substitute("]","").trim().toInt();
        MzTabString p;
        p.fromCellString(cells[2]);
        mz_tab_metadata.study_variable_description[n] = p;
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
        std::cout << "#" << cells[i] << "#" << std::endl;
        if (cells[i] == "accession")
        {
          protein_accession_index = i;
        } else if (cells[i] == "description")
        {
          protein_description_index = i;
        } else if (cells[i] == "taxid")
        {
          protein_taxid_index = i;
        } else if (cells[i] == "species")
        {
          protein_species_index = i;
        } else if (cells[i] == "database")
        {
          protein_database_index = i;
        } else if (cells[i] == "database_version")
        {
          protein_database_version_index = i;
        } else if (cells[i] =="search_engine")
        {
          protein_search_engine_index = i;
        } else if (cells[i].hasPrefix("best_search_engine_score["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("best_search_engine_score[", "").substitute("]","").trim().toInt();
          protein_best_search_engine_score_to_column_index[n] = i;
        } else if (cells[i].hasPrefix("search_engine_score["))
        {
          std::pair<Size, Size> pair = extractIndexPairsFromBrackets_(cells[i].toQString());
          protein_column_index_to_score_runs_pair[i] = pair;
        } else if (cells[i] == "reliability")
        {
          protein_reliability_index = i;
        } else if (cells[i].hasPrefix("num_psms_ms_run["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("num_psms_ms_run[", "").substitute("]","").trim().toInt();
          protein_num_psms_ms_run_indices[n] = i;
        } else if (cells[i].hasPrefix("num_peptides_distinct_ms_run["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("num_peptides_distinct_ms_run[", "").substitute("]","").trim().toInt();
          protein_num_peptides_distinct_ms_run_indices[n] = i;
        } else if (cells[i].hasPrefix("num_peptides_unique_ms_run["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("num_peptides_unique_ms_run[", "").substitute("]","").trim().toInt();
          protein_num_peptides_unique_ms_run_indices[n] = i;
        } else if (cells[i] == "ambiguity_members")
        {
          protein_ambiguity_members_index = i;
        } else if (cells[i] == "modifications")
        {
          protein_modifications_index = i;
        } else if (cells[i] == "uri")
        {
          protein_uri_index = i;
        } else if (cells[i] == "go_terms")
        {
          protein_go_terms_index = i;
        } else if (cells[i] == "protein_coverage")
        {
          protein_coverage_index = i;
        } else if (cells[i].hasPrefix("protein_abundance_assay["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("protein_abundance_assay[", "").substitute("]","").trim().toInt();
          protein_abundance_assay_indices[n] = i;
        } else if (cells[i].hasPrefix("protein_abundance_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("protein_abundance_study_variable[", "").substitute("]","").trim().toInt();
          protein_abundance_stdev_study_variable_indices[n] = i;
        }  else if (cells[i].hasPrefix("protein_abundance_stdev_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("protein_abundance_stdev_study_variable[", "").substitute("]","").trim().toInt();
          protein_abundance_stdev_study_variable_indices[n] = i;
        }  else if (cells[i].hasPrefix("protein_abundance_std_err_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("protein_abundance_std_err_study_variable[", "").substitute("]","").trim().toInt();
          protein_abundance_std_err_study_variable_indices[n] = i;
        } else if (cells[i].hasPrefix("opt_"))
        {
          protein_custom_opt_columns[cells[i]] = i;
        }
      }
      continue;
    }

    // parse protein section
    if (section == "PRT")
    {
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

      if (protein_column_index_to_score_runs_pair.size() == 0 && mz_tab_metadata.mz_tab_mode.toCellString() == "Complete")
      {
        cout << "Error: mandatory protein search_engine_score_ms_run column(s) missing" << endl;
        // TODO: check if number of ms_runs and search_engine_scores match meta data section information
      }

      if (protein_num_psms_ms_run_indices.size() == 0 && mz_tab_metadata.mz_tab_mode.toCellString() == "Complete"
          && mz_tab_metadata.mz_tab_type.toCellString() == "Identification")
      {
        cout << "Error: mandatory protein num_psms_ms_run column(s) missing" << endl;
      }

      if (protein_num_peptides_distinct_ms_run_indices.size() == 0 && mz_tab_metadata.mz_tab_mode.toCellString() == "Complete"
          && mz_tab_metadata.mz_tab_type.toCellString() == "Identification")
      {
        cout << "Error: mandatory protein num_peptides_distinct_ms_run column(s) missing" << endl;
      }

      if (protein_num_peptides_unique_ms_run_indices.size() == 0 && mz_tab_metadata.mz_tab_mode.toCellString() == "Complete"
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

      if (protein_abundance_assay_indices.size() == 0 && mz_tab_metadata.mz_tab_mode.toCellString() == "Complete"
          && mz_tab_metadata.mz_tab_type.toCellString() == "Quantification")
      {
        cout << "Error: mandatory protein protein_abundance_assay column(s) missing" << endl;
      }

      if (protein_abundance_study_variable_indices.size() == 0 && mz_tab_metadata.mz_tab_type.toCellString() == "Quantification")
      {
        cout << "Error: mandatory protein_abundance_study_variable column(s) missing" << endl;
      }

      if (protein_abundance_stdev_study_variable_indices.size() == 0 && mz_tab_metadata.mz_tab_type.toCellString() == "Quantification")
      {
        cout << "Error: mandatory protein_abundance_stdev_study_variable column(s) missing" << endl;
      }

      if (protein_abundance_std_err_study_variable_indices.size() == 0 && mz_tab_metadata.mz_tab_type.toCellString() == "Quantification")
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

      row.protein_coverage.fromCellString(cells[protein_coverage_index]);

      // quantification data
      for (map<Size, Size>::const_iterator it = protein_abundance_assay_indices.begin(); it != protein_abundance_assay_indices.end(); ++it)
      {
        row.protein_abundance_assay[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = protein_abundance_study_variable_indices.begin(); it != protein_abundance_study_variable_indices.end(); ++it)
      {
        row.protein_abundance_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = protein_abundance_stdev_study_variable_indices.begin(); it != protein_abundance_stdev_study_variable_indices.end(); ++it)
      {
        row.protein_abundance_stdev_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = protein_abundance_std_err_study_variable_indices.begin(); it != protein_abundance_std_err_study_variable_indices.end(); ++it)
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
        } else if (cells[i] == "accession")
        {
          peptide_accession_index = i;
        } else if (cells[i] == "unique")
        {
          peptide_unique_index = i;
        } else if (cells[i] == "database")
        {
          peptide_database_index = i;
        } else if (cells[i] == "database_version")
        {
          peptide_database_version_index = i;
        } else if (cells[i] == "search_engine")
        {
          peptide_search_engine_index = i;
        } else if (cells[i].hasPrefix("best_search_engine_score["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("best_search_engine_score[", "").substitute("]","").trim().toInt();
          peptide_best_search_engine_score_to_column_index[n] = i;
        } else if (cells[i].hasPrefix("search_engine_score["))
        {
          std::pair<Size, Size> pair = extractIndexPairsFromBrackets_(cells[i].toQString());
          peptide_column_index_to_score_runs_pair[i] = pair;
        } else if (cells[i] == "reliability")
        {
          peptide_reliability_index = i;
        } else if (cells[i] == "modifications")
        {
          peptide_modifications_index = i;
        } else if (cells[i] == "retention_time")
        {
          peptide_retention_time_index = i;
        } else if (cells[i] == "retention_time_window")
        {
          peptide_retention_time_window_index = i;
        }  else if (cells[i] == "charge")
        {
          peptide_charge_index = i;
        }  else if (cells[i] == "mass_to_charge")
        {
          peptide_mass_to_charge_index = i;
        } else if (cells[i] == "uri")
        {
          protein_uri_index = i;
        } else if (cells[i] == "spectra_ref")
        {
          peptide_spectra_ref_index = i;
        }
        else if (cells[i].hasPrefix("peptide_abundance_assay["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("peptide_abundance_assay[", "").substitute("]","").trim().toInt();
          peptide_abundance_assay_indices[n] = i;
        } else if (cells[i].hasPrefix("peptide_abundance_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("peptide_abundance_study_variable[", "").substitute("]","").trim().toInt();
          peptide_abundance_stdev_study_variable_indices[n] = i;
        }  else if (cells[i].hasPrefix("peptide_abundance_stdev_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("peptide_abundance_stdev_study_variable[", "").substitute("]","").trim().toInt();
          peptide_abundance_stdev_study_variable_indices[n] = i;
        }  else if (cells[i].hasPrefix("peptide_abundance_std_err_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("peptide_abundance_std_err_study_variable[", "").substitute("]","").trim().toInt();
          peptide_abundance_std_err_study_variable_indices[n] = i;
        } else if (cells[i].hasPrefix("opt_"))
        {
          peptide_custom_opt_columns[cells[i]] = i;
        }
      }
      continue;
    }

    // parse peptide section
    if (section == "PEP")
    {
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

      for (map<Size, std::pair<Size, Size> >::const_iterator it = peptide_column_index_to_score_runs_pair.begin(); it != peptide_column_index_to_score_runs_pair.end(); ++it)
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

      if (peptide_uri_index != 0)
      {
        row.uri.fromCellString(cells[peptide_uri_index]);
      }

      row.spectra_ref.fromCellString(cells[peptide_spectra_ref_index]);

      // quantification data
      for (map<Size, Size>::const_iterator it = peptide_abundance_assay_indices.begin(); it != peptide_abundance_assay_indices.end(); ++it)
      {
        row.peptide_abundance_assay[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = peptide_abundance_study_variable_indices.begin(); it != peptide_abundance_study_variable_indices.end(); ++it)
      {
        row.peptide_abundance_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = peptide_abundance_stdev_study_variable_indices.begin(); it != peptide_abundance_stdev_study_variable_indices.end(); ++it)
      {
        row.peptide_abundance_stdev_study_variable[it->first].fromCellString(cells[it->second]);
      }

      for (map<Size, Size>::const_iterator it = peptide_abundance_std_err_study_variable_indices.begin(); it != peptide_abundance_std_err_study_variable_indices.end(); ++it)
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
        } else if (cells[i] == "PSM_ID")
        {
          psm_psm_id_index = i;
        } else if (cells[i] == "accession")
        {
          psm_accession_index = i;
        } else if (cells[i] == "unique")
        {
          psm_unique_index = i;
        } else if (cells[i] == "database")
        {
          psm_database_index = i;
        } else if (cells[i] == "database_version")
        {
          psm_database_version_index = i;
        } else if (cells[i] == "search_engine")
        {
          psm_search_engine_index = i;
        } else if (cells[i].hasPrefix("search_engine_score["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("search_engine_score[", "").substitute("]","").trim().toInt();
          psm_search_engine_score_to_column_index[n] = i;
        } else if (cells[i].hasPrefix("reliability"))
        {
          psm_reliability_index = i;
        } else if (cells[i] == "modifications")
        {
          psm_modifications_index = i;
        } else if (cells[i] == "retention_time")
        {
          psm_retention_time_index = i;
        } else if (cells[i] == "charge")
        {
          psm_charge_index = i;
        }  else if (cells[i] == "exp_mass_to_charge")
        {
          psm_exp_mass_to_charge_index = i;
        } else if (cells[i] == "calc_mass_to_charge")
        {
          psm_calc_mass_to_charge_index = i;
        } else if (cells[i] == "uri")
        {
          protein_uri_index = i;
        } else if (cells[i] == "spectra_ref")
        {
          psm_spectra_ref_index = i;
        } else if (cells[i] == "pre")
        {
          psm_pre_index = i;
        } else if (cells[i] == "post")
        {
          psm_post_index = i;
        } else if (cells[i] == "start")
        {
          psm_start_index = i;
        } else if (cells[i] == "end")
        {
          psm_end_index = i;
        }  else if (cells[i] == "opt_")
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

      if (psm_uri_index != 0)
      {
        row.uri.fromCellString(cells[psm_uri_index]);
      }

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
        if (cells[i].hasPrefix("identifier"))
        {
          smallmolecule_identifier_index = i;
        } else if (cells[i].hasPrefix("chemical_formula"))
        {
          smallmolecule_chemical_formula_index = i;
        } else if (cells[i].hasPrefix("smiles"))
        {
          smallmolecule_smiles_index = i;
        } else if (cells[i].hasPrefix("inchi_key"))
        {
          smallmolecule_inchi_key_index = i;
        } else if (cells[i].hasPrefix("description"))
        {
          smallmolecule_description_index = i;
        } else if (cells[i].hasPrefix("exp_mass_to_charge"))
        {
          smallmolecule_exp_mass_to_charge_index = i;
        } else if (cells[i].hasPrefix("calc_mass_to_charge"))
        {
          smallmolecule_calc_mass_to_charge_index = i;
        } else if (cells[i].hasPrefix("charge"))
        {
          smallmolecule_charge_index = i;
        } else if (cells[i].hasPrefix("retention_time"))
        {
          smallmolecule_retention_time_index = i;
        } else if (cells[i].hasPrefix("taxid"))
        {
          smallmolecule_taxid_index = i;
        } else if (cells[i].hasPrefix("species"))
        {
          smallmolecule_species_index = i;
        } else if (cells[i].hasPrefix("database"))
        {
          smallmolecule_database_index = i;
        } else if (cells[i].hasPrefix("database_version"))
        {
          smallmolecule_database_version_index = i;
        } else if (cells[i].hasPrefix("reliability"))
        {
          smallmolecule_reliability_index = i;
        } else if (cells[i].hasPrefix("uri"))
        {
          smallmolecule_uri_index = i;
        } else if (cells[i].hasPrefix("spectra_ref"))
        {
          smallmolecule_spectra_ref_index = i;
        } else if (cells[i].hasPrefix("search_engine"))
        {
          smallmolecule_search_engine_index = i;
        } else if (cells[i].hasPrefix("best_search_engine_score["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("best_search_engine_score[", "").substitute("]","").trim().toInt();
          smallmolecule_best_search_engine_score_to_column_index[n] = i;
        } else if (cells[i].hasPrefix("search_engine_score["))
        {
          std::pair<Size, Size> pair = extractIndexPairsFromBrackets_(cells[i].toQString());
          smallmolecule_column_index_to_score_runs_pair[i] = pair;
        } else if (cells[i].hasPrefix("modifications"))
        {
          smallmolecule_modifications_index = i;
        } else if (cells[i].hasPrefix("smallmolecule_abundance_assay["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("smallmolecule_abundance_assay[", "").substitute("]","").trim().toInt();
          smallmolecule_abundance_assay_indices[n] = i;
        } else if (cells[i].hasPrefix("smallmolecule_abundance_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("smallmolecule_abundance_study_variable[", "").substitute("]","").trim().toInt();
          smallmolecule_abundance_study_variable_indices[n] = i;
        }  else if (cells[i].hasPrefix("smallmolecule_abundance_stdev_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("smallmolecule_abundance_stdev_study_variable[", "").substitute("]","").trim().toInt();
          smallmolecule_abundance_stdev_study_variable_indices[n] = i;
        }  else if (cells[i].hasPrefix("smallmolecule_abundance_std_err_study_variable["))
        {
          String s = cells[i];
          Size n = (Size)s.substitute("smallmolecule_abundance_std_err_study_variable[", "").substitute("]","").trim().toInt();
          smallmolecule_abundance_std_err_study_variable_indices[n] = i;
        } else if (cells[i].hasPrefix("opt_"))
        {
          smallmolecule_custom_opt_columns[cells[i]] = i;
        }
      }
      continue;
    }

    // parse small molecule section
    if (section == "SML")
    {
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

      for (map<Size, Size>::const_iterator it = smallmolecule_abundance_std_err_study_variable_indices.begin(); it != smallmolecule_abundance_std_err_study_variable_indices.end(); ++it)
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
    mz_tab.setMetaData(mz_tab_metadata);
    mz_tab.setProteinSectionRows(mz_tab_protein_section_data);
    mz_tab.setPeptideSectionRows(mz_tab_peptide_section_data);
    mz_tab.setPSMSectionRows(mz_tab_psm_section_data);
    mz_tab.setSmallMoleculeSectionRows(mz_tab_small_molecule_section_data);
}

void MzTabFile::generateMzTabMetaDataSection_(const MzTabMetaData& md, StringList& sl) const
{
  sl.push_back(String("MTD\tmzTab-version\t") + md.mz_tab_version.toCellString());
  sl.push_back(String("MTD\tmzTab-mode\t") + md.mz_tab_mode.toCellString());
  sl.push_back(String("MTD\tmzTab-type\t") + md.mz_tab_type.toCellString());

  if (!md.mz_tab_id.isNull())
  {
    String s = String("MTD\tmzTab-ID\t") + md.mz_tab_id.toCellString();
    sl.push_back(s);
  }

  if (!md.title.isNull())
  {
    String s = String("MTD\ttitle\t") + md.title.toCellString();
    sl.push_back(s);
  }

  sl.push_back(String("MTD\tdescription\t") + md.description.toCellString());

  for (map<Size, MzTabParameterList>::const_iterator it = md.sample_processing.begin(); it != md.sample_processing.end(); ++it)
  {
    String s = "MTD\tsample_processing[" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.instrument_name.begin(); it != md.instrument_name.end(); ++it)
  {
    String s = "MTD\tinstrument[" + String(it->first) + "]-name\t" +  it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.instrument_source.begin(); it != md.instrument_source.end(); ++it)
  {
    String s = "MTD\tinstrument[" + String(it->first) + "]-source\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.instrument_analyzer.begin(); it != md.instrument_analyzer.end(); ++it)
  {
    String s = "MTD\tinstrument[" + String(it->first) + "]-analyzer\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.instrument_detector.begin(); it != md.instrument_detector.end(); ++it)
  {
    String s = "MTD\tinstrument[" + String(it->first) + "]-detector\t" + it->second.toCellString();
    sl.push_back(s);
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

  for (map<Size, MzTabParameter>::const_iterator it = md.search_engine_score.begin(); it != md.search_engine_score.end(); ++it)
  {
    String s = "MTD\tsearch_engine_score[\t" + String(it->first) + "]\t" + it->second.toCellString();
    sl.push_back(s);
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

  for (map<Size, MzTabString>::const_iterator it = md.contact_name.begin(); it != md.contact_name.end(); ++it)
  {
    String s = "MTD\tcontact[" + String(it->first) + "]-name\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.contact_affiliation.begin(); it != md.contact_affiliation.end(); ++it)
  {
    String s = "MTD\tcontact[" + String(it->first) + "]-affiliation\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.contact_email.begin(); it != md.contact_email.end(); ++it)
  {
    String s = "MTD\tcontact[" + String(it->first) + "]-email\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.uri.begin(); it != md.uri.end(); ++it)
  {
    String s = "MTD\turi[" + String(it->first) + String("]\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.fixed_mod.begin(); it != md.fixed_mod.end(); ++it)
  {
    String s = "MTD\tfixed_mod[" + String(it->first) + String("]\t")+ it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.fixed_mod_site.begin(); it != md.fixed_mod_site.end(); ++it)
  {
    String s = "MTD\tfixed_mod[" + String(it->first) + String("]-site\t")+ it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.fixed_mod_position.begin(); it != md.fixed_mod_position.end(); ++it)
  {
    String s = "MTD\tfixed_mod[" + String(it->first) + String("]-position\t")+ it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabParameter>::const_iterator it = md.variable_mod.begin(); it != md.variable_mod.end(); ++it)
  {
    String s = "MTD\tvariable_mod[" + String(it->first) + String("]\t")+ it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.variable_mod_site.begin(); it != md.variable_mod_site.end(); ++it)
  {
    String s = "MTD\tvariable_mod[" + String(it->first) + String("]-site\t")+ it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.variable_mod_position.begin(); it != md.variable_mod_position.end(); ++it)
  {
    String s = "MTD\tvariable_mod[" + String(it->first) + String("]-position\t")+ it->second.toCellString();
    sl.push_back(s);
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

  // ms_run[1-n]-format
  for (map<Size, MzTabParameter>::const_iterator it = md.ms_run_format.begin(); it != md.ms_run_format.end(); ++it)
  {
    String s = "MTD\tms_run[" + String(it->first) + "]-format\t" + it->second.toCellString();
    sl.push_back(s);
  }

  // ms_run[1-n]-location
  for (map<Size, MzTabString>::const_iterator it = md.ms_run_location.begin(); it != md.ms_run_location.end(); ++it)
  {
    String s = "MTD\tms_run[" + String(it->first) + "]-location\t" + it->second.toCellString();
    sl.push_back(s);
  }

  // ms_run[1-n]-id_format
  for (map<Size, MzTabParameter>::const_iterator it = md.ms_run_id_format.begin(); it != md.ms_run_id_format.end(); ++it)
  {
    String s = "MTD\tms_run[" + String(it->first) + "]-id_format\t" + it->second.toCellString();
    sl.push_back(s);
  }

  // ms_run[1-n]-fragmentation_method
  for (map<Size, MzTabParameter>::const_iterator it = md.ms_run_id_format.begin(); it != md.ms_run_id_format.end(); ++it)
  {
    String s = "MTD\tms_run[" + String(it->first) + "]-fragmentation_method\t" + it->second.toCellString();
    sl.push_back(s);
  }

  // custom
  for (map<Size, MzTabParameterList>::const_iterator it = md.custom.begin(); it != md.custom.end(); ++it)
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


  for (map<Size, MzTabParameter>::const_iterator it = md.assay_quantification_reagent.begin(); it != md.assay_quantification_reagent.end(); ++it)
  {
    String s = "MTD\tassay[" + String(it->first) + "]-quantification_reagent\t" + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, map<Size, MzTabParameter> >::const_iterator it = md.assay_quantification_mod.begin(); it != md.assay_quantification_mod.end(); ++it)
  {
    for (map<Size, MzTabParameter>::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
    {
      String s = "MTD\tassay[" + String(it->first) + String("]-quantification_mod[") + String(jt->first) + String("]\t") + jt->second.toCellString();
      sl.push_back(s);
    }
  }

  for (map<Size, map<Size, MzTabString> >::const_iterator it = md.assay_quantification_mod_site.begin(); it != md.assay_quantification_mod_site.end(); ++it)
  {
    for (map<Size, MzTabString>::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
    {
      String s = "MTD\tassay[" + String(it->first) + String("]-quantification_mod[") + String(jt->first) + String("]-site\t") + jt->second.toCellString();
      sl.push_back(s);
    }
  }

  for (map<Size, map<Size, MzTabString> >::const_iterator it = md.assay_quantification_mod_position.begin(); it != md.assay_quantification_mod_position.end(); ++it)
  {
    for (map<Size, MzTabString>::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
    {
      String s = "MTD\tassay[" + String(it->first) + String("]-quantification_mod[") + String(jt->first) + String("]-position\t") + jt->second.toCellString();
      sl.push_back(s);
    }
  }

  for (map<Size, MzTabString>::const_iterator it = md.assay_sample_ref.begin(); it != md.assay_sample_ref.end(); ++it)
  {
    String s = "MTD\tassay[" + String(it->first) + String("]-sample_ref\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.assay_ms_run_ref.begin(); it != md.assay_ms_run_ref.end(); ++it)
  {
    String s = "MTD\tassay[" + String(it->first) + String("]-ms_run_ref\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabIntegerList>::const_iterator it = md.study_variable_assay_refs.begin(); it != md.study_variable_assay_refs.end(); ++it)
  {
    String s = "MTD\tstudy_variable[" + String(it->first) + "]-assay_refs\t" + it->second.toCellString();
    sl.push_back(s);
  }


  for (map<Size, MzTabIntegerList>::const_iterator it = md.study_variable_sample_refs.begin(); it != md.study_variable_sample_refs.end(); ++it)
  {
    String s = "MTD\tstudy_variable[" + String(it->first) + String("]-sample_refs\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.study_variable_description.begin(); it != md.study_variable_description.end(); ++it)
  {
    String s = "MTD\tstudy_variable[" + String(it->first) + String("]-descriptiont\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.cv_label.begin(); it != md.cv_label.end(); ++it)
  {
    String s = "MTD\tcv[" + String(it->first) + String("]-label\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.cv_full_name.begin(); it != md.cv_full_name.end(); ++it)
  {
    String s = "MTD\tcv[" + String(it->first) + String("]-full_name\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.cv_version.begin(); it != md.cv_version.end(); ++it)
  {
    String s = "MTD\tcv[" + String(it->first) + String("]-version\t") + it->second.toCellString();
    sl.push_back(s);
  }

  for (map<Size, MzTabString>::const_iterator it = md.cv_url.begin(); it != md.cv_url.end(); ++it)
  {
    String s = "MTD\tcv[" + String(it->first) + String("]-url\t") + it->second.toCellString();
    sl.push_back(s);
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
  sl.push_back(String("\n"));
}

String MzTabFile::generateMzTabProteinHeader_(Size search_ms_runs, Size num_psms_ms_runs, Size num_peptides_distinct_ms_runs, Size num_peptides_unique_ms_run, Size assays, Size study_variables, const std::vector<String>& optional_columns) const
{
  StringList header;
  header.push_back("PRH");
  header.push_back("accession");
  header.push_back("description");
  header.push_back("taxid");
  header.push_back("species");
  header.push_back("database");
  header.push_back("database_version");
  header.push_back("search_engine");
  header.push_back("best_search_engine_score");

  for (Size i = 0; i != search_ms_runs; ++i)
  {
    header.push_back(String("search_engine_score_ms_run[") + String(i) + String("]"));
  }

  header.push_back("reliability");

  for (Size i = 0; i != num_psms_ms_runs; ++i)
  {
    header.push_back(String("num_psms_ms_run[") + String(i) + String("]"));
  }

  for (Size i = 0; i != num_peptides_distinct_ms_runs; ++i)
  {
    header.push_back(String("num_peptides_distinct_ms_run[") + String(i) + String("]"));
  }

  for (Size i = 0; i != num_peptides_unique_ms_run; ++i)
  {
    header.push_back(String("num_peptides_unique_ms_run[") + String(i) + String("]"));
  }

  header.push_back("ambiguity_members");
  header.push_back("modifications");
  header.push_back("uri");
  header.push_back("go_terms");
  header.push_back("protein_coverage");

  for (Size i = 0; i != assays; ++i)
  {
    header.push_back(String("protein_abundance_assay[") + String(i) + String("]"));
  }

  for (Size i = 0; i != study_variables; ++i)
  {
    header.push_back(String("protein_abundance_study_variable[") + String(i) + String("]"));
    header.push_back(String("protein_abundance_stdev_study_variable[") + String(i) + String("]"));
    header.push_back(String("protein_abundance_std_error_study_variable[") + String(i) + String("]"));
  }

  std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));

  return ListUtils::concatenate(header, "\t");
}

String MzTabFile::generateMzTabProteinSectionRow_(const MzTabProteinSectionRow& row) const
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

  s.push_back(row.reliability.toCellString());

  for (std::map<Size, MzTabInteger>::const_iterator it = row.num_psms_ms_run.begin(); it != row.num_psms_ms_run.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  for (std::map<Size, MzTabInteger>::const_iterator it = row.num_peptides_distinct_ms_run.begin(); it != row.num_peptides_distinct_ms_run.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  String MzTabFile::generateMzTabProteinHeader_(Int n_subsamples, const std::vector<String>& optional_protein_columns) const
  for (std::map<Size, MzTabInteger>::const_iterator it = row.num_peptides_unique_ms_run.begin(); it != row.num_peptides_unique_ms_run.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  s.push_back(row.ambiguity_members.toCellString());
  s.push_back(row.modifications.toCellString());
  s.push_back(row.uri.toCellString());
  s.push_back(row.go_terms.toCellString());
  s.push_back(row.protein_coverage.toCellString());

  // quantification columns
  for (std::map<Size, MzTabDouble>::const_iterator it = row.protein_abundance_assay.begin(); it != row.protein_abundance_assay.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  std::map<Size, MzTabDouble>::const_iterator sv_it = row.protein_abundance_study_variable.begin();
  std::map<Size, MzTabDouble>::const_iterator sv_stdev_it = row.protein_abundance_stdev_study_variable.begin();
  std::map<Size, MzTabDouble>::const_iterator sv_error_it = row.protein_abundance_std_error_study_variable.begin();

  for (;
       sv_it != row.protein_abundance_study_variable.end()
       && sv_stdev_it != row.protein_abundance_stdev_study_variable.end()
       && sv_error_it != row.protein_abundance_std_error_study_variable.end();
       ++sv_it, ++sv_stdev_it, ++sv_error_it)
  {
    s.push_back(sv_it->second.toCellString());
    s.push_back(sv_stdev_it->second.toCellString());
    s.push_back(sv_error_it->second.toCellString());
  }

  // print optional columns
  for (Size i = 0; i != row.opt_.size(); ++i)
  {
    s.push_back(row.opt_[i].second.toCellString());
  }

  return ListUtils::concatenate(s, "\t");
}

void MzTabFile::generateMzTabProteinSection_(const MzTabProteinSectionRows& rows, StringList& sl) const
{
  for (MzTabProteinSectionRows::const_iterator it = rows.begin(); it != rows.end(); ++it)
  {
    sl.push_back(generateMzTabProteinSectionRow_(*it));
  }
  sl.push_back(String("\n"));
}

void MzTabFile::generateMzTabPeptideSection_(const MzTabPeptideSectionRows& rows, StringList& sl) const
{
  for (MzTabPeptideSectionRows::const_iterator it = rows.begin(); it != rows.end(); ++it)
  {
    sl.push_back(generateMzTabPeptideSectionRow_(*it));
  }
  sl.push_back(String("\n"));
}

void MzTabFile::generateMzTabSmallMoleculeSection_(const MzTabSmallMoleculeSectionRows& rows, StringList& sl) const
{
  for (MzTabSmallMoleculeSectionRows::const_iterator it = rows.begin(); it != rows.end(); ++it)
  {
    sl.push_back(generateMzTabSmallMoleculeSectionRow_(*it));
  }
}

String MzTabFile::generateMzTabPeptideHeader_(Size search_ms_runs, Size assays, Size study_variables, const vector<String>& optional_columns) const
{
  StringList header;
  header.push_back("PEH");
  header.push_back("sequence");
  header.push_back("accession");
  header.push_back("unique");
  header.push_back("database");
  header.push_back("database_version");
  header.push_back("search_engine");
  header.push_back("best_search_engine_score");

  for (Size i = 0; i != search_ms_runs; ++i)
  {
    header.push_back(String("search_engine_score_ms_run[") + String(i) + String("]"));
  }

  header.push_back("reliability");
  header.push_back("modifications");
  header.push_back("retention_time");
  header.push_back("retention_time_window");
  header.push_back("charge");
  header.push_back("mass_to_charge");
  header.push_back("uri");
  header.push_back("spectra_ref");

  for (Size i = 0; i != assays; ++i)
  {
    header.push_back(String("peptide_abundance_assay[") + String(i) + String("]"));
  }

  for (Size i = 0; i != study_variables; ++i)
  {
    header.push_back(String("peptide_abundance_study_variable[") + String(i) + String("]"));
    header.push_back(String("peptide_abundance_stdev_study_variable[") + String(i) + String("]"));
    header.push_back(String("peptide_abundance_std_error_study_variable[") + String(i) + String("]"));
  }

  std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));

  return ListUtils::concatenate(header, "\t");
}

String MzTabFile::generateMzTabPSMHeader_(const vector<String>& optional_columns) const
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
  header.push_back("search_engine_score");
  header.push_back("reliability");
  header.push_back("modifications");
  header.push_back("retention_time");
  header.push_back("charge");
  header.push_back("exp_mass_to_charge");
  header.push_back("calc_mass_to_charge");
  header.push_back("uri");
  header.push_back("spectra_ref");
  header.push_back("pre");
  header.push_back("post");
  header.push_back("start");
  header.push_back("end");

  std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));

  return ListUtils::concatenate(header, "\t");
}

String MzTabFile::generateMzTabPeptideSectionRow_(const MzTabPeptideSectionRow& row) const
{
  StringList s;
  s.push_back("PEP");
  s.push_back(row.sequence.toCellString());
  s.push_back(row.accession.toCellString());
  s.push_back(row.unique.toCellString());
  s.push_back(row.database.toCellString());
  s.push_back(row.database_version.toCellString());
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

  s.push_back(row.reliability.toCellString());
  s.push_back(row.modifications.toCellString());
  s.push_back(row.retention_time.toCellString());
  s.push_back(row.charge.toCellString());
  s.push_back(row.mass_to_charge.toCellString());
  s.push_back(row.uri.toCellString());
  s.push_back(row.spectra_ref.toCellString());

  // quantification columns
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

  // print optional columns
  for (Size i = 0; i != row.opt_.size(); ++i)
  {
    s.push_back(row.opt_[i].second.toCellString());
  }

  return ListUtils::concatenate(s, "\t");
}

void MzTabFile::generateMzTabPSMSection_(const MzTabPSMSectionRows& rows, StringList& sl) const
{
  for (MzTabPSMSectionRows::const_iterator it = rows.begin(); it != rows.end(); ++it)
  {
    sl.push_back(generateMzTabPSMSectionRow_(*it));
  }
  sl.push_back(String("\n"));
}

String MzTabFile::generateMzTabPSMSectionRow_(const MzTabPSMSectionRow& row) const
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

  for (map<Size, MzTabDouble>::const_iterator it = row.search_engine_score.begin(); it != row.search_engine_score.end(); ++it)
  {
    s.push_back(it->second.toCellString());
  }

  s.push_back(row.reliability.toCellString());
  s.push_back(row.modifications.toCellString());
  s.push_back(row.retention_time.toCellString());
  s.push_back(row.charge.toCellString());
  s.push_back(row.exp_mass_to_charge.toCellString());
  s.push_back(row.calc_mass_to_charge.toCellString());
  s.push_back(row.uri.toCellString());
  s.push_back(row.spectra_ref.toCellString());
  s.push_back(row.pre.toCellString());
  s.push_back(row.post.toCellString());
  s.push_back(row.start.toCellString());
  s.push_back(row.end.toCellString());

  // print optional columns
  for (Size i = 0; i != row.opt_.size(); ++i)
  {
    s.push_back(row.opt_[i].second.toCellString());
  }

  return ListUtils::concatenate(s, "\t");
}

String MzTabFile::generateMzTabSmallMoleculeHeader_(Size ms_runs, Size assays, Size study_variables, const vector<String>& optional_smallmolecule_columns) const
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
  header.push_back("reliability");
  header.push_back("uri");
  header.push_back("spectra_ref");
  header.push_back("search_engine");
  header.push_back("best_search_engine_score");

  for (Size i = 0; i != ms_runs; ++i)
  {
    header.push_back(String("search_engine_score_ms_run[") + String(i) + String("]"));
  }

  header.push_back("modifications");

  for (Size i = 0; i != assays; ++i)
  {
    header.push_back(String("smallmolecule_abundance_assay[") + String(i) + String("]"));
  }

  for (Size i = 0; i != study_variables; ++i)
  {
    header.push_back(String("smallmolecule_abundance_study_variable[") + String(i) + String("]"));
    header.push_back(String("smallmolecule_abundance_stdev_study_variable[") + String(i) + String("]"));
    header.push_back(String("smallmolecule_abundance_std_error_study_variable[") + String(i) + String("]"));
  }

  // copy optional column names to header
  std::copy(optional_smallmolecule_columns.begin(), optional_smallmolecule_columns.end(), std::back_inserter(header));

  return ListUtils::concatenate(header, "\t");
}

String MzTabFile::generateMzTabSmallMoleculeSectionRow_(const MzTabSmallMoleculeSectionRow& row) const
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
  s.push_back(row.reliability.toCellString());
  s.push_back(row.uri.toCellString());
  s.push_back(row.spectra_ref.toCellString());
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

  // print optional columns
  for (Size i = 0; i != row.opt_.size(); ++i)
  {
    s.push_back(row.opt_[i].second.toCellString());
  }

  return ListUtils::concatenate(s, "\t");
}

void MzTabFile::store(const String& filename, const MzTab& mz_tab) const
{
  TextFile out;

  generateMzTabMetaDataSection_(mz_tab.getMetaData(), out);

  bool complete = (mz_tab.getMetaData().mz_tab_mode.toCellString() == "Complete");
  bool quantification = (mz_tab.getMetaData().mz_tab_type.toCellString() == "Quantification");
  Size ms_runs = mz_tab.getMetaData().ms_run_location.size();

  const MzTabProteinSectionRows& protein_section = mz_tab.getProteinSectionRows();
  const MzTabPeptideSectionRows& peptide_section = mz_tab.getPeptideSectionRows();
  const MzTabPSMSectionRows& psm_section = mz_tab.getPSMSectionRows();
  const MzTabSmallMoleculeSectionRows& smallmolecule_section = mz_tab.getSmallMoleculeSectionRows();

  if (!protein_section.empty())
  {
    Size search_ms_runs;
    Size num_psms_ms_runs;
    Size num_peptides_distinct_ms_runs;
    Size num_peptides_unique_ms_run;
    Size assays = protein_section[0].protein_abundance_assay.size();
    Size study_variables = protein_section[0].protein_abundance_study_variable.size();;

    if (complete)
    {
      // all ms_runs mandatory
      search_ms_runs = ms_runs;

      if (!quantification) // mandatory for complete identification files
      {
        num_psms_ms_runs = ms_runs;
        num_peptides_distinct_ms_runs = ms_runs;
        num_peptides_unique_ms_run = ms_runs;
      } else
      {
        num_psms_ms_runs = protein_section[0].num_psms_ms_run.size();
        num_peptides_distinct_ms_runs = protein_section[0].num_peptides_distinct_ms_run.size();
        num_peptides_unique_ms_run = protein_section[0].num_peptides_unique_ms_run.size();
      }
    } else // only user provided ones are reported
    {
      search_ms_runs = protein_section[0].search_engine_score_ms_run.size();
      num_psms_ms_runs = protein_section[0].num_psms_ms_run.size();
      num_peptides_distinct_ms_runs = protein_section[0].num_peptides_distinct_ms_run.size();
      num_peptides_unique_ms_run = protein_section[0].num_peptides_unique_ms_run.size();
    }

    out.push_back(generateMzTabProteinHeader_(ms_runs, num_psms_ms_runs, num_peptides_distinct_ms_runs, num_peptides_unique_ms_run, assays, study_variables, mz_tab.getProteinOptionalColumnNames()));
    generateMzTabProteinSection_(protein_section, out);
  }

  if (!peptide_section.empty())
  {
    Size assays = peptide_section[0].peptide_abundance_assay.size();
    Size study_variables = peptide_section[0].peptide_abundance_study_variable.size();
    Size search_ms_runs;
    if (complete)
    {
      // all ms_runs mandatory
      search_ms_runs = ms_runs;
    } else // only user provided ones are reported
    {
      search_ms_runs = mz_tab.getPeptideSectionRows()[0].search_engine_score_ms_run.size();
    }

    out.push_back(generateMzTabPeptideHeader_(search_ms_runs, assays, study_variables, mz_tab.getPeptideOptionalColumnNames()));
    generateMzTabPeptideSection_(mz_tab.getPeptideSectionRows(), out);
  }

  if (!psm_section.empty())
  {
    out.push_back(generateMzTabPSMHeader_(mz_tab.getPSMOptionalColumnNames()));
    generateMzTabPSMSection_(mz_tab.getPSMSectionRows(), out);
  }

  if (!smallmolecule_section.empty())
  {
    Size assays = smallmolecule_section[0].smallmolecule_abundance_assay.size();
    Size study_variables = smallmolecule_section[0].smallmolecule_abundance_study_variable.size();
    Size search_ms_runs;
    if (complete)
    {
      // all ms_runs mandatory
      search_ms_runs = ms_runs;
    } else // only user provided ones are reported
    {
      search_ms_runs = mz_tab.getPeptideSectionRows()[0].search_engine_score_ms_run.size();
    }

    out.push_back(generateMzTabSmallMoleculeHeader_(ms_runs, assays, study_variables, mz_tab.getSmallMoleculeOptionalColumnNames()));
    generateMzTabSmallMoleculeSection_(mz_tab.getSmallMoleculeSectionRows(), out);
  }

  out.store(filename);
}

void MzTabFile::store(const String& filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids, String in, String document_id) const
{
  // copy data as we are going to apply some filters
  vector<ProteinIdentification> prot_ids = protein_ids;
  vector<PeptideIdentification> pep_ids = peptide_ids;

  // create tab separated output stream for MzTab file
  ofstream txt_out(filename.c_str());
  SVOutStream output(txt_out, "\t", "_", String::NONE);

  // pre-filter for best PSM
  sortPSM_(pep_ids.begin(), pep_ids.end());
  keepFirstPSM_(pep_ids.begin(), pep_ids.end());

  // warn if no coverage information is available
  bool has_coverage = true;
  try // might throw Exception::MissingInformation() if no protein sequence information is added
  {
    for (Size i = 0; i < prot_ids.size(); ++i)
    {
      prot_ids[i].computeCoverage(pep_ids);
    }
  }
  catch (Exception::MissingInformation& e)
  {
    cout << e.what() << "\n";
    has_coverage = false;
  }

  // partition into runs (as this could be a merged idXML)
  map<String, vector<PeptideIdentification> > map_run_to_pepids;
  map<String, vector<ProteinIdentification> > map_run_to_proids;
  partitionIntoRuns_(pep_ids, prot_ids, map_run_to_pepids, map_run_to_proids);
  // cout << "idXML contains: " << map_run_to_pepids.size() << " runs." << endl;

  MapAccPepType map_run_accesion_to_peptides;
  createProteinToPeptideLinks_(map_run_to_pepids, map_run_accesion_to_peptides);

  // every ProteinIdentification corresponds to a search engine run and contains protein hits

  // write meta data for each run
  bool meta_info_printed = false;
  Size run_count = 0;
  map<String, vector<ProteinIdentification> >::const_iterator mprot_it = map_run_to_proids.begin();
  map<String, vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();

  for (; mprot_it != map_run_to_proids.end(); ++mprot_it) // iterate over runs
  {
    // extract ProteinIdentifications of this run (should be only 1)
    const vector<ProteinIdentification>& prot_ids = mprot_it->second;
    for (vector<ProteinIdentification>::const_iterator it = prot_ids.begin(); it != prot_ids.end(); ++it)
    {
      String UNIT_ID = "OpenMS_" + String(run_count);
      String title = document_id;
      if (title != "")
      {
        output << "MOD" << UNIT_ID + "-title" << title;
        meta_info_printed = true;
      }
    }
    run_count++;
  }

  if (meta_info_printed)
  {
    output << endl;
  }

  // write protein table header
  if (meta_info_printed)
  {
    output << endl;
  }

  // determine the number of sub samples in each run from protein ids (it is assumed that peptide ids don't introduce new sub sample categories)
  map<String, Size> map_run_to_n_subsamples = extractNumberOfSubSamples_(map_run_to_proids);
  writeProteinHeader_(output, map_run_to_n_subsamples);

  // write protein table data
  run_count = 0;
  mprot_it = map_run_to_proids.begin();
  mpep_it = map_run_to_pepids.begin();
  for (; mprot_it != map_run_to_proids.end(); ++mprot_it, ++mpep_it) // iterate over runs
  {
    // extract ProteinIdentifications of this run (should be only 1)
    const vector<ProteinIdentification>& prot_ids = mprot_it->second;
    for (vector<ProteinIdentification>::const_iterator prot_id_it = prot_ids.begin();
         prot_id_it != prot_ids.end(); ++prot_id_it, ++run_count)
    {
      writeProteinData_(output, *prot_id_it, run_count, in, has_coverage, map_run_accesion_to_peptides, map_run_to_n_subsamples);
    }
  }

  // write peptide header
  output << endl;

  writePeptideHeader_(output, map_run_to_n_subsamples);

  run_count = 0;
  mprot_it = map_run_to_proids.begin();
  mpep_it = map_run_to_pepids.begin();
  for (; mpep_it != map_run_to_pepids.end(); ++mprot_it, ++mpep_it, ++run_count) // iterate over runs
  {
    // extract ProteinIdentifications of this run (should be only 1)
    const vector<ProteinIdentification>& prot_ids = mprot_it->second;
    // extract PeptideIdentifications of this run
    const vector<PeptideIdentification>& pep_ids = mpep_it->second;
    // iterate over runs of peptide identifications
    for (vector<PeptideIdentification>::const_iterator pep_id_it = pep_ids.begin();
         pep_id_it != pep_ids.end(); ++pep_id_it)
    {
      // TODO: check if bad design of Protein/PeptideIdentification as search engine parameters are stored in prot.
      String openms_search_engine_name = prot_ids[0].getSearchEngine();
      String search_engine_cvParams = mapSearchEngineToCvParam_(openms_search_engine_name);

      const ProteinIdentification::SearchParameters& sp = prot_ids[0].getSearchParameters();
      String UNIT_ID_String = "OpenMS_" + String(run_count);
      String database_String = (sp.db != "" ? sp.db : "--");
      database_String = "file://" + database_String;
      String database_version_String = (sp.db_version != "" ? sp.db_version : "--");

      for (vector<PeptideHit>::const_iterator peptide_hit_it = pep_id_it->getHits().begin();
           peptide_hit_it != pep_id_it->getHits().end(); ++peptide_hit_it)
      {
        String sequence = peptide_hit_it->getSequence().toString();
        String accession = extractProteinAccession_(*peptide_hit_it);
        String unit_id = UNIT_ID_String;
        String unique;
        String database = database_String;
        String database_version = database_version_String;
        String search_engine = search_engine_cvParams;
        String search_engine_score = mapSearchEngineScoreToCvParam_(openms_search_engine_name,
                                                                    peptide_hit_it->getScore(),
                                                                    pep_id_it->getScoreType());
        String modifications = extractPeptideModifications_(*peptide_hit_it); //TODO: check if terminal mods work

        String spectra_ref = "--";

        if (peptide_hit_it->metaValueExists("mzTab:unique"))
        {
          bool is_unique = peptide_hit_it->getMetaValue("mzTab:unique").toBool();
          if (is_unique)
          {
            unique = "1";
          }
          else
          {
            unique = "0";
          }
        }
        else // no uniqueness annotation
        {
          // if unique protein is present peptide can be assigned
          if (peptide_hit_it->getProteinAccessions().size() == 1)
          {
            unique = "1";
          }
          else
          {
            unique = "0";
          }
        }

        String retention_time;
        if (pep_id_it->metaValueExists("RT")) // Note: RT stored on pep_id_it not on hit
        {
          retention_time = String::number(String(pep_id_it->getMetaValue("RT")).toDouble(), 2);
        }
        else
        {
          retention_time = "--";
        }

        String mass_to_charge;
        if (pep_id_it->metaValueExists("MZ")) // Note: MZ stored on pep_id_it not on hit
        {
          mass_to_charge = String::number(String(pep_id_it->getMetaValue("MZ")).toDouble(), 10);
        }
        else
        {
          mass_to_charge = "--";
        }

        String charge = "NA";
        if (peptide_hit_it->getCharge() != 0)
        {
          charge = peptide_hit_it->getCharge();
        }

        String uri = "file://" + in;

        output << "PEP" << sequence << accession << unit_id << unique << database
               << database_version << search_engine << search_engine_score
               << modifications << retention_time << charge
               << mass_to_charge << uri << spectra_ref;

        // get number of sub samples for this run
        map<String, Size>::const_iterator sub_it = map_run_to_n_subsamples.find(mpep_it->first);

        Size n_subsamples = 0;
        if (sub_it != map_run_to_n_subsamples.end())
        {
          n_subsamples = sub_it->second;
        }

        for (Size n = 1; n <= n_subsamples; ++n)
        {
          {
            String key = "mzTab:peptide_abundance_sub[" + String(n) + "]";
            String abundancy_value = "--";
            if (peptide_hit_it->metaValueExists(key))
            {
              abundancy_value = peptide_hit_it->getMetaValue(key);
            }
            output << abundancy_value;
          }
          {
            String key = "mzTab:peptide_abundance_stdev_sub[" + String(n) + "]";
            String abundancy_value = "--";
            if (peptide_hit_it->metaValueExists(key))
            {
              abundancy_value = peptide_hit_it->getMetaValue(key);
            }
            output << abundancy_value;
          }
          {
            String key = "mzTab:peptide_abundance_std_error_sub[" + String(n) + "]";
            String abundancy_value = "--";
            if (peptide_hit_it->metaValueExists(key))
            {
              abundancy_value = peptide_hit_it->getMetaValue(key);
            }
            output << abundancy_value;
          }
        }

        output << endl;

      }
    }
  }
  txt_out.close();
}

/// Extract protein and peptide identifications for each run. maps are assumed empty.
void MzTabFile::partitionIntoRuns_(const vector<PeptideIdentification>& pep_ids,
                                   const vector<ProteinIdentification>& pro_ids,
                                   map<String, vector<PeptideIdentification> >& map_run_to_pepids,
                                   map<String, vector<ProteinIdentification> >& map_run_to_proids
                                   )
{
  {
    // idXML can be merged so we have to deal with several runs
    // Extract protein and peptide identifications for each run
    for (Size i = 0; i != pro_ids.size(); ++i)
    {
      map_run_to_proids[pro_ids[i].getIdentifier()].push_back(pro_ids[i]);
    }

    for (Size i = 0; i != pep_ids.size(); ++i)
    {
      map_run_to_pepids[pep_ids[i].getIdentifier()].push_back(pep_ids[i]);
    }

      // perform some sanity check (run ids should be identical)
      assert(map_run_to_proids.size() == map_run_to_pepids.size());
      std::map<String, std::vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();
      std::map<String, std::vector<ProteinIdentification> >::const_iterator mprot_it = map_run_to_proids.begin();
      for (; mpep_it != map_run_to_pepids.end(); ++mpep_it, ++mprot_it)
      {
        assert(mpep_it->first == mprot_it->first);
      }
    }
  }
}

void MzTabFile::createProteinToPeptideLinks_(const map<String, vector<PeptideIdentification> >& map_run_to_pepids, MapAccPepType& map_run_accession_to_pephits)
{
  // create links for each run
  map<String, vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();

  for (; mpep_it != map_run_to_pepids.end(); ++mpep_it)
  {
    const String& run = mpep_it->first;
    const vector<PeptideIdentification>& pids = mpep_it->second;
    for (Size i = 0; i != pids.size(); ++i)
    {
      const std::vector<PeptideHit>& phits = pids[i].getHits();
      if (phits.size() == 1)
      {
        const std::vector<String>& accessions = phits[0].getProteinAccessions();
        for (Size k = 0; k != accessions.size(); ++k)
        {
          const String& accession = accessions[k];
          std::pair<String, String> key = make_pair(run, accession);
          map_run_accession_to_pephits[key].push_back(phits[0]);
        }
      }
    }
  }
}

/// Extracts, if possible a unique protein accession for a peptide hit in mzTab format. Otherwise NA is returned
String MzTabFile::extractProteinAccession_(const PeptideHit& peptide_hit)
{
  // if unique protein is present peptide can be assigned
  String accession;
  if (peptide_hit.getProteinAccessions().size() == 1)
  {
    accession = peptide_hit.getProteinAccessions()[0];
  }
  else
  {
    accession = "NA"; // no unique accession
  }
  return accession;
}

String MzTabFile::extractPeptideModifications_(const PeptideHit& peptide_hit)
{
  String mods_string;

  const AASequence& aa_seq = peptide_hit.getSequence();
  bool first = true;

  // check terminal modifications
  if (aa_seq.hasNTerminalModification())
  {
    if (!first)
    {
      mods_string +=  ",";
    }
    else
    {
      first = false;
    }
    String position = "0";
    String unimod_name = aa_seq.getNTerminalModification();
    String unimod_accession =  ModificationsDB::getInstance()->getModification(unimod_name).getUniModAccession();
    mods_string += position  + "-" + unimod_accession;
  }

  if (aa_seq.hasCTerminalModification())
  {
    if (!first)
    {
      mods_string +=  ",";
    }
    else
    {
      first = false;
    }
    String position = String(aa_seq.size() + 1);
    String unimod_name = aa_seq.getCTerminalModification();
    String unimod_accession =  ModificationsDB::getInstance()->getModification(unimod_name).getUniModAccession();
    mods_string += position  + "-" + unimod_accession;
  }

  // check internal modifications
  for (Size i = 0; i != aa_seq.size(); ++i)
  {
    if (aa_seq[i].isModified())
    {
      if (!first)
      {
        mods_string +=  ",";
      }
      else
      {
        first = false;
      }
      String position = String(i + 1);
      // find all modifications with the given name (but different residue/term specifity)
      std::set<const ResidueModification*> modis;
      ModificationsDB::getInstance()->searchModifications(modis, aa_seq[i].getModification(), ResidueModification::ANYWHERE);
      if (!modis.empty())
      {
        // all have the same unimod accession (=record_id) so just take the first one
        set<const ResidueModification*>::const_iterator mit = modis.begin();
        String unimod_accession = (*mit)->getUniModAccession();
        mods_string += position  + "-" + unimod_accession.c_str();
      }
    }
  }

  if (mods_string.length() == 0)
  {
    return "--";
  }

  return mods_string;
}

String MzTabFile::mapSearchEngineToCvParam_(const String& openms_search_engine_name)
{
  String s = openms_search_engine_name;
  s.toUpper();

  if (s == "OMSSA")
  {
    return "[MS,MS:1001475,OMSSA,]";
  }
  else if (s == "MASCOT")
  {
    return "[MS,MS:1001207,MASCOT,]";
  }
  else if (s == "XTANDEM")
  {
    return "[MS,MS:1001476,xtandem,]";
  }
  else if (s == "SEQUEST")
  {
    return "[MS,MS:1001208,Sequest,]";
  }
  else if (s == "COMPNOVO")
  {
    return "[,,CompNovo,]";
  }
  else if (s == "PROTEINPROPHET")
  {
    return "[,,ProteinProphet,]";
  }
  else
  {
    return "NA";
  }
  /*
    TODO:
    additional search engine strings in OpenMS:
    OpenMS/ConsensusID
        InsPecT
        PILIS
        In-silico digestion
        PepNovo
        TurboSEQUEST SEQUEST ???
     */
}

map<String, Size> MzTabFile::extractNumberOfSubSamples_(const map<String, vector<ProteinIdentification> >& map_run_to_proids)
{
  map<String, set<Size> > map_run_to_subsamples_id;

  // for each run...
  for (map<String, vector<ProteinIdentification> >::const_iterator run_it = map_run_to_proids.begin();
       run_it != map_run_to_proids.end(); ++run_it)
  {
    String run = run_it->first;
    const vector<ProteinIdentification>& protein_ids = run_it->second;
    // note: per run there should only exist one protein identification
    for (vector<ProteinIdentification>::const_iterator prot_it = protein_ids.begin();
         prot_it != protein_ids.end(); ++prot_it)
    {
      const ProteinIdentification& protein_id = *prot_it;
      const vector<ProteinHit>& protein_hits = protein_id.getHits();
      // for each ProteinHit...
      for (vector<ProteinHit>::const_iterator pit = protein_hits.begin(); pit != protein_hits.end(); ++pit)
      {
        vector<String> metainfo_keys;
        pit->getKeys(metainfo_keys);
        // find meta values starting with mzTab:protein_abundance_sub
        for (vector<String>::const_iterator s_it = metainfo_keys.begin(); s_it != metainfo_keys.end(); ++s_it)
        {
          //cout << *s_it << endl;
          if (s_it->hasPrefix("mzTab:protein_abundance_sub"))
          {
            String s = *s_it;
            s = s.substitute("mzTab:protein_abundance_sub", "").remove('[').remove(']');
            Size subsample_number = (Size) s.toInt();
            map_run_to_subsamples_id[run].insert(subsample_number);
          }
        }
      }
    }
  }

  // count and return subsample set sizes
  map<String, Size> map_run_to_nsubsamples;
  for (map<String, set<Size> >::const_iterator sub_it = map_run_to_subsamples_id.begin();
       sub_it != map_run_to_subsamples_id.end(); ++sub_it)
  {
    map_run_to_nsubsamples[sub_it->first] = sub_it->second.size();
  }

  return map_run_to_nsubsamples;
}

void MzTabFile::writePeptideHeader_(SVOutStream& output, map<String, Size> n_sub_samples)
{
  output << "PEH" << "sequence" << "accession" << "unit_id" << "unique" << "database"
         << "database_version" << "search_engine" << "search_engine_score"
         << "modifications" << "retention_time" << "charge"
         << "mass_to_charge" << "uri" << "spectra_ref";

  // to generate sufficient number of columns the maximum of sub samples in all runs is used
  Size max_subsamples = 0;
  for (map<String, Size>::const_iterator run_it = n_sub_samples.begin();
       run_it != n_sub_samples.end(); ++run_it)
  {
    if (run_it->second > max_subsamples)
    {
      max_subsamples = run_it->second;
    }
  }

  // print column headers
  for (Size i = 1; i <= max_subsamples; ++i)
  {
    output << String("peptide_abundance_sub[") + String(i) + String("]") <<
              String("peptide_abundance_stdev_sub[") + String(i) + String("]") <<
              String("peptide_abundance_std_error_sub[") + String(i) + String("]");
  }

  output << endl;
}

void MzTabFile::writeProteinHeader_(SVOutStream& output, map<String, Size> n_sub_samples)
{
  output << "PRH" << "accession" << "unit_id" << "description" << "taxid"
         << "species" << "database" << "database_version" << "search_engine"
         << "search_engine_score" << "reliability" << "num_peptides" << "num_peptides_distinct"
         << "num_peptides_unambiguous" << "ambiguity_members" << "modifications" << "uri"
         << "go_terms" << "protein_coverage";

  // to generate sufficient number of columns the maximum of sub samples in all runs is used
  Size max_subsamples = 0;
  for (map<String, Size>::const_iterator run_it = n_sub_samples.begin();
       run_it != n_sub_samples.end(); ++run_it)
  {
    if (run_it->second > max_subsamples)
    {
      max_subsamples = run_it->second;
    }
  }

  // print column headers
  for (Size i = 1; i <= max_subsamples; ++i)
  {
    output << String("protein_abundance_sub[") + String(i) + String("]") <<
              String("protein_abundance_stdev_sub[") + String(i) + String("]") <<
              String("protein_abundance_std_error_sub[") + String(i) + String("]");
  }

  output << endl;
}

// same as distinct but additional constraint of uniquenes (=maps to exactly one Protein)
String MzTabFile::extractNumPeptidesUnambiguous_(String common_identifier, String protein_accession,
                                                 const MapAccPepType& map_run_accesion_to_peptides)
{
  std::pair<String, String> key = make_pair(common_identifier, protein_accession);
  String ret = "0";
  MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
  if (it != map_run_accesion_to_peptides.end())
  {
    const std::vector<PeptideHit>& peptide_hits = it->second;

    // mzTab unambigous peptides are all peptides with different AA sequence OR Modifications
    std::set<String> sequences;
    for (vector<PeptideHit>::const_iterator pet = peptide_hits.begin(); pet != peptide_hits.end(); ++pet)
    {
      // only add sequences of unique peptides
      if (pet->getProteinAccessions().size() == 1)
      {
        sequences.insert(pet->getSequence().toString()); // AASequence with Modifications
      }
    }
    ret = String(sequences.size());
  }
  return ret;
}

void MzTabFile::writeProteinData_(SVOutStream& output,
                                  const ProteinIdentification& prot_id,
                                  Size run_count,
                                  String input_filename,
                                  bool has_coverage,
                                  const MapAccPepType& map_run_accesion_to_peptides,
                                  const map<String, Size>& map_run_to_num_sub
                                  )
{
  // TODO: maybe save these ProteinIdentification run properties in meta data
  // it->getScoreType()
  // it->isHigherScoreBetter())
  // it->getDateTime().toString(Qt::ISODate).toStdString()
  // it->getSearchEngineVersion();

  // search parameters
  const ProteinIdentification::SearchParameters& sp = prot_id.getSearchParameters();
  // TODO: maybe save these SearchParameters properties in a user param
  // String charges; ///< The allowed charges for the search
  // PeakMassType mass_type; ///< Mass type of the peaks
  // std::vector<String> fixed_modifications; ///< Used fixed modifications
  // std::vector<String> variable_modifications; ///< Allowed variable modifications
  // ProteinIdentification::NamesOfDigestionEnzyme[sp.enzyme]
  // UInt missed_cleavages; ///< The number of allowed missed cleavages
  // double peak_mass_tolerance; ///< Mass tolerance of fragment ions (Dalton)
  // double precursor_tolerance; ///< Mass tolerance of precursor ions (Dalton)

  // in OpenMS global to a ProteinIdentification
  String UNIT_ID_String = "OpenMS_" + String(run_count);
  String database_String = (sp.db != "" ? sp.db : "--");
  database_String = "file://" +  database_String;
  String database_version_String = (sp.db_version != "" ? sp.db_version : "--");
  String species_String =  (sp.taxonomy == "0" || sp.taxonomy == "" ? "--" : sp.taxonomy);
  String search_engine_cvParams = mapSearchEngineToCvParam_(prot_id.getSearchEngine());
  String openms_search_engine_name = prot_id.getSearchEngine();
  //
  for (vector<ProteinHit>::const_iterator protein_hit_it = prot_id.getHits().begin();
       protein_hit_it != prot_id.getHits().end(); ++protein_hit_it)
  {
    String accession = protein_hit_it->getAccession();
    String unit_id = UNIT_ID_String; // run specific in OpenMS
    String description = "--"; // TODO: support description in protein hit
    String taxid = "--"; // TODO: mapping to NCBI taxid needed
    String species = species_String; // run specific in OpenMS
    String database = database_String; // run specific in OpenMS
    String database_version = database_version_String; // run specific in OpenMS
    String search_engine = search_engine_cvParams;
    String search_engine_score = mapSearchEngineScoreToCvParam_(openms_search_engine_name,
                                                                protein_hit_it->getScore(),
                                                                prot_id.getScoreType());

    String reliability = "--";
    String num_peptides;
    String num_peptides_distinct = extractNumPeptidesDistinct_(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
    String num_peptides_unambiguous = extractNumPeptidesUnambiguous_(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
    String ambiguity_members = "NA"; //TODO
    String modifications = "NA"; // TODO
    String uri = "file://" + input_filename;
    String go_terms = "--";
    String protein_coverage;

    if (has_coverage)
    {
      protein_coverage = String(protein_hit_it->getCoverage() / 100.0);
    }
    else
    {
      protein_coverage = "NA";
    }

    if (protein_hit_it->metaValueExists("num_peptides"))
    {
      num_peptides = protein_hit_it->getMetaValue("num_peptides");
    }
    else
    {
      num_peptides = extractNumPeptides_(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
    }

    output << "PRT" << accession << unit_id << description << taxid
           << species << database << database_version << search_engine
           << search_engine_score << reliability << num_peptides << num_peptides_distinct
           << num_peptides_unambiguous << ambiguity_members << modifications << uri
           << go_terms << protein_coverage;

    // get number of sub samples for this run
    map<String, Size>::const_iterator sub_it = map_run_to_num_sub.find(prot_id.getIdentifier());
    Size n_subsamples = 0;
    if (sub_it != map_run_to_num_sub.end())
    {
      n_subsamples = sub_it->second;
    }

    for (Size n = 1; n <= n_subsamples; ++n)
    {
      {
        String key = "mzTab:protein_abundance_sub[" + String(n) + "]";
        String abundancy_value = "--";
        if (protein_hit_it->metaValueExists(key))
        {
          abundancy_value = protein_hit_it->getMetaValue(key);
        }
        output << abundancy_value;
      }
      {
        String key = "mzTab:protein_abundance_stdev_sub[" + String(n) + "]";
        String abundancy_value = "--";
        if (protein_hit_it->metaValueExists(key))
        {
          abundancy_value = protein_hit_it->getMetaValue(key);
        }
        output << abundancy_value;
      }
      {
        String key = "mzTab:protein_abundance_std_error_sub[" + String(n) + "]";
        String abundancy_value = "--";
        if (protein_hit_it->metaValueExists(key))
        {
          abundancy_value = protein_hit_it->getMetaValue(key);
        }
        output << abundancy_value;
      }
    }
    output << endl;
  }
}

String MzTabFile::extractNumPeptidesDistinct_(String common_identifier, String protein_accession,
                                              const MapAccPepType& map_run_accesion_to_peptides)
{
  std::pair<String, String> key = make_pair(common_identifier, protein_accession);
  String ret = "0";
  MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
  if (it != map_run_accesion_to_peptides.end())
  {
    const vector<PeptideHit>& peptide_hits = it->second;

    // mzTab unambigous peptides are all peptides with different AA sequence OR Modifications
    std::set<String> sequences;
    for (vector<PeptideHit>::const_iterator pet = peptide_hits.begin(); pet != peptide_hits.end(); ++pet)
    {
      sequences.insert(pet->getSequence().toString()); // AASequence including Modifications
    }

    ret = String(sequences.size());
  }

  return ret;
}

String MzTabFile::extractNumPeptides_(const String& common_identifier, const String& protein_accession,
                                      const MapAccPepType& map_run_accesion_to_peptides)
{
  pair<String, String> key = make_pair(common_identifier, protein_accession);
  String ret = "0";
  MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
  if (it != map_run_accesion_to_peptides.end())
  {
    const vector<PeptideHit>& peptide_hits = it->second;
    ret = String(peptide_hits.size());
  }
  return ret;
}

String MzTabFile::mapSearchEngineScoreToCvParam_(const String& openms_search_engine_name, double score, String score_type)
{
  String s;

  if (score_type.hasSubstring("Consensus"))
  {
    s = "[,,Consensus:score,";
  }
  else if (score_type == "q-value")
  {
    s = "[MS,MS:1001364,pep:global FDR,";
  }
  else if (score_type == "FDR")
  {
    s = "[MS,MS:1001364,pep:global FDR,";
  }
  else if (score_type == "Posterior Error Probability")
  {
    s = "[,,PEP,";
  }
  else if (score_type == "PhosphoScore")
  {
    s = "[,,PhosphoScore,";
  }
  else if (openms_search_engine_name == "OMSSA")
  {
    s = "[MS,MS:1001328,OMSSA:evalue,";
  }
  else if (openms_search_engine_name == "Mascot")
  {
    s = "[MS,MS:1001171,MASCOT:score,";
  }
  else if (openms_search_engine_name == "XTandem")
  {
    s = "[MS,MS:1001330,X!Tandem:expect,";
  }
  else if (openms_search_engine_name == "SEQUEST")
  {
    s = "[MS,MS:1001155,Sequest:xcorr,";
  }
  else if (openms_search_engine_name == "CompNovo")
  {
    s = "[,,CompNovo,";
  }
  else if (score_type == "ProteinProphet probability")
  {
    s = "[,,ProteinProphet,";
  }
  else
  {
    return "NA";
  }

  s += String::number(score, 8) + "]";
  return s;
  /*
    TODO:
    additional search engine strings in OpenMS:
    OpenMS/ConsensusID
        InsPecT
        PILIS und PILIS-E-value
        In-silico digestion
        PepNovo
        TurboSEQUEST SEQUEST ???
        InterProphet probability
        ProteinProphet probability
     */
}

void MzTabFile::sortPSM_(vector<PeptideIdentification>::iterator begin, vector<PeptideIdentification>::iterator end)
{
  for (vector<PeptideIdentification>::iterator pep_id_it = begin; pep_id_it != end; ++pep_id_it)
  {
    pep_id_it->assignRanks();
  }
}

void MzTabFile::keepFirstPSM_(vector<PeptideIdentification>::iterator begin, vector<PeptideIdentification>::iterator end)
{
  IDFilter id_filter;
  for (vector<PeptideIdentification>::iterator pep_id_it = begin; pep_id_it != end; ++pep_id_it)
  {
    PeptideIdentification new_pep_id;
    id_filter.filterIdentificationsByBestHits(*pep_id_it, new_pep_id, false);
    pep_id_it->setHits(new_pep_id.getHits());
  }
}

}
