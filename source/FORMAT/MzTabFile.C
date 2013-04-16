// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
#include <fstream>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

  MzTabFile::MzTabFile()
  {

  }

  MzTabFile::~MzTabFile()
  {

  }

  void MzTabFile::generateMzTabMetaDataSection_(const MzTabMetaData& map, StringList& sl) const
  {
    for (MzTabMetaData::const_iterator it = map.begin(); it != map.end(); ++it)
    {
      const String& unit_id = it->first;
      const MzTabUnitIdMetaData& md = it->second;

      //	{UNIT_ID}-title
      if (!md.title.empty())
      {
        String s = String("MTD\t") + unit_id + "-title\t" + md.title;
        sl << s;
      }

      // {UNIT_ID}-description
      if (!md.description.empty())
      {
        String s = String("MTD\t") + unit_id + "-description\t" + md.description;
        sl << s;
      }

      // {UNIT_ID}-sample_processing[1-n]
      for (Size i = 0; i != md.sample_processing.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-sample_processing[" + String(i + 1) + "]\t" + md.sample_processing[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-instrument[1-n]-name
      for (Size i = 0; i != md.instrument_name.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-instrument[" + String(i + 1) + "]-name\t" + md.instrument_name[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-instrument[1-n]-source
      for (Size i = 0; i != md.instrument_source.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-instrument[" + String(i + 1) + "]-source\t" + md.instrument_source[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-instrument[1-n]-analyzer
      for (Size i = 0; i != md.instrument_analyzer.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-instrument[" + String(i + 1) + "]-analyzer\t" + md.instrument_analyzer[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-instrument[1-n]-detector
      for (Size i = 0; i != md.instrument_detector.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-instrument[" + String(i + 1) + "]-detector\t" + md.instrument_detector[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-software[1-n]
      for (Size i = 0; i != md.software.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-software[" + String(i + 1) + "]\t" + md.software[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-software[1-n]-setting. Note that lines for a single software (=index) MAY occur multiple times.
      for (Size i = 0; i != md.software_setting.size(); ++i)
      {
        for (Size j = 0; j != md.software_setting[i].size(); ++j)
        {
          String s = "MTD\t" + unit_id + "-software[" + String(i + 1) + "]-setting\t" + md.software_setting[i][j];
          sl << s;
        }
      }

      // {UNIT_ID}-false_discovery_rate
      for (Size i = 0; i != md.false_discovery_rate.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-false_discovery_rate\t" + md.false_discovery_rate[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-publication
      for (Size i = 0; i != md.publication.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-publication\t" + md.publication[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-contact[1-n]-name
      for (Size i = 0; i != md.contact_name.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-contact[" + String(i + 1) + "]-name\t" + md.contact_name[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-contact[1-n]-affiliation
      for (Size i = 0; i != md.contact_affiliation.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-contact[" + String(i + 1) + "]-affiliation\t" + md.contact_name[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-contact[1-n]-email
      for (Size i = 0; i != md.contact_email.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-contact[" + String(i + 1) + "]-email\t" + md.contact_name[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-uri
      for (Size i = 0; i != md.uri.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-uri\t" + md.uri[i];
        sl << s;
      }

      // {UNIT_ID}-mod
      if (!md.mod.isNull())
      {
        String s = "MTD\t" + unit_id + "-mod\t" + md.mod.toCellString();
        sl << s;
      }

      // {UNIT_ID}-quantification_method
      if (!md.quantification_method.isNull())
      {
        String s = "MTD\t" + unit_id + "-quantification_method\t" + md.quantification_method.toCellString();
        sl << s;
      }

      // {UNIT_ID}-protein-quantification_unit
      if (!md.protein_quantification_unit.isNull())
      {
        String s = "MTD\t" + unit_id + "-protein-quantification_unit\t" + md.protein_quantification_unit.toCellString();
        sl << s;
      }

      // {UNIT_ID}-peptide-quantification_unit
      if (!md.peptide_quantification_unit.isNull())
      {
        String s = "MTD\t" + unit_id + "-peptide_quantification-unit\t" + md.peptide_quantification_unit.toCellString();
        sl << s;
      }

      // {UNIT_ID}-small_molecule-quantification_unit
      if (!md.small_molecule_quantification_unit.isNull())
      {
        String s = "MTD\t" + unit_id + "-small_molecule-quantification_unit\t" + md.small_molecule_quantification_unit.toCellString();
        sl << s;
      }

      // {UNIT_ID}-ms_file[1-n]-format
      for (Size i = 0; i != md.ms_file_format.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-ms_file[" + String(i + 1) + "]-format\t" + md.ms_file_format[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-ms_file[1-n]-location
      for (Size i = 0; i != md.ms_file_location.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-ms_file[" + String(i + 1) + "]-location\t" + md.ms_file_location[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-ms_file[1-n]-id_format
      for (Size i = 0; i != md.ms_file_id_format.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-ms_file[" + String(i + 1) + "]-id_format\t" + md.ms_file_id_format[i].toCellString();
        sl << s;
      }

      // {UNIT_ID}-custom
      for (Size i = 0; i != md.custom.size(); ++i)
      {
        String s = "MTD\t" + unit_id + "-custom\t" + md.custom[i].toCellString();
        sl << s;
      }

      for (Size i = 0; i != md.sub_id_data.size(); ++i)
      {
        const MzTabSubIdMetaData& submd = md.sub_id_data[i];

        String sub_id = submd.sub_id;

        // format is UNIT_ID-SUB_ID... if SUB_ID is given
        if (!sub_id.empty())
        {
          sub_id = String("-") + sub_id;
        }

        // {UNIT_ID}(-{SUB_ID})-species[1-n]
        for (Size j = 0; j != submd.species.size(); ++j)
        {
          String s = "MTD\t" + unit_id + sub_id + "-species[" + String(j + 1) + "]\t" + submd.species[j].toCellString();
          sl << s;
        }

        // {UNIT_ID}(-{SUB_ID})-tissue[1-n]
        for (Size j = 0; j != submd.tissue.size(); ++j)
        {
          String s = "MTD\t" + unit_id + sub_id + "-tissue[" + String(j + 1) + "]\t" + submd.tissue[j].toCellString();
          sl << s;
        }

        // {UNIT_ID}(-{SUB_ID})-cell_type[1-n]
        for (Size j = 0; j != submd.cell_type.size(); ++j)
        {
          String s = "MTD\t" + unit_id + sub_id + "-cell_type[" + String(j + 1) + "]\t" + submd.cell_type[j].toCellString();
          sl << s;
        }

        // {UNIT_ID}(-{SUB_ID})-disease[1-n]
        for (Size j = 0; j != submd.disease.size(); ++j)
        {
          String s = "MTD\t" + unit_id + sub_id + "-disease[" + String(j + 1) + "]\t" + submd.disease[j].toCellString();
          sl << s;
        }

        // {UNIT_ID}(-{SUB_ID})-description[1-n]
        for (Size j = 0; j != submd.description.size(); ++j)
        {
          String s = String("MTD\t") + unit_id + sub_id + "-description\t" + submd.description[j];
          sl << s;
        }

        // {UNIT_ID}-{SUB_ID}-quantification_reagent
        for (Size j = 0; j != submd.quantification_reagent.size(); ++j)
        {
          String s = String("MTD\t") + unit_id + sub_id + "-quantification_reagent\t" + submd.quantification_reagent[j].toCellString();
          sl << s;
        }

        // {UNIT_ID}-{SUB_ID}-custom
        for (Size j = 0; j != submd.custom.size(); ++j)
        {
          String s = String("MTD\t") + unit_id + sub_id + "-custom" + submd.custom[j].toCellString();
          sl << s;
        }
      }

      // {UNIT_ID}-colunit-protein
      for (Size i = 0; i != md.colunit_protein.size(); ++i)
      {
        String s = String("MTD\t") + unit_id + "-colunit-protein" + md.colunit_protein[i];
        sl << s;
      }

      // {UNIT_ID}-colunit-peptide
      for (Size i = 0; i != md.colunit_peptide.size(); ++i)
      {
        String s = String("MTD\t") + unit_id + "-colunit-peptide" + md.colunit_peptide[i];
        sl << s;
      }

      // {UNIT_ID}-colunit-small_molecule
      for (Size i = 0; i != md.colunit_small_molecule.size(); ++i)
      {
        String s = String("MTD\t") + unit_id + "-colunit-small_molecule" + md.colunit_small_molecule[i];
        sl << s;
      }
    }
  }

  void MzTabFile::generateProteinHeader_(Int n_subsamples, const vector<String>& optional_protein_columns, StringList& sl) const
  {
    StringList header;
    header << "PRH"
        <<  "accession" << "unit_id" << "description" << "taxid" << "species"
        << "database" << "database_version" << "search_engine"
        << "search_engine_score" << "reliability" << "num_peptides"
        << "num_peptides_distinct" << "num_peptides_unambiguous"
        << "ambiguity_members" << "modifications" << "uri" << "go_terms" <<"protein_coverage";

    for (Int i = 1; i <= n_subsamples; ++i)
    {
      header << String("protein_abundance_sub[") + String(i) + String("]")
             << String("protein_abundance_stdev_sub[") + String(i) + String("]")
             << String("protein_abundance_std_error_sub[") + String(i) + String("]");
    }

    std::copy(optional_protein_columns.begin(), optional_protein_columns.end(), std::back_inserter(header));

    sl << (header.concatenate("\t"));
  }

  String MzTabFile::generateMzTabProteinSectionRow_(const MzTabProteinSectionRow& row, const String& unit_id) const
  {
    StringList s;
    s << "PRT"
      << row.accession.toCellString() << unit_id << row.description.toCellString()
      << row.taxid.toCellString() << row.species.toCellString() << row.database.toCellString() << row.database_version.toCellString()
      << row.search_engine.toCellString() << row.search_engine_score.toCellString() << row.reliability.toCellString()
      << row.num_peptides.toCellString() << row.num_peptides_distinct.toCellString() << row.num_peptides_unambiguous.toCellString()
      << row.modifications.toCellString() << row.uri.toCellString() << row.go_terms.toCellString() << row.protein_coverage.toCellString();

    // quantification columns
    Size nsub = row.protein_abundance_sub.size();
    if (nsub != row.protein_abundance_stdev_sub.size() || nsub != row.protein_abundance_std_error_sub.size())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Different number of columns for sub quantification.'"));
    }

    for (Size i = 0; i != nsub; ++i)
    {
      s << row.protein_abundance_sub[i]
        << row.protein_abundance_stdev_sub[i]
        << row.protein_abundance_std_error_sub[i];
    }

    // print optional columns
    for (Size i = 0; i <= row.opt_.size(); ++i)
    {
      s << row.opt_[i].second.toCellString();
    }

    return s.concatenate("\t");
  }

  void MzTabFile::generateMzTabProteinSection_(const MzTabProteinSectionData& map, StringList& sl) const
  {
    for (MzTabProteinSectionData::const_iterator it = map.begin(); it != map.end(); ++it)
    {
      String unit_id = it->first;
      for (MzTabProteinSectionRows::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
      {
        sl << generateMzTabProteinSectionRow_(*jt, unit_id);
      }
    }
  }

  void MzTabFile::generateMzTabPeptideSection_(const MzTabPeptideSectionData& map, StringList& sl) const
  {
    for (MzTabPeptideSectionData::const_iterator it = map.begin(); it != map.end(); ++it)
    {
      String unit_id = it->first;
      for (MzTabPeptideSectionRows::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
      {
        sl << generateMzTabPeptideSectionRow_(*jt, unit_id);
      }
    }
  }

  void MzTabFile::generateMzTabSmallMoleculeSection_(const MzTabSmallMoleculeSectionData & map, StringList& sl) const
  {
    for (MzTabSmallMoleculeSectionData::const_iterator it = map.begin(); it != map.end(); ++it)
    {
      String unit_id = it->first;
      for (MzTabSmallMoleculeSectionRows::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
      {
        sl << generateMzTabSmallMoleculeSectionRow_(*jt, unit_id);
      }
    }
  }

  String MzTabFile::generateMzTabPeptideHeader_(Int n_subsamples, const vector<String>& optional_protein_columns) const
  {
    StringList header;
    header << "PEH"
        << "sequence" <<  "accession" << "unit_id" << "unique"
        << "database" << "database_version" << "search_engine"
        << "search_engine_score" << "reliability" << "modifications"
        << "retention_time" << "charge"
        << "mass_to_charge" << "uri" << "spectra_ref";

    for (Int i = 1; i <= n_subsamples; ++i)
    {
      header << String("peptide_abundance_sub[") + String(i) + String("]")
          << String("peptide_abundance_stdev_sub[") + String(i) + String("]")
          << String("peptide_abundance_std_error_sub[") + String(i) + String("]");
    }

    std::copy(optional_protein_columns.begin(), optional_protein_columns.end(), std::back_inserter(header));

    return header.concatenate("\t");
  }

  String MzTabFile::generateMzTabPeptideSectionRow_(const MzTabPeptideSectionRow& row, const String& unit_id) const
  {
    StringList s;
    s << "PEP"
      << row.sequence.toCellString() << row.accession.toCellString() << unit_id
      << row.unique.toCellString() << row.database.toCellString() << row.database_version.toCellString()
      << row.search_engine.toCellString() << row.search_engine_score.toCellString() << row.reliability.toCellString()
      << row.modifications.toCellString() << row.retention_time.toCellString() << row.charge.toCellString() << row.mass_to_charge.toCellString()
      << row.uri.toCellString() << row.spectra_ref.toCellString();

    // quantification columns
    Size nsub = row.peptide_abundance_sub.size();
    if (nsub != row.peptide_abundance_stdev_sub.size() || nsub != row.peptide_abundance_std_error_sub.size())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Different number of columns for sub quantification.'"));
    }

    for (Size i = 0; i != nsub; ++i)
    {
      s << row.peptide_abundance_sub[i]
        << row.peptide_abundance_stdev_sub[i]
        << row.peptide_abundance_std_error_sub[i];
    }

    // print optional columns
    for (Size i = 0; i <= row.opt_.size(); ++i)
    {
      s << row.opt_[i].second.toCellString();
    }

    return s.concatenate("\t");
  }

  String MzTabFile::generateMzTabSmallMoleculeHeader_(Int n_subsamples, const vector<String>& optional_smallmolecule_columns) const
  {
    StringList header;
    header << "SMH"
        << "identifier" << "unit_id" << "chemical_formula"
        << "smiles" << "inchi_key" << "description"
        << "mass_to_charge" << "charge" << "retention_time"
        << "taxid" << "species" << "database" << "database_version"
        << "reliability" << "uri" << "spectra_ref"
        << "search_engine" << "search_engine_score" << "modifications";

    for (Int i = 1; i <= n_subsamples; ++i)
    {
      header << String("smallmolecule_abundance_sub[") + String(i) + String("]")
          << String("smallmolecule_abundance_stdev_sub[") + String(i) + String("]")
          << String("smallmolecule_std_error_sub[") + String(i) + String("]");
    }

    // copy optional column names to header
    std::copy(optional_smallmolecule_columns.begin(), optional_smallmolecule_columns.end(), std::back_inserter(header));

    return header.concatenate("\t");
  }

  String MzTabFile::generateMzTabSmallMoleculeSectionRow_(const MzTabSmallMoleculeSectionRow& row, const String& unit_id) const
  {
    StringList s;
    s << "SML"
      << row.identifier.toCellString() << unit_id << row.chemical_formula.toCellString()
      << row.smiles.toCellString() << row.inchi_key.toCellString() << row.description.toCellString()
      << row.mass_to_charge.toCellString() << row.charge.toCellString() << row.retention_time.toCellString()
      << row.taxid.toCellString() << row.species.toCellString() << row.database.toCellString()
      << row.database_version.toCellString() << row.reliability.toCellString() << row.uri.toCellString()
      << row.spectra_ref.toCellString() << row.search_engine.toCellString() << row.search_engine_score.toCellString()
      << row.modifications.toCellString();

    // quantification columns
    Size nsub = row.smallmolecule_abundance_sub.size();
    if (nsub != row.smallmolecule_abundance_stdev_sub.size() || nsub != row.smallmolecule_abundance_std_error_sub.size())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Different number of columns for sub quantification.'"));
    }

    for (Size i = 0; i != nsub; ++i)
    {
      s << row.smallmolecule_abundance_sub[i].toCellString()
        << row.smallmolecule_abundance_stdev_sub[i].toCellString()
        << row.smallmolecule_abundance_std_error_sub[i].toCellString();
    }

    // print optional columns
    for (Size i = 0; i <= row.opt_.size(); ++i)
    {
      s << row.opt_[i].second.toCellString();
    }

    return s.concatenate("\t");
  }

  void MzTabFile::store(const String & filename, const MzTab& mz_tab) const
  {
    TextFile out;
    generateMzTabMetaDataSection_(mz_tab.getMetaData(), out);
    generateMzTabProteinSection_(mz_tab.getProteinSectionData(), out);
    generateMzTabPeptideSection_(mz_tab.getPeptideSectionData(), out);
    generateMzTabSmallMoleculeSection_(mz_tab.getSmallMoleculeSectionData(), out);
    out.store(filename);
  }

  void MzTabFile::store(const String & filename, const std::vector<ProteinIdentification> & protein_ids, const std::vector<PeptideIdentification> & peptide_ids, String in, String document_id) const
  {
    // copy data as we are going to apply some filters
    vector<ProteinIdentification> prot_ids = protein_ids;
    vector<PeptideIdentification> pep_ids = peptide_ids;

    // create tab separated output stream for MzTab file
    ofstream txt_out(filename.c_str());
    SVOutStream output(txt_out, "\t", "_", String::NONE);

    Size num_runs = 0;

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
    catch (Exception::MissingInformation & e)
    {
      cout << e.what() << "\n";
      has_coverage = false;
    }

    // partition into runs (as this could be a merged idXML)
    map<String, vector<PeptideIdentification> > map_run_to_pepids;
    map<String, vector<ProteinIdentification> > map_run_to_proids;
    partitionIntoRuns(pep_ids, prot_ids, map_run_to_pepids, map_run_to_proids);
    num_runs = map_run_to_pepids.size();
    cout << "idXML contains: " << num_runs << " runs." << endl;

    MapAccPepType map_run_accesion_to_peptides;
    createProteinToPeptideLinks(map_run_to_pepids, map_run_accesion_to_peptides);

    // every ProteinIdentification corresponds to a search engine run and contains protein hits

    // write meta data for each run
    bool meta_info_printed = false;
    Size run_count = 0;
    map<String, vector<ProteinIdentification> >::const_iterator mprot_it = map_run_to_proids.begin();
    map<String, vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();

    for (; mprot_it != map_run_to_proids.end(); ++mprot_it) // iterate over runs
    {
      // extract ProteinIdentifications of this run (should be only 1)
      const vector<ProteinIdentification> & prot_ids = mprot_it->second;
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
      const vector<ProteinIdentification> & prot_ids = mprot_it->second;
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
      const vector<ProteinIdentification> & prot_ids = mprot_it->second;
      // extract PeptideIdentifications of this run
      const vector<PeptideIdentification> & pep_ids = mpep_it->second;
      // iterate over runs of peptide identifications
      for (vector<PeptideIdentification>::const_iterator pep_id_it = pep_ids.begin();
           pep_id_it != pep_ids.end(); ++pep_id_it)
      {
        // TODO: check if bad design of Protein/PeptideIdentification as search engine parameters are stored in prot.
        String openms_search_engine_name = prot_ids[0].getSearchEngine();
        String search_engine_cvParams = mapSearchEngineToCvParam_(openms_search_engine_name);

        const ProteinIdentification::SearchParameters & sp = prot_ids[0].getSearchParameters();
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
  void MzTabFile::partitionIntoRuns(const vector<PeptideIdentification> & pep_ids,
                                    const vector<ProteinIdentification> & pro_ids,
                                    map<String, vector<PeptideIdentification> > & map_run_to_pepids,
                                    map<String, vector<ProteinIdentification> > & map_run_to_proids
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
      map<String, vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();
      map<String, vector<ProteinIdentification> >::const_iterator mprot_it = map_run_to_proids.begin();
      for (; mpep_it != map_run_to_pepids.end(); ++mpep_it, ++mprot_it)
      {
        assert(mpep_it->first == mprot_it->first);
      }
    }
  }

  void MzTabFile::createProteinToPeptideLinks(const map<String, vector<PeptideIdentification> > & map_run_to_pepids, MapAccPepType & map_run_accession_to_pephits)
  {
    // create links for each run
    map<String, vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();

    for (; mpep_it != map_run_to_pepids.end(); ++mpep_it)
    {
      const String & run = mpep_it->first;
      const vector<PeptideIdentification> & pids = mpep_it->second;
      for (Size i = 0; i != pids.size(); ++i)
      {
        const std::vector<PeptideHit> & phits = pids[i].getHits();
        if (phits.size() == 1)
        {
          const std::vector<String> & accessions = phits[0].getProteinAccessions();
          for (Size k = 0; k != accessions.size(); ++k)
          {
            const String & accession = accessions[k];
            std::pair<String, String> key = make_pair(run, accession);
            map_run_accession_to_pephits[key].push_back(phits[0]);
          }
        }
      }
    }
  }

/// Extracts, if possible a unique protein accession for a peptide hit in mzTab format. Otherwise NA is returned
  String MzTabFile::extractProteinAccession_(const PeptideHit & peptide_hit)
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

  String MzTabFile::extractPeptideModifications_(const PeptideHit & peptide_hit)
  {
    String mods_string;

    const AASequence & aa_seq = peptide_hit.getSequence();
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
        std::set<const ResidueModification *> modis;
        ModificationsDB::getInstance()->searchModifications(modis, aa_seq[i].getModification(), ResidueModification::ANYWHERE);
        if (!modis.empty())
        {
          // all have the same unimod accession (=record_id) so just take the first one
          set<const ResidueModification *>::const_iterator mit = modis.begin();
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

  String MzTabFile::mapSearchEngineToCvParam_(const String & openms_search_engine_name)
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

  map<String, Size> MzTabFile::extractNumberOfSubSamples_(const map<String, vector<ProteinIdentification> > & map_run_to_proids)
  {
    map<String, set<Size> > map_run_to_subsamples_id;

    // for each run...
    for (map<String, vector<ProteinIdentification> >::const_iterator run_it = map_run_to_proids.begin();
         run_it != map_run_to_proids.end(); ++run_it)
    {
      String run = run_it->first;
      const vector<ProteinIdentification> & protein_ids = run_it->second;
      // note: per run there should only exist one protein identification
      for (vector<ProteinIdentification>::const_iterator prot_it = protein_ids.begin();
           prot_it != protein_ids.end(); ++prot_it)
      {
        const ProteinIdentification & protein_id = *prot_it;
        const vector<ProteinHit> & protein_hits = protein_id.getHits();
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

  void MzTabFile::writePeptideHeader_(SVOutStream & output, map<String, Size> n_sub_samples)
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

  void MzTabFile::writeProteinHeader_(SVOutStream & output, map<String, Size> n_sub_samples)
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
  String MzTabFile::extractNumPeptidesUnambiguous(String common_identifier, String protein_accession,
                                                  const MapAccPepType & map_run_accesion_to_peptides)
  {
    std::pair<String, String> key = make_pair(common_identifier, protein_accession);
    String ret = "0";
    MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
    if (it != map_run_accesion_to_peptides.end())
    {
      const std::vector<PeptideHit> & peptide_hits = it->second;

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

  void MzTabFile::writeProteinData_(SVOutStream & output,
                                    const ProteinIdentification & prot_id,
                                    Size run_count,
                                    String input_filename,
                                    bool has_coverage,
                                    const MapAccPepType & map_run_accesion_to_peptides,
                                    const map<String, Size> & map_run_to_num_sub
                                    )
  {
    // TODO: maybe save these ProteinIdentification run properties in meta data
    // it->getScoreType()
    // it->isHigherScoreBetter())
    // it->getDateTime().toString(Qt::ISODate).toStdString()
    // it->getSearchEngineVersion();

    // search parameters
    const ProteinIdentification::SearchParameters & sp = prot_id.getSearchParameters();
    // TODO: maybe save these SearchParameters properties in a user param
    // String charges; ///< The allowed charges for the search
    // PeakMassType mass_type; ///< Mass type of the peaks
    // std::vector<String> fixed_modifications; ///< Used fixed modifications
    // std::vector<String> variable_modifications; ///< Allowed variable modifications
    // ProteinIdentification::NamesOfDigestionEnzyme[sp.enzyme]
    // UInt missed_cleavages; ///< The number of allowed missed cleavages
    // DoubleReal peak_mass_tolerance; ///< Mass tolerance of fragment ions (Dalton)
    // DoubleReal precursor_tolerance; ///< Mass tolerance of precursor ions (Dalton)

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
      String num_peptides_distinct = extractNumPeptidesDistinct(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
      String num_peptides_unambiguous = extractNumPeptidesUnambiguous(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
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
        num_peptides = extractNumPeptides(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
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

  String MzTabFile::extractNumPeptidesDistinct(String common_identifier, String protein_accession,
                                               const MapAccPepType & map_run_accesion_to_peptides)
  {
    std::pair<String, String> key = make_pair(common_identifier, protein_accession);
    String ret = "0";
    MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
    if (it != map_run_accesion_to_peptides.end())
    {
      const vector<PeptideHit> & peptide_hits = it->second;

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

  String MzTabFile::extractNumPeptides(const String & common_identifier, const String & protein_accession,
                                       const MapAccPepType & map_run_accesion_to_peptides)
  {
    pair<String, String> key = make_pair(common_identifier, protein_accession);
    String ret = "0";
    MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
    if (it != map_run_accesion_to_peptides.end())
    {
      const vector<PeptideHit> & peptide_hits = it->second;
      ret = String(peptide_hits.size());
    }
    return ret;
  }

  String MzTabFile::mapSearchEngineScoreToCvParam_(const String & openms_search_engine_name, DoubleReal score, String score_type)
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
