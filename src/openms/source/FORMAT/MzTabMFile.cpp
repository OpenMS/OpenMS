// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabMFile.h>

#include <OpenMS/FORMAT/TextFile.h>

#include <boost/regex.hpp>

using namespace std;

namespace OpenMS
{

  MzTabMFile::MzTabMFile(){}

  MzTabMFile::~MzTabMFile(){}

  void MzTabMFile::generateMzTabMMetaDataSection_(const MzTabMMetaData& md, StringList& sl) const
  {
    std::cout << "run metadata generation" << std::endl;
    sl.push_back(String("MTD\tmzTab-version\t") + md.mz_tab_version.toCellString()); // mandatory
    sl.push_back(String("MTD\tmzTab-ID\t") + md.mz_tab_id.toCellString()); // mandatory

    if (!md.title.isNull())
    {
      String s = String("MTD\ttitle\t") + md.title.toCellString();
      sl.push_back(s);
    }

    if(!md.description.isNull())
    {
      String s = String("MTD\tdescription\t") + md.description.toCellString();
      sl.push_back(s);
    }

    for (map<Size, MzTabParameterList>::const_iterator it = md.sample_processing.begin(); it != md.sample_processing.end(); ++it)
    {
      String s = "MTD\tsample_processing[" + String(it->first) + "]\t" + it->second.toCellString();
      sl.push_back(s);
    }

    for (const auto& sp : md.sample_processing)
    {
      String s = "MTD\tsample_processing[" + String(sp.first) + "]\t" + sp.second.toCellString();
      sl.push_back(s);
    }

    for (const auto& inst : md.instrument)
    {
      const MzTabInstrumentMetaData& imd = inst.second;

      if (!imd.name.isNull())
      {
        String s = "MTD\tinstrument[" + String(inst.first) + "]-name\t" +  imd.name.toCellString();
        sl.push_back(s);
      }

      if (!imd.source.isNull())
      {
        String s = "MTD\tinstrument[" + String(inst.first) + "]-source\t" + imd.source.toCellString();
        sl.push_back(s);
      }

      for (const auto& mit : imd.analyzer)
      {
        if (!mit.second.isNull())
        {
          String s = "MTD\tinstrument[" + String(inst.first) + "]-analyzer[" + String(mit.first) + "]\t" + mit.second.toCellString();
          sl.push_back(s);
        }
      }

      if (!imd.detector.isNull())
      {
        String s = "MTD\tinstrument[" + String(inst.first) + "]-detector\t" + imd.detector.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto & sw : md.software)
    {
      MzTabSoftwareMetaData msmd = sw.second;

      String s = "MTD\tsoftware[" + String(sw.first) + "]\t" + msmd.software.toCellString(); // mandatory
      sl.push_back(s);

      for (const auto& setting : msmd.setting)
      {
        String s = "MTD\tsoftware[" + String(sw.first) + "]-setting[" + String(setting.first) + String("]\t") + setting.second.toCellString();
        sl.push_back(s);
      }
    }

    for (auto const& pub : md.publication)
    {
      String s = "MTD\tpublication[" + String(pub.first) + "]\t" + pub.second.toCellString();
      sl.push_back(s);
    }

    for (const auto& contact : md.contact)
    {
      const MzTabContactMetaData& md = contact.second;
      if (!md.name.isNull())
      {
        String s = "MTD\tcontact[" + String(contact.first) + "]-name\t" + md.name.toCellString();
        sl.push_back(s);
      }

      if (!md.affiliation.isNull())
      {
        String s = "MTD\tcontact[" + String(contact.first) + "]-affiliation\t" + md.affiliation.toCellString();
        sl.push_back(s);
      }

      if (!md.email.isNull())
      {
        String s = "MTD\tcontact[" + String(contact.first) + "]-email\t" + md.email.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto& uri : md.uri)
    {
      String s = "MTD\turi[" + String(uri.first) + "]\t" + uri.second.toCellString();
      sl.push_back(s);
    }

    for (const auto& ext_study : md.external_study_uri)
    {
      String s = "MTD\texternal_study_uri[" + String(ext_study.first) + "]\t" + ext_study.second.toCellString();
      sl.push_back(s);
    }

    String s = String("MTD\tquantification_method\t") + md.quantification_method.toCellString(); // mandatory
    sl.push_back(s);

    for (const auto& sample : md.sample)
    {
      const MzTabSampleMetaData& msmd = sample.second;

      if (!msmd.description.isNull())
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-description\t" + msmd.description.toCellString();
        sl.push_back(s);
      }

      for (const auto& species : msmd.species)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-species[" + String(species.first) + "]\t" + species.second.toCellString();
        sl.push_back(s);
      }

      for (const auto& tissue : msmd.tissue)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-tissue[" + String(tissue.first) + "]\t" + tissue.second.toCellString();
        sl.push_back(s);
      }

      for (const auto& cell_type : msmd.cell_type)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-cell_type[" + String(cell_type.first) + "]\t" + cell_type.second.toCellString();
        sl.push_back(s);
      }

      for (auto const& disease : msmd.disease)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-disease[" + String(disease.first) + "]\t" + disease.second.toCellString();
        sl.push_back(s);
      }

      for (const auto& custom : msmd.custom)
      {
        String s = "MTD\tsample[" + String(sample.first) + "]-custom[" + String(custom.first) + "]\t" + custom.second.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto& ms_run : md.ms_run)
    {
      const MzTabMMSRunMetaData& rmmd = ms_run.second;

      String s = "MTD\tms_run[" + String(ms_run.first) + "]-location\t" + rmmd.location.toCellString(); // mandatory
      sl.push_back(s);

      if (!rmmd.instrument_ref.isNull())
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-instrument_ref\t" + rmmd.instrument_ref.toCellString();
        sl.push_back(s);
      }

      if (!rmmd.format.isNull())
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-format\t" + rmmd.format.toCellString();
        sl.push_back(s);
      }

      for (const auto& fragmentation_method : rmmd.fragmentation_method)
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-fragmentation_method[" + String(fragmentation_method.first) + "]\t" + fragmentation_method.second.toCellString();
        sl.push_back(s);
      }

      for (const auto& scan_polarity : rmmd.scan_polarity)
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-scan_polarity[" + String(scan_polarity.first) + "]\t" + scan_polarity.second.toCellString(); // mandatory
        sl.push_back(s);
      }

      if (!rmmd.hash.isNull())
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-hash\t" + rmmd.hash.toCellString();
        sl.push_back(s);
      }

      if (!rmmd.hash_method.isNull())
      {
        String s = "MTD\tms_run[" + String(ms_run.first) + "]-hash_method\t" + rmmd.hash_method.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto& assay : md.assay)
    {
      const MzTabMAssayMetaData& amd = assay.second;

      String name = "MTD\tassay[" + String(assay.first) + "]-name\t" + amd.name.toCellString(); // mandatory
      sl.push_back(name);

      for (const auto& custom : amd.custom)
      {
        String s = "MTD\tms_run[" + String(assay.first) + "]-custom[" + String(custom.first) + "]\t" + custom.second.toCellString();
        sl.push_back(s);
      }

      if (!amd.external_uri.isNull())
      {
        String s = "MTD\tassay[" + String(assay.first) + "]-external_uri\t" + amd.external_uri.toCellString();
        sl.push_back(s);
      }

      if (!amd.sample_ref.isNull())
      {
        String s = "MTD\tassay[" + String(assay.first) + "]-sample_ref\tsample[" + amd.sample_ref.toCellString() + "]";
        sl.push_back(s);
      }

      String ms_run_ref = "MTD\tassay[" + String(assay.first) + "]-ms_run_ref\tms_run[" + amd.ms_run_ref.toCellString() + "]"; // mandatory
      sl.push_back(ms_run_ref);
    }

    for (const auto& sv : md.study_variable)
    {
      const MzTabMStudyVariableMetaData& svmd = sv.second;

      String name = "MTD\tstudy_variable[" + String(sv.first) + "]-name\t" + svmd.name.toCellString(); // mandatory
      sl.push_back(name);

      std::vector<MzTabString> strings;
      MzTabStringList refs_string;
      for (const auto& ref : svmd.assay_refs)
      {
        strings.emplace_back(MzTabString("assay[" + std::to_string(ref) + ']'));
      }
      refs_string.set(strings);
      String refs = "MTD\tstudy_variable[" + String(sv.first) + "]-assay_refs\t" + refs_string.toCellString(); // mandatory
      sl.push_back(refs);

      if (!svmd.average_function.isNull())
      {
        String s = "MTD\tstudy_variable[" + String(sv.first) + "]-average_function\t" + svmd.average_function.toCellString();
        sl.push_back(s);
      }

      if (!svmd.variation_function.isNull())
      {
        String s = "MTD\tstudy_variable[" + String(sv.first) + "]-variation_function\t" + svmd.variation_function.toCellString();
        sl.push_back(s);
      }

      String description = "MTD\tstudy_variable[" + String(sv.first) + "]-description\t" + svmd.description.toCellString(); // mandatory
      sl.push_back(description);

      for (const auto& factor : svmd.factors.get())
      {
        String s = "MTD\tstudy_variable[" + String(sv.first) + "]-factors\t" + factor.toCellString();
        sl.push_back(s);
      }
    }

    for (const auto& custom : md.custom)
    {
      String s = "MTD\tcustom[" + String(custom.first) + "]\t" + custom.second.toCellString();
      sl.push_back(s);
    }

    for (const auto& cv : md.cv)
    {
      const MzTabCVMetaData& cvmd = cv.second;

      String label = "MTD\tcv[" + String(cv.first) + "]-label\t" + cvmd.label.toCellString(); // mandatory
      sl.push_back(label);

      String full_name = "MTD\tcv[" + String(cv.first) + "]-full_name\t" + cvmd.full_name.toCellString();  // mandatory
      sl.push_back(full_name);

      String version = "MTD\tcv[" + String(cv.first) + "]-version\t" + cvmd.version.toCellString();  // mandatory
      sl.push_back(version);

      String url = "MTD\tcv[" + String(cv.first) + "]-uri\t" + cvmd.url.toCellString();  // mandatory
      sl.push_back(url);
    }

    for (const auto& db : md.database)
    {
      MzTabMDatabaseMetaData dbmd = db.second;

      String database = "MTD\tdatabase[" + String(db.first) + "]-database\t" + dbmd.database.toCellString();  // mandatory
      sl.push_back(database);

      String prefix = "MTD\tdatabase[" + String(db.first) + "]-prefix\t" + dbmd.prefix.toCellString();  // mandatory
      sl.push_back(prefix);

      String version = "MTD\tdatabase[" + String(db.first) + "]-version\t" + dbmd.version.toCellString();  // mandatory
      sl.push_back(version);

      String uri = "MTD\tdatabase[" + String(db.first) + "]-uri\t" + dbmd.uri.toCellString();  // mandatory
      sl.push_back(uri);
    }

    for (const auto& agent : md.derivatization_agent)
    {
      String s = "MTD\tderivatization_agent[" + String(agent.first) + "]-uri\t" + agent.second.toCellString();
      sl.push_back(s);
    }

    sl.push_back(String("MTD\tsmall_molecule_quantification_unit\t") + md.small_molecule_quantification_unit.toCellString()); // mandatory
    // TODO: Have to check if feature section is reported?
    sl.push_back(String("MTD\tsmall_molecule_feature_quantification_unit\t") + md.small_molecule_feature_quantification_unit.toCellString()); // mandatory (feature section)

    if (!md.small_molecule_identification_reliability.isNull())
    {
      sl.push_back(String("MTD\tsmall_molecule_identification_reliability\t") + md.small_molecule_identification_reliability.toCellString()); // mandatory
    }

    for (const auto& id_conf : md.id_confidence_measure)
    {
      String s = "MTD\tid_confidence_measure[" + String(id_conf.first) + "]\t" + id_conf.second.toCellString(); // mandatory
      sl.push_back(s);
    }

    for (const auto& csm : md.colunit_small_molecule)
    {
      String s = "MTD\tcolunit_small_molecule\t" + csm.toCellString();
      sl.push_back(s);
    }

    for (const auto& csmf : md.colunit_small_molecule_feature)
    {
      String s = "MTD\tcolunit_small_molecule\t" + csmf.toCellString();
      sl.push_back(s);
    }

    for (const auto& csme : md.colunit_small_molecule_evidence)
    {
      String s = "MTD\tcolunit_small_molecule\t" + csme.toCellString();
      sl.push_back(s);
    }
  }

  String MzTabMFile::generateMzTabSmallMoleculeHeader_(const MzTabMMetaData& meta, const std::vector<String>& optional_columns, size_t& n_columns) const
  {
    StringList header;
    header.emplace_back("SMH");
    header.emplace_back("SML_ID");
    header.emplace_back("SMF_ID_REFS");
    header.emplace_back("database_identifier");
    header.emplace_back("chemical_formula");
    header.emplace_back("smiles");
    header.emplace_back("inchi");
    header.emplace_back("chemical_name");
    header.emplace_back("uri");
    header.emplace_back("theoretical_neutral_mass");
    header.emplace_back("adduct_ions");
    header.emplace_back("reliability");
    header.emplace_back("best_id_confidence_measure");
    header.emplace_back("best_id_confidence_value");

    for (const auto& a : meta.assay)
    {
      header.emplace_back(String("abundance_assay[") + String(a.first) + String("]"));
    }

    for (const auto& a : meta.study_variable)
    {
      header.emplace_back(String("abundance_study_variable[") + String(a.first) + String("]"));
    }

    for (const auto& a : meta.study_variable)
    {
      header.emplace_back(String("abundance_variation_study_variable[") + String(a.first) + String("]"));
    }


    // abundance_assay[1-n]
    // abundance_study_variable[1-n]
    // abundance_variation_study_variable [1-n]
    // opt_{identifier}_*



    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

//  String MzTabMFile::generateMzTabSmallMoleculeSectionRow_() const
//  {
//
//  }

  String MzTabMFile::generateMzTabSmallMoleculeFeatureHeader_(const MzTabMSmallMoleculeFeatureSectionRow& row, const std::vector<String>& optional_columns, const MzTabMMetaData& meta, size_t& n_columns) const
  {
   StringList header;
   header.emplace_back("SFH ");
   header.emplace_back("SMF_ID ");
   header.emplace_back("SME_ID_REFS ");
   header.emplace_back("SME_ID_REF_ambiguity_code ");
   header.emplace_back("adduct_ion ");
   header.emplace_back("isotopomer ");
   header.emplace_back("exp_mass_to_charge ");
   header.emplace_back("charge ");
   header.emplace_back("retention_time_in_seconds ");
   header.emplace_back("retention_time_in_seconds_start ");
   header.emplace_back("retention_time_in_seconds_end ");
   // abundance_assay[1-n]
   // opt_{identifier}_*




   std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
   n_columns = header.size();
   return ListUtils::concatenate(header, "\t");
  }

  //String MzTabMFile::generateMzTabSmallMoleculeFeatureSectionRow_() const;

  String MzTabMFile::generateMzTabSmallMoleculeEvidenceHeader_(const MzTabMSmallMoleculeEvidenceSectionRow & row, const std::vector<String>& optional_columns, const MzTabMMetaData& meta, size_t& n_columns) const
  {
    StringList header;
    header.emplace_back("SEH");
    header.emplace_back("SME_ID");
    header.emplace_back("evidence_input_id");
    header.emplace_back("database_identifier");
    header.emplace_back("chemical_formula");
    header.emplace_back("smiles");
    header.emplace_back("inchi");
    header.emplace_back("chemical_name");
    header.emplace_back("uri");
    header.emplace_back("derivatized_form");
    header.emplace_back("adduct_ion");
    header.emplace_back("exp_mass_to_charge");

    header.emplace_back("charge");
    header.emplace_back("theoretical_mass_to_charge");
    header.emplace_back("spectra_ref");

    header.emplace_back("identification_method");
    header.emplace_back("ms_level");
    // id_confidence_measure[1-n]

    header.emplace_back("rank");
    // opt_{identifier}_*



    std::copy(optional_columns.begin(), optional_columns.end(), std::back_inserter(header));
    n_columns = header.size();
    return ListUtils::concatenate(header, "\t");
  }

  //String MzTabMFile::generateMzTabSmallMoleculeEvidenceSectionRow_() const;

  void MzTabMFile::store(const String& filename, const MzTabM& mztab_m) const
  {
    OPENMS_LOG_INFO << "exporting identification data: \"" << filename << "\" to MzTab-M: " << std::endl;

    if (!(FileHandler::hasValidExtension(filename, FileTypes::MZTAB) || FileHandler::hasValidExtension(filename, FileTypes::TSV)))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '"
                                                                                                    + FileTypes::typeToName(FileTypes::MZTAB) + "' or '" + FileTypes::typeToName(FileTypes::TSV) + "'");
    }

    StringList out;
    generateMzTabMMetaDataSection_(mztab_m.getMetaData(), out);

    const MzTabMSmallMoleculeSectionRows& sm_section = mztab_m.getMSmallMoleculeSectionRows();
    const MzTabMSmallMoleculeFeatureSectionRows& feature_section = mztab_m.getMSmallMoleculeFeatureSectionRows();
    const MzTabMSmallMoleculeEvidenceSectionRows& evidence_section = mztab_m.getMSmallMoleculeEvidenceSectionRows();



    // TODO: add to later on - empty rows / comment rows
    TextFile tmp_out;
    for (TextFile::ConstIterator it = out.begin(); it != out.end(); ++it)
    {
      tmp_out.addLine(*it);
    }
    tmp_out.store(filename);

  }

}

#pragma clang diagnostic pop
