// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <unordered_map>

using namespace std;

namespace OpenMS
{

  IdXMLFile::IdXMLFile() :
    XMLHandler("", "1.5"),
    XMLFile("/SCHEMAS/IdXML_1_5.xsd", "1.5"),
    last_meta_(nullptr),
    document_id_(),
    prot_id_in_run_(false)
  {
  }

  void IdXMLFile::load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids)
  {
    String document_id;
    load(filename, protein_ids, peptide_ids, document_id);
  }

  void IdXMLFile::load(const String& filename, std::vector<ProteinIdentification>& protein_ids,
                       std::vector<PeptideIdentification>& peptide_ids, String& document_id)
  {
    startProgress(0, 0, "Loading idXML");
    //Filename for error messages in XMLHandler
    file_ = filename;

    protein_ids.clear();
    peptide_ids.clear();

    prot_ids_ = &protein_ids;
    pep_ids_ = &peptide_ids;
    document_id_ = &document_id;

    parse_(filename, this);

    //reset members
    prot_ids_ = nullptr;
    pep_ids_ = nullptr;
    last_meta_ = nullptr;
    parameters_.clear();
    param_ = ProteinIdentification::SearchParameters();
    id_ = "";
    prot_id_ = ProteinIdentification();
    pep_id_ = PeptideIdentification();
    prot_hit_ = ProteinHit();
    pep_hit_ = PeptideHit();
    proteinid_to_accession_.clear();

    endProgress();
  }

  void IdXMLFile::store(const String& filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids, const String& document_id)
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::IDXML))
    {
      throw Exception::UnableToCreateFile(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          filename,
          "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::IDXML) + "'");
    }

    //set filename for the handler. Just in case (e.g. when fatalError function is used).
    file_ = filename;

    //open stream
    std::ofstream os(filename.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    startProgress(0, peptide_ids.size(), "Storing idXML");

    os.precision(writtenDigits<double>(0.0));

    // write header
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    os << "<?xml-stylesheet type=\"text/xsl\" href=\"https://www.openms.de/xml-stylesheet/IdXML.xsl\" ?>\n";
    os << "<IdXML version=\"" << getVersion() << "\"";
    if (!document_id.empty())
    {
      os << " id=\"" << document_id << "\"";
    }
    os << " xsi:noNamespaceSchemaLocation=\"https://www.openms.de/xml-schema/IdXML_1_5.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

    // look up different search parameters
    std::vector<ProteinIdentification::SearchParameters> params;
    for (std::vector<ProteinIdentification>::const_iterator it = protein_ids.begin(); it != protein_ids.end(); ++it)
    {
      if (find(params.begin(), params.end(), it->getSearchParameters()) == params.end())
      {
        params.push_back(it->getSearchParameters());
      }
    }

    // write search parameters
    for (Size i = 0; i != params.size(); ++i)
    {
      os << "\t<SearchParameters "
         << "id=\"SP_" << i << "\" "
         << "db=\"" << writeXMLEscape(params[i].db) << "\" "
         << "db_version=\"" << writeXMLEscape(params[i].db_version) << "\" "
         << "taxonomy=\"" << writeXMLEscape(params[i].taxonomy) << "\" ";
      if (params[i].mass_type == ProteinIdentification::MONOISOTOPIC)
      {
        os << "mass_type=\"monoisotopic\" ";
      }
      else if (params[i].mass_type == ProteinIdentification::AVERAGE)
      {
        os << "mass_type=\"average\" ";
      }
      os << "charges=\"" << params[i].charges << "\" ";
      String enzyme_name = params[i].digestion_enzyme.getName();
      os << "enzyme=\"" << enzyme_name.toLower() << "\" ";
      String precursor_unit = params[i].precursor_mass_tolerance_ppm ? "true" : "false";
      String peak_unit = params[i].fragment_mass_tolerance_ppm ? "true" : "false";

      os << "missed_cleavages=\"" << params[i].missed_cleavages << "\" "
         << "precursor_peak_tolerance=\"" << params[i].precursor_mass_tolerance << "\" ";
      os << "precursor_peak_tolerance_ppm=\"" << precursor_unit << "\" ";
      os << "peak_mass_tolerance=\"" << params[i].fragment_mass_tolerance << "\" ";
      os << "peak_mass_tolerance_ppm=\"" << peak_unit << "\" ";
      os << ">\n";

      //modifications
      for (Size j = 0; j != params[i].fixed_modifications.size(); ++j)
      {
        os << "\t\t<FixedModification name=\"" << writeXMLEscape(params[i].fixed_modifications[j]) << "\" />\n";
        //Add MetaInfo, when modifications has it (Andreas)
      }
      for (Size j = 0; j != params[i].variable_modifications.size(); ++j)
      {
        os << "\t\t<VariableModification name=\"" << writeXMLEscape(params[i].variable_modifications[j]) << "\" />\n";
        //Add MetaInfo, when modifications has it (Andreas)
      }

      writeUserParam_("UserParam", os, params[i], 4);
      if (params[i].enzyme_term_specificity != EnzymaticDigestion::SPEC_UNKNOWN)
      {
        os << "\t\t\t\t<UserParam name=\"EnzymeTermSpecificity\" type=\"string\" value=\"" << EnzymaticDigestion::NamesOfSpecificity[params[i].enzyme_term_specificity] << "\" />\n";
      }

      os << "\t</SearchParameters>\n";
    }

    //empty search parameters
    if (params.empty())
    {
      os << "<SearchParameters charges=\"+0, +0\" id=\"ID_1\" db_version=\"0\" mass_type=\"monoisotopic\" peak_mass_tolerance=\"0.0\" precursor_peak_tolerance=\"0.0\" db=\"Unknown\"/>\n";
    }

    // throws if protIDs are not unique, i.e. PeptideIDs will be randomly assigned (bad!)
    checkUniqueIdentifiers_(protein_ids);

    UInt prot_count = 0;
    std::unordered_map<string, UInt> accession_to_id;
    size_t protein_count{0};
    for (const auto& pi : protein_ids)
    {
      protein_count += pi.getHits().size();
    }
    accession_to_id.reserve(protein_count); // expect this many keys (avoid rehashing)

    // identifiers of protein identifications that are already written
    std::vector<String> done_identifiers;

    // write ProteinIdentification Runs
    for (Size i = 0; i < protein_ids.size(); ++i)
    {
      done_identifiers.push_back(protein_ids[i].getIdentifier());

      os << "\t<IdentificationRun ";
      os << "date=\"" << protein_ids[i].getDateTime().getDate() << "T" << protein_ids[i].getDateTime().getTime() << "\" ";
      os << "search_engine=\"" << writeXMLEscape(protein_ids[i].getSearchEngine()) << "\" ";
      os << "search_engine_version=\"" << writeXMLEscape(protein_ids[i].getSearchEngineVersion()) << "\" ";
      // identifier
      for (Size j = 0; j != params.size(); ++j)
      {
        if (params[j] == protein_ids[i].getSearchParameters())
        {
          os << "search_parameters_ref=\"SP_" << j << "\" ";
          break;
        }
      }
      os << ">\n";
      os << "\t\t<ProteinIdentification ";
      os << "score_type=\"" << writeXMLEscape(protein_ids[i].getScoreType()) << "\" ";
      if (protein_ids[i].isHigherScoreBetter())
      {
        os << "higher_score_better=\"true\" ";
      }
      else
      {
        os << "higher_score_better=\"false\" ";
      }
      os << "significance_threshold=\"" << protein_ids[i].getSignificanceThreshold() << "\" >\n";

      // write protein hits
      size_t hit_count { protein_ids[i].getHits().size() };
      for (Size j = 0; j < hit_count; ++j)
      {
        os << "\t\t\t<ProteinHit "
           << "id=\"PH_" << String(prot_count) << "\" "
           << "accession=\"" << writeXMLEscape(protein_ids[i].getHits()[j].getAccession()) << "\" "
           << "score=\"" << String(protein_ids[i].getHits()[j].getScore()) << "\" ";
        accession_to_id[protein_ids[i].getHits()[j].getAccession()] = prot_count;
        ++prot_count;

        double coverage = protein_ids[i].getHits()[j].getCoverage();
        if (coverage != ProteinHit::COVERAGE_UNKNOWN)
        {
          os << "coverage=\"" << String(coverage) << "\" ";
        }

        os << "sequence=\"" << writeXMLEscape(protein_ids[i].getHits()[j].getSequence()) << "\" >\n";
        writeUserParam_("UserParam", os, protein_ids[i].getHits()[j], 4);
        os << "\t\t\t</ProteinHit>\n";
      }

      // add ProteinGroup info to metavalues (hack)
      MetaInfoInterface meta = protein_ids[i];
      addProteinGroups_(meta, protein_ids[i].getProteinGroups(),
                        "protein_group", accession_to_id, STORE);
      addProteinGroups_(meta, protein_ids[i].getIndistinguishableProteins(),
                        "indistinguishable_proteins", accession_to_id, STORE);
      writeUserParam_("UserParam", os, meta, 3);

      os << "\t\t</ProteinIdentification>\n";

      //write PeptideIdentifications

      Size count_wrong_id(0);
      Size count_empty(0);

      for (Size l = 0; l < peptide_ids.size(); ++l)
      {
        setProgress(l);

        if (peptide_ids[l].getIdentifier() != protein_ids[i].getIdentifier())
        {
          ++count_wrong_id;
          continue;
        }
        else if (peptide_ids[l].getHits().empty())
        {
          ++count_empty;
          continue;
        }

        os << "\t\t<PeptideIdentification "
           << "score_type=\"" << writeXMLEscape(peptide_ids[l].getScoreType()) << "\" ";
        if (peptide_ids[l].isHigherScoreBetter())
        {
          os << "higher_score_better=\"true\" ";
        }
        else
        {
          os << "higher_score_better=\"false\" ";
        }
        os << "significance_threshold=\"" << String(peptide_ids[l].getSignificanceThreshold()) << "\" ";
        // mz
        if (peptide_ids[l].hasMZ())
        {
          os << "MZ=\"" << String(peptide_ids[l].getMZ()) << "\" ";
        }
        // rt
        if (peptide_ids[l].hasRT())
        {
          os << "RT=\"" << String(peptide_ids[l].getRT()) << "\" ";
        }
        // spectrum_reference
        const DataValue& dv = peptide_ids[l].getMetaValue("spectrum_reference");
        if (dv != DataValue::EMPTY)
        {
          os << "spectrum_reference=\"" << writeXMLEscape(dv.toString()) << "\" ";
        }
        os << ">\n";

        // write peptide hits
        std::vector<String> protein_accessions;

        // copy current hit
        PeptideIdentification pep_id = peptide_ids[l];

        // sort by score
        pep_id.sort();
        const vector<PeptideHit>& pep_hits = pep_id.getHits();

        for (const PeptideHit& p_hit : pep_hits)
        {
          os << "\t\t\t<PeptideHit"
             << " score=\"" << String(p_hit.getScore()) << "\""
             << " sequence=\"" << writeXMLEscape(p_hit.getSequence().toString()) << "\""
             << " charge=\"" << String(p_hit.getCharge()) << "\"";

          const std::vector<PeptideEvidence>& pes = p_hit.getPeptideEvidences();

          createFlankingAAXMLString_(pes, os);
          createPositionXMLString_(pes, os);

          // Extract all protein accessions.
          // Note: protein accessions correspond to neighboring AAs and start/end
          // positions, so we have to keep the same order and allow duplicates
          // (for peptides matching multiple times in the same protein)

          protein_accessions.clear();
          for (vector<PeptideEvidence>::const_iterator pe = pes.begin(); pe != pes.end(); ++pe)
          {
            const String& protein_accession = pe->getProteinAccession();

            // empty accessions are not written out (legacy code)
            if (!protein_accession.empty())
            {
              const auto acc = accession_to_id.find(protein_accession);
              if (acc != accession_to_id.end())
              {
                protein_accessions.emplace_back("PH_" + String(acc->second));
              }
              else
              {
                throw Exception::ElementNotFound(
                    __FILE__,
                    __LINE__,
                    OPENMS_PRETTY_FUNCTION,
                    "No accession " + protein_accession + " found in run '" + protein_ids[i].getIdentifier() +
                    "' for PSM " + p_hit.getSequence().toString() + "_" + String(p_hit.getCharge()) +
                    ". Please contact the maintainer of this tool e.g. on GitHub as this should not happen.");
              }
            }
          }

          if (!protein_accessions.empty())
          {
            os << " protein_refs=\"" << ListUtils::concatenate(protein_accessions, " ") << "\"";
          }

          os << " >\n";
          writeFragmentAnnotations_("UserParam", os, p_hit.getPeakAnnotations(), 4);
          writeUserParam_("UserParam", os, p_hit, 4);

          // write out the (optional) peptide prophet / interprophet results as UserParams
          {
            int k = 0;
            for (std::vector<PeptideHit::PepXMLAnalysisResult>::const_iterator ar_it = p_hit.getAnalysisResults().begin();
                ar_it != p_hit.getAnalysisResults().end(); ++ar_it, ++k)
            {
              os << "\t\t\t\t<UserParam type=\"string\" name=\"_ar_" << String(k) << "_score_type\" value=\"" << ar_it->score_type << "\"/>" << "\n";
              os << "\t\t\t\t<UserParam type=\"float\" name=\"_ar_" << String(k) << "_score\" value=\"" << String(ar_it->main_score) << "\"/>" << "\n";
              if (!ar_it->sub_scores.empty())
              {
                for (std::map<String, double>::const_iterator subscore_it = ar_it->sub_scores.begin();
                    subscore_it != ar_it->sub_scores.end(); ++subscore_it)
                {
                  os << "\t\t\t\t<UserParam type=\"float\" name=\"_ar_" << String(k) << "_subscore_" << subscore_it->first <<"\" value=\"" << String(subscore_it->second) << "\"/>" << "\n";
                }
              }
            }

          }
          os << "\t\t\t</PeptideHit>\n";
        }

        // do not write "spectrum_reference" since it is written as attribute already
        pep_id.removeMetaValue("spectrum_reference");
        writeUserParam_("UserParam", os, pep_id, 3);
        os << "\t\t</PeptideIdentification>\n";
      }

      os << "\t</IdentificationRun>\n";

      // on more than one protein Ids (=runs) there must be wrong mappings and the message would be useless. However, a single run should not have wrong mappings!
      if (count_wrong_id && protein_ids.size() == 1) OPENMS_LOG_WARN << "Omitted writing of " << count_wrong_id << " peptide identifications due to wrong protein mapping." << std::endl;
      if (count_empty) OPENMS_LOG_WARN << "Omitted writing of " << count_empty << " peptide identifications due to empty hits." << std::endl;
    }

    // empty protein ids  parameters
    if (protein_ids.empty())
    {
      os << "<IdentificationRun date=\"1900-01-01T01:01:01.0Z\" search_engine=\"Unknown\" search_parameters_ref=\"ID_1\" search_engine_version=\"0\"/>\n";
    }

    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      if (find(done_identifiers.begin(), done_identifiers.end(), peptide_ids[i].getIdentifier()) == done_identifiers.end())
      {
        warning(STORE, String("Omitting peptide identification because of missing ProteinIdentification with identifier '") + peptide_ids[i].getIdentifier() + "' while writing '" + filename + "'!");
      }
    }
    // write footer
    os << "</IdXML>\n";

    // close stream
    os.close();

    endProgress();

    //reset members
    prot_ids_ = nullptr;
    pep_ids_ = nullptr;
    last_meta_ = nullptr;
    parameters_.clear();
    param_ = ProteinIdentification::SearchParameters();
    id_ = "";
    prot_id_ = ProteinIdentification();
    pep_id_ = PeptideIdentification();
    prot_hit_ = ProteinHit();
    pep_hit_ = PeptideHit();
    proteinid_to_accession_.clear();
  }

  void IdXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
    String tag = sm_.convert(qname);

    //START
    if (tag == "IdXML")
    {
      //check file version against schema version
      String file_version = "";
      prot_id_in_run_ = false;

      optionalAttributeAsString_(file_version, attributes, "version");
      if (file_version.empty())
      {
        file_version = "1.0";  //default version is 1.0
      }
      if (file_version.toDouble() > version_.toDouble())
      {
        warning(LOAD, "The XML file (" + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefined program behavior.");
      }

      //document id
      String document_id = "";
      optionalAttributeAsString_(document_id, attributes, "id");
      (*document_id_) = document_id;
    }
    //SEARCH PARAMETERS
    else if (tag == "SearchParameters")
    {
      //store id
      id_ =  attributeAsString_(attributes, "id");

      //reset parameters
      param_ = ProteinIdentification::SearchParameters();

      //load parameters
      param_.db = attributeAsString_(attributes, "db");
      param_.db_version = attributeAsString_(attributes, "db_version");

      optionalAttributeAsString_(param_.taxonomy, attributes, "taxonomy");
      param_.charges = attributeAsString_(attributes, "charges");
      optionalAttributeAsUInt_(param_.missed_cleavages, attributes, "missed_cleavages");
      param_.fragment_mass_tolerance = attributeAsDouble_(attributes, "peak_mass_tolerance");

      String peak_unit;
      optionalAttributeAsString_(peak_unit, attributes, "peak_mass_tolerance_ppm");
      param_.fragment_mass_tolerance_ppm = peak_unit == "true" ? true : false;

      param_.precursor_mass_tolerance = attributeAsDouble_(attributes, "precursor_peak_tolerance");
      String precursor_unit;
      optionalAttributeAsString_(precursor_unit, attributes, "precursor_peak_tolerance_ppm");
      param_.precursor_mass_tolerance_ppm = precursor_unit == "true" ? true : false;

      //mass type
      String mass_type = attributeAsString_(attributes, "mass_type");
      if (mass_type == "monoisotopic")
      {
        param_.mass_type = ProteinIdentification::MONOISOTOPIC;
      }
      else if (mass_type == "average")
      {
        param_.mass_type = ProteinIdentification::AVERAGE;
      }
      //enzyme
      String enzyme;
      optionalAttributeAsString_(enzyme, attributes, "enzyme");
      if (ProteaseDB::getInstance()->hasEnzyme(enzyme))
      {
        param_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme));
      }
      last_meta_ = &param_;
    }
    else if (tag == "FixedModification")
    {
      param_.fixed_modifications.push_back(attributeAsString_(attributes, "name"));
      //change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
      last_meta_ = nullptr;
    }
    else if (tag == "VariableModification")
    {
      param_.variable_modifications.push_back(attributeAsString_(attributes, "name"));
      //change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
      last_meta_ = nullptr;
    }
    // RUN
    else if (tag == "IdentificationRun")
    {
      pep_id_ = PeptideIdentification();
      prot_id_ = ProteinIdentification();

      prot_id_.setSearchEngine(attributeAsString_(attributes, "search_engine"));
      prot_id_.setSearchEngineVersion(attributeAsString_(attributes, "search_engine_version"));

      //search parameters
      String ref = attributeAsString_(attributes, "search_parameters_ref");
      if (parameters_.find(ref) == parameters_.end())
      {
        fatalError(LOAD, String("Invalid search parameters reference '") + ref + "'");
      }
      prot_id_.setSearchParameters(parameters_[ref]);

      //date
      prot_id_.setDateTime(DateTime::fromString(attributeAsString_(attributes, "date")));

      // set identifier (with UID to make downstream merging of prot_ids possible)
      // Note: technically, it would be preferable to prefix the UID for faster string comparison, but this results in random write-orderings during file store (breaks tests)
      prot_id_.setIdentifier(prot_id_.getSearchEngine() + '_' + attributeAsString_(attributes, "date") + '_' + String(UniqueIdGenerator::getUniqueId()));
    }
    //PROTEINS
    else if (tag == "ProteinIdentification")
    {
      prot_id_.setScoreType(attributeAsString_(attributes, "score_type"));

      //optional significance threshold
      double tmp(0.0);
      optionalAttributeAsDouble_(tmp, attributes, "significance_threshold");
      if (tmp != 0.0)
      {
        prot_id_.setSignificanceThreshold(tmp);
      }

      //score orientation
      prot_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes, "higher_score_better")));

      last_meta_ = &prot_id_;
    }
    else if (tag == "ProteinHit")
    {
      prot_hit_ = ProteinHit();
      String accession = attributeAsString_(attributes, "accession");
      prot_hit_.setAccession(accession);
      prot_hit_.setScore(attributeAsDouble_(attributes, "score"));

      // coverage
      double coverage = -std::numeric_limits<double>::max();
      optionalAttributeAsDouble_(coverage, attributes, "coverage");
      if (coverage != -std::numeric_limits<double>::max())
      {
        prot_hit_.setCoverage(coverage);
      }

      // sequence
      String tmp;
      optionalAttributeAsString_(tmp, attributes, "sequence");
      prot_hit_.setSequence(std::move(tmp));

      last_meta_ = &prot_hit_;

      // insert id and accession to map
      proteinid_to_accession_[attributeAsString_(attributes, "id")] = accession;
    }
    // PEPTIDES
    else if (tag == "PeptideIdentification")
    {
      // check whether a prot id has been given, add "empty" one to list else
      if (!prot_id_in_run_)
      {
        prot_ids_->push_back(prot_id_);
        prot_id_in_run_ = true; // set to true, cause we have created one; will be reset for next run
      }

      //set identifier
      pep_id_.setIdentifier(prot_ids_->back().getIdentifier());

      pep_id_.setScoreType(attributeAsString_(attributes, "score_type"));

      //optional significance threshold
      double tmp(0.0);
      optionalAttributeAsDouble_(tmp, attributes, "significance_threshold");
      if (tmp != 0.0)
      {
        pep_id_.setSignificanceThreshold(tmp);
      }

      //score orientation
      pep_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes, "higher_score_better")));

      //MZ
      double tmp2 = -std::numeric_limits<double>::max();
      optionalAttributeAsDouble_(tmp2, attributes, "MZ");
      if (tmp2 != -std::numeric_limits<double>::max())
      {
        pep_id_.setMZ(tmp2);
      }
      //RT
      tmp2 = -std::numeric_limits<double>::max();
      optionalAttributeAsDouble_(tmp2, attributes, "RT");
      if (tmp2 != -std::numeric_limits<double>::max())
      {
        pep_id_.setRT(tmp2);
      }
      String tmp3;
      optionalAttributeAsString_(tmp3, attributes, "spectrum_reference");
      if (!tmp3.empty())
      {
        pep_id_.setSpectrumReference( tmp3);
      }

      last_meta_ = &pep_id_;
    }
    else if (tag == "PeptideHit")
    {
      pep_hit_ = PeptideHit();
      peptide_evidences_.clear();

      pep_hit_.setCharge(attributeAsInt_(attributes, "charge"));
      pep_hit_.setScore(attributeAsDouble_(attributes, "score"));
      pep_hit_.setSequence(AASequence::fromString(String(attributeAsString_(attributes, "sequence"))));

      //parse optional protein ids to determine accessions
      const XMLCh* refs = attributes.getValue(sm_.convert("protein_refs").c_str());
      if (refs != nullptr)
      {
        String accession_string = sm_.convert(refs);
        accession_string.trim();
        std::vector<String> accessions;
        accession_string.split(' ', accessions);
        if (!accession_string.empty() && accessions.empty())
        {
          accessions.push_back(accession_string);
        }

        for (std::vector<String>::const_iterator it = accessions.begin(); it != accessions.end(); ++it)
        {
          const auto it2 = proteinid_to_accession_.find(*it);
          if (it2 != proteinid_to_accession_.end())
          {
            PeptideEvidence pe;
            pe.setProteinAccession(it2->second);
            peptide_evidences_.push_back(std::move(pe));
          }
          else
          {
            fatalError(LOAD, String("Invalid protein reference '") + *it + "'");
          }
        }
      }

      //aa_before
      String tmp;
      optionalAttributeAsString_(tmp, attributes, "aa_before");

      if (!tmp.empty())
      {
        std::vector<String> parts;
        tmp.split(' ', parts);
        if (peptide_evidences_.size() < parts.size())
        {
          peptide_evidences_.resize(parts.size());
        }

        for (Size i = 0; i != parts.size(); ++i)
        {
          peptide_evidences_[i].setAABefore(parts[i][0]);
        }
      }

      //aa_after
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "aa_after");
      if (!tmp.empty())
      {
        std::vector<String> parts;
        tmp.split(' ', parts);
        if (peptide_evidences_.size() < parts.size())
        {
          peptide_evidences_.resize(parts.size());
        }

        for (Size i = 0; i != parts.size(); ++i)
        {
          peptide_evidences_[i].setAAAfter(parts[i][0]);
        }
      }

      //start
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "start");

      if (!tmp.empty())
      {
        std::vector<String> parts;
        tmp.split(' ', parts);
        if (peptide_evidences_.size() < parts.size())
        {
          peptide_evidences_.resize(parts.size());
        }

        for (Size i = 0; i != parts.size(); ++i)
        {
          peptide_evidences_[i].setStart(parts[i].toInt());
        }
      }

      //end
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "end");
      if (!tmp.empty())
      {
        std::vector<String> parts;
        tmp.split(' ', parts);
        if (peptide_evidences_.size() < parts.size())
        {
          peptide_evidences_.resize(parts.size());
        }

        for (Size i = 0; i != parts.size(); ++i)
        {
          peptide_evidences_[i].setEnd(parts[i].toInt());
        }
      }

      last_meta_ = &pep_hit_;
    }
    // USERPARAM
    else if (tag == "UserParam")
    {
      if (last_meta_ == nullptr)
      {
        fatalError(LOAD, "Unexpected tag 'UserParam'!");
      }

      String name = attributeAsString_(attributes, "name");
      String type = attributeAsString_(attributes, "type");

      // Handle specially encoded pepXML analysis results
      if (name.hasPrefix("_ar_"))
      {
        // must be in PeptideHit (indicated by special _ar_ prefix)
        String sfx = name.substr(4, name.size());
        String val_name = sfx.substr(sfx.find("_") + 1, sfx.size());
        if (val_name.hasPrefix("subscore"))
        {
          String score_name = val_name.substr(val_name.find("_") + 1, val_name.size());
          current_analysis_result_.sub_scores[score_name] = attributeAsDouble_(attributes, "value");
        }
        else if (val_name == "score_type")
        {
          if (!current_analysis_result_.score_type.empty())
          {
            pep_hit_.addAnalysisResults(current_analysis_result_);
          }
          current_analysis_result_.score_type = attributeAsString_(attributes, "value");
        }
        else if (val_name == "score")
        {
          current_analysis_result_.main_score = attributeAsDouble_(attributes, "value");
        }
        return;
      }

      if (type == "int")
      {
        last_meta_->setMetaValue(name, attributeAsInt_(attributes, "value"));
      }
      else if (type == "float")
      {
        last_meta_->setMetaValue(name, attributeAsDouble_(attributes, "value"));
      }
      else if (type == "string")
      {
        String value = (String)attributeAsString_(attributes, "value");

        // TODO: check if we are parsing a peptide hit
        if (name == Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM)
        {
          std::vector<PeptideHit::PeakAnnotation> annotations;
          parseFragmentAnnotation_(value, annotations);
          pep_hit_.setPeakAnnotations(annotations);
          return;
      }
        last_meta_->setMetaValue(name, value);
      }
      else if (type == "intList")
      {
        last_meta_->setMetaValue(name, attributeAsIntList_(attributes, "value"));
      }
      else if (type == "floatList")
      {
        last_meta_->setMetaValue(name, attributeAsDoubleList_(attributes, "value"));
      }
      else if (type == "stringList")
      {
        last_meta_->setMetaValue(name, attributeAsStringList_(attributes, "value"));
      }
      else
      {
        fatalError(LOAD, String("Invalid UserParam type '") + type + "' of parameter '" + name + "'");
      }
    }
  }

  void IdXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    String tag = sm_.convert(qname);

    // START
    if (tag == "IdXML")
    {
      prot_id_in_run_ = false;
    }
    // SEARCH PARAMETERS
    else if (tag == "SearchParameters")
    {
      if (last_meta_->metaValueExists("EnzymeTermSpecificity"))
      {
        String spec = last_meta_->getMetaValue("EnzymeTermSpecificity");
        if (spec != "unknown")
        {
          param_.enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(EnzymaticDigestion::getSpecificityByName(spec));
        }
      }
      last_meta_ = nullptr;
      parameters_[id_] = param_;
    }
    else if (tag == "FixedModification")
    {
      last_meta_ = &param_;
    }
    else if (tag == "VariableModification")
    {
      last_meta_ = &param_;
    }
    // PROTEIN IDENTIFICATIONS
    else if (tag == "ProteinIdentification")
    {
      // post processing of ProteinGroups (hack)
      getProteinGroups_(prot_id_.getProteinGroups(), "protein_group");
      getProteinGroups_(prot_id_.getIndistinguishableProteins(),
                        "indistinguishable_proteins");

      prot_ids_->push_back(prot_id_);
      prot_id_ = ProteinIdentification();
      last_meta_  = nullptr;
      prot_id_in_run_ = true;
    }
    else if (tag == "IdentificationRun")
    {
      if (prot_ids_->empty())
      {
        // add empty <ProteinIdentification> if there was none so far (that's where the IdentificationRun parameters are stored)
        prot_ids_->emplace_back(std::move(prot_id_));
      }
      prot_id_ = ProteinIdentification();
      last_meta_ = nullptr;
      prot_id_in_run_ = false;
    }
    else if (tag == "ProteinHit")
    {
      prot_id_.insertHit(std::move(prot_hit_));
      last_meta_ = &prot_id_;
    }
    //PEPTIDES
    else if (tag == "PeptideIdentification")
    {
      pep_ids_->emplace_back(std::move(pep_id_));
      pep_id_ = PeptideIdentification();
      last_meta_ = nullptr;
    }
    else if (tag == "PeptideHit")
    {
      pep_hit_.setPeptideEvidences(std::move(peptide_evidences_));
      peptide_evidences_.clear(); // clear will reset the vector to a valid, known state

      if (!current_analysis_result_.score_type.empty())
      {
        pep_hit_.addAnalysisResults(current_analysis_result_);
      }
      current_analysis_result_ = PeptideHit::PepXMLAnalysisResult();
      pep_id_.insertHit(std::move(pep_hit_));
      last_meta_ = &pep_id_;
    }
  }

  void IdXMLFile::addProteinGroups_(
    MetaInfoInterface& meta, const std::vector<ProteinIdentification::ProteinGroup>&
    groups, const String& group_name, const std::unordered_map<string, UInt>& accession_to_id,
    XMLHandler::ActionMode mode)
  {
    for (Size g = 0; g < groups.size(); ++g)
    {
      String name = group_name + "_" + String(g);
      if (meta.metaValueExists(name))
      {
        warning(mode, String("Metavalue '") + name + "' already exists. Overwriting...");
      }
      String accessions;
      for (StringList::const_iterator acc_it = groups[g].accessions.begin();
           acc_it != groups[g].accessions.end(); ++acc_it)
      {
        if (acc_it != groups[g].accessions.begin())
        {
          accessions += ",";
        }
        const auto pos = accession_to_id.find(*acc_it);
        if (pos != accession_to_id.end())
        {
          accessions += "PH_" + String(pos->second);
        }
        else
        {
          fatalError(mode, String("Invalid protein reference '") + *acc_it + "'");
        }
      }
      String value = String(groups[g].probability) + "," + accessions;
      meta.setMetaValue(name, value);
    }
  }

  void IdXMLFile::getProteinGroups_(std::vector<ProteinIdentification::ProteinGroup>&
                                    groups, const String& group_name)
  {
    groups.clear();
    Size g_id = 0;
    String current_meta = group_name + "_" + String(g_id);
    StringList values;
    while (last_meta_->metaValueExists(current_meta)) // assumes groups have incremental g_IDs
    {
      // convert to proper ProteinGroup
      ProteinIdentification::ProteinGroup g;
      String(last_meta_->getMetaValue(current_meta)).split(',', values);
      if (values.size() < 2)
      {
        fatalError(LOAD, String("Invalid UserParam for ProteinGroups (not enough values)'"));
      }
      g.probability = values[0].toDouble();
      for (Size i_ind = 1; i_ind < values.size(); ++i_ind)
      {
        g.accessions.push_back(proteinid_to_accession_[values[i_ind]]);
      }
      groups.push_back(std::move(g));
      last_meta_->removeMetaValue(current_meta);
      current_meta = group_name + "_" + String(++g_id);
    }
  }

  std::ostream& IdXMLFile::createFlankingAAXMLString_(const std::vector<PeptideEvidence> & pes, std::ostream& os)
  {
    // Check if information on previous/following aa available. If not, we will not write it out
    bool has_aa_before_information(false);
    bool has_aa_after_information(false);
    String aa_string;

    for (const PeptideEvidence& it : pes)
    {
      if (it.getAABefore() != PeptideEvidence::UNKNOWN_AA)
      {
        has_aa_before_information = true;
      }
      if (it.getAAAfter() != PeptideEvidence::UNKNOWN_AA)
      {
        has_aa_after_information = true;
      }
    }

    if (has_aa_before_information)
    {
      os << " aa_before=\"" << pes.begin()->getAABefore();;
      for (std::vector<PeptideEvidence>::const_iterator it = pes.begin() + 1; it != pes.end(); ++it)
      {
        os << " " << it->getAABefore();
      }
      os << "\"";
    }

    if (has_aa_after_information)
    {
      os << " aa_after=\"" << pes.begin()->getAAAfter();
      for (std::vector<PeptideEvidence>::const_iterator it = pes.begin() + 1; it != pes.end(); ++it)
      {
        os << " " << it->getAAAfter();
      }
      os << "\"";
    }
    return os;
  }

  std::ostream& IdXMLFile::createPositionXMLString_(const std::vector<PeptideEvidence>& pes, std::ostream& os)
  {
    bool has_aa_start_information(false);
    bool has_aa_end_information(false);

    for (const PeptideEvidence& it : pes)
    {
      if (it.getStart() != PeptideEvidence::UNKNOWN_POSITION)
      {
        has_aa_start_information = true;
      }
      if (it.getEnd() != PeptideEvidence::UNKNOWN_POSITION)
      {
        has_aa_end_information = true;
      }
    }

    if (has_aa_start_information || has_aa_end_information)
    {
      if (has_aa_start_information)
      {
        os << " start=\"" << String(pes.begin()->getStart());
        for (std::vector<PeptideEvidence>::const_iterator it = pes.begin() + 1; it != pes.end(); ++it)
        {
          os << " " << String(it->getStart());
        }
        os << "\"";
      }

      if (has_aa_end_information)
      {
        os << " end=\"" << String(pes.begin()->getEnd());
        for (std::vector<PeptideEvidence>::const_iterator it = pes.begin() + 1; it != pes.end(); ++it)
        {
          os << " " << String(it->getEnd());
        }
        os << "\"";
      }
    }
    return os;
  }

  void IdXMLFile::writeFragmentAnnotations_(const String & tag_name, std::ostream & os,
                                            const std::vector<PeptideHit::PeakAnnotation>& annotations, UInt indent)
  {
    String val;
    PeptideHit::PeakAnnotation::writePeakAnnotationsString_(val, annotations);
    if (!val.empty())
    {
      os << String(indent, '\t') << "<" << writeXMLEscape(tag_name) << R"( type="string" name="fragment_annotation" value=")" << writeXMLEscape(val) << "\"/>" << "\n";
    }
  }

  void IdXMLFile::parseFragmentAnnotation_(const String& s, std::vector<PeptideHit::PeakAnnotation> & annotations)
  {
    if (s.empty()) { return; }
    StringList as;
    s.split_quoted('|', as);

    // for each peak annotation: split string and fill fragment annotation entries
    StringList fields;
    for (const auto& pa : as)
    {
      pa.split_quoted(',', fields);
      if (fields.size() != 4)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "Invalid fragment annotation. Four comma-separated fields required. String is: '" + pa + "'");
      }
      PeptideHit::PeakAnnotation fa;
      fa.mz = fields[0].toDouble();
      fa.intensity = fields[1].toDouble();
      fa.charge = fields[2].toInt();
      fa.annotation = fields[3].unquote();
      annotations.push_back(fa);
    }
  }
} // namespace OpenMS
