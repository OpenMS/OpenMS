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
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <fstream>
#include <numeric>
#include <sstream>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

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

  void IdXMLFile::load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids)
  {
    std::string document_id;
    load(filename, protein_ids, peptide_ids, document_id);
  }

  void IdXMLFile::load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids,
                       PeptideIdentificationList& peptide_ids, std::string& document_id)
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

  void IdXMLFile::store(const std::string& filename, const std::vector<ProteinIdentification>& protein_ids, const PeptideIdentificationList& peptide_ids, const std::string& document_id)
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

    // look up different search parameters. SearchParameters::operator== ignores meta values, so runs
    // collapse onto one block; the block carries the modification definitions of all of them.
    const auto definitions = ModificationDefinitionIO::collect(protein_ids, peptide_ids);
    std::vector<ProteinIdentification::SearchParameters> params;
    std::vector<std::set<const ResidueModification*>> params_defs;
    for (std::vector<ProteinIdentification>::const_iterator it = protein_ids.begin(); it != protein_ids.end(); ++it)
    {
      const Size idx = static_cast<Size>(find(params.begin(), params.end(), it->getSearchParameters()) - params.begin());
      if (idx == params.size())
      {
        params.push_back(it->getSearchParameters());
        params_defs.emplace_back();
      }
      const auto d = definitions.find(it->getIdentifier());
      if (d != definitions.end()) params_defs[idx].insert(d->second.begin(), d->second.end());
    }
    for (Size i = 0; i != params.size(); ++i)
    {
      ModificationDefinitionIO::attach(params[i], params_defs[i]);
    }

    // write search parameters
    for (Size i = 0; i != params.size(); ++i)
    {
      os << "\t<SearchParameters "
         << "id=\"SP_" << i << "\" "
         << "db=\"" << writeXMLEscape(params[i].db) << "\" "
         << "db_version=\"" << writeXMLEscape(params[i].db_version) << "\" "
         << "taxonomy=\"" << writeXMLEscape(params[i].taxonomy) << "\" ";
      if (params[i].mass_type == ProteinIdentification::PeakMassType::MONOISOTOPIC)
      {
        os << "mass_type=\"monoisotopic\" ";
      }
      else if (params[i].mass_type == ProteinIdentification::PeakMassType::AVERAGE)
      {
        os << "mass_type=\"average\" ";
      }
      os << "charges=\"" << params[i].charges << "\" ";
      std::string enzyme_name = params[i].digestion_enzyme.getName();
      os << "enzyme=\"" << StringUtils::toLower(enzyme_name) << "\" ";
      std::string precursor_unit = params[i].precursor_mass_tolerance_ppm ? "true" : "false";
      std::string peak_unit = params[i].fragment_mass_tolerance_ppm ? "true" : "false";

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
      std::set<const ResidueModification*> all_defs;
      for (const auto& d : definitions) all_defs.insert(d.second.begin(), d.second.end());
      if (all_defs.empty())
      {
        os << "<SearchParameters charges=\"+0, +0\" id=\"ID_1\" db_version=\"0\" mass_type=\"monoisotopic\" peak_mass_tolerance=\"0.0\" precursor_peak_tolerance=\"0.0\" db=\"Unknown\"/>\n";
      }
      else
      {
        os << "<SearchParameters charges=\"+0, +0\" id=\"ID_1\" db_version=\"0\" mass_type=\"monoisotopic\" peak_mass_tolerance=\"0.0\" precursor_peak_tolerance=\"0.0\" db=\"Unknown\">\n";
        ProteinIdentification::SearchParameters sp;
        ModificationDefinitionIO::attach(sp, all_defs);
        writeUserParam_("UserParam", os, sp, 4);
        os << "\t</SearchParameters>\n";
      }
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

    // assign the peptide identifications to their runs in one pass, keeping the input order
    std::unordered_map<std::string, Size> run_of_identifier;
    run_of_identifier.reserve(protein_ids.size());
    for (Size i = 0; i < protein_ids.size(); ++i)
    {
      run_of_identifier.emplace(protein_ids[i].getIdentifier(), i);
    }
    std::vector<std::vector<Size>> run_peptide_ids(protein_ids.size()); // the ones with hits, written
    std::vector<Size> run_empty_count(protein_ids.size(), 0); // the ones without hits, omitted
    std::vector<Size> without_run; // the ones whose identifier names no run, omitted
    for (Size l = 0; l < peptide_ids.size(); ++l)
    {
      const auto run = run_of_identifier.find(peptide_ids[l].getIdentifier());
      if (run == run_of_identifier.end())
      {
        without_run.push_back(l);
      }
      else if (peptide_ids[l].getHits().empty())
      {
        ++run_empty_count[run->second];
      }
      else
      {
        run_peptide_ids[run->second].push_back(l);
      }
    }

    // write ProteinIdentification Runs
    for (Size i = 0; i < protein_ids.size(); ++i)
    {
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
      
      double significance_threshold = protein_ids[i].getSignificanceThreshold();
      os << "significance_threshold=\"" << StringUtils::toStr(significance_threshold) << "\" >\n";

      // write protein hits
      size_t hit_count { protein_ids[i].getHits().size() };
      for (Size j = 0; j < hit_count; ++j)
      {
        os << "\t\t\t<ProteinHit "
           << "id=\"PH_" << StringUtils::toStr(prot_count) << "\" "
           << "accession=\"" << writeXMLEscape(protein_ids[i].getHits()[j].getAccession()) << "\" "
           << "score=\"" << StringUtils::toStr(protein_ids[i].getHits()[j].getScore()) << "\" ";
        accession_to_id[protein_ids[i].getHits()[j].getAccession()] = prot_count;
        ++prot_count;

        double coverage = protein_ids[i].getHits()[j].getCoverage();
        if (coverage != ProteinHit::COVERAGE_UNKNOWN)
        {
          os << "coverage=\"" << StringUtils::toStr(coverage) << "\" ";
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
      //
      // The peptide identifications of the run are formatted in parallel, in blocks of consecutive ones, each block into
      // a string of its own, and the blocks are written in input order as soon as they are ready: the same bytes as
      // formatting them one after the other into os, while at most one block per thread is held in memory. Meta values
      // are read by registry index and their names cached per thread (MetaInfoRegistry takes a process-wide lock for
      // every name), and the hits are visited in the order of PeptideIdentification::sort() instead of sorting a copy.
      // If formatting fails, neither the failing block nor any later one is written, and the error of the first failing
      // block in input order is rethrown.

      const std::vector<Size>& to_write = run_peptide_ids[i];
      const Size count_empty = run_empty_count[i];
      const Size count_wrong_id = peptide_ids.size() - to_write.size() - count_empty;

      const std::string& run_identifier = protein_ids[i].getIdentifier();
      const MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
      // written as attributes of PeptideIdentification, not as UserParam (UInt(-1) if never registered)
      const UInt spectrum_reference_index = registry.getIndex("spectrum_reference");
      const UInt significance_threshold_index = registry.getIndex(Constants::UserParam::SIGNIFICANCE_THRESHOLD);

      // buffers of one thread
      struct Scratch
      {
        std::vector<std::string> names; // registry index -> name; registered names are not empty, so "" is not resolved yet
        std::vector<UInt> keys;
        std::vector<Size> order;
        std::vector<PeptideHit> sorted_hits;
        std::vector<std::string> protein_accessions;
      };

      // writeUserParam_(), leaving out the meta values with index skip_a or skip_b
      const auto write_user_params = [](std::ostream& out, const MetaInfoInterface& meta_info, const std::string& tag_start,
                                        Scratch& scratch, UInt skip_a, UInt skip_b)
      {
        scratch.keys.clear();
        meta_info.getKeys(scratch.keys);
        for (const UInt key : scratch.keys)
        {
          if (key == skip_a || key == skip_b) continue;
          if (key >= scratch.names.size()) scratch.names.resize(static_cast<Size>(key) + 1);
          std::string& name = scratch.names[key];
          if (name.empty()) name = MetaInfoInterface::metaRegistry().getName(key);
          writeUserParamValue_(out, tag_start, name, meta_info.getMetaValue(key));
        }
      };

      const auto write_peptide_identification = [&](std::ostream& out, const PeptideIdentification& pep_id, Scratch& scratch)
      {
        out << "\t\t<PeptideIdentification "
            << "score_type=\"" << writeXMLEscape(pep_id.getScoreType()) << "\" ";
        if (pep_id.isHigherScoreBetter())
        {
          out << "higher_score_better=\"true\" ";
        }
        else
        {
          out << "higher_score_better=\"false\" ";
        }
        out << "significance_threshold=\"" << StringUtils::toStr(pep_id.getSignificanceThreshold()) << "\" ";

        // mz
        if (pep_id.hasMZ())
        {
          out << "MZ=\"" << StringUtils::toStr(pep_id.getMZ()) << "\" ";
        }
        // rt
        if (pep_id.hasRT())
        {
          out << "RT=\"" << StringUtils::toStr(pep_id.getRT()) << "\" ";
        }
        // spectrum_reference
        const DataValue& dv = pep_id.getMetaValue(spectrum_reference_index);
        if (dv != DataValue::EMPTY)
        {
          out << "spectrum_reference=\"" << writeXMLEscape(dv.toString()) << "\" ";
        }
        out << ">\n";

        // write peptide hits, in the order of PeptideIdentification::sort() (a stable sort by score)
        const vector<PeptideHit>* pep_hits = &pep_id.getHits();
        const auto comparator = PeptideIdentification::getScoreComparator(pep_id.isHigherScoreBetter());
        scratch.order.resize(pep_hits->size());
        std::iota(scratch.order.begin(), scratch.order.end(), Size(0));
        if (std::any_of(pep_hits->begin(), pep_hits->end(), [](const PeptideHit& hit) { return std::isnan(hit.getScore()); }))
        {
          // NaN scores are not a strict weak order, so the result depends on the sorting algorithm, which may differ
          // between sorting indices and sorting the hits (e.g. libc++): sort a copy of the hits, as sort() does
          scratch.sorted_hits = *pep_hits;
          std::stable_sort(scratch.sorted_hits.begin(), scratch.sorted_hits.end(), comparator);
          pep_hits = &scratch.sorted_hits;
        }
        else
        {
          std::stable_sort(scratch.order.begin(), scratch.order.end(), [&](Size x, Size y) { return comparator((*pep_hits)[x], (*pep_hits)[y]); });
        }

        for (const Size h : scratch.order)
        {
          const PeptideHit& p_hit = (*pep_hits)[h];
          out << "\t\t\t<PeptideHit"
              << " score=\"" << StringUtils::toStr(p_hit.getScore()) << "\""
              << " sequence=\"" << writeXMLEscape(p_hit.getSequence().toString()) << "\""
              << " charge=\"" << StringUtils::toStr(p_hit.getCharge()) << "\"";

          const std::vector<PeptideEvidence>& pes = p_hit.getPeptideEvidences();

          createFlankingAAXMLString_(pes, out);
          createPositionXMLString_(pes, out);

          // Extract all protein accessions.
          // Note: protein accessions correspond to neighboring AAs and start/end
          // positions, so we have to keep the same order and allow duplicates
          // (for peptides matching multiple times in the same protein)

          scratch.protein_accessions.clear();
          for (vector<PeptideEvidence>::const_iterator pe = pes.begin(); pe != pes.end(); ++pe)
          {
            const std::string& protein_accession = pe->getProteinAccession();

            // empty accessions are not written out (legacy code)
            if (!protein_accession.empty())
            {
              const auto acc = accession_to_id.find(protein_accession);
              if (acc != accession_to_id.end())
              {
                scratch.protein_accessions.emplace_back("PH_" + StringUtils::toStr(acc->second));
              }
              else
              {
                // constructing an OpenMS exception sets the process-wide GlobalExceptionHandler: one thread at a time
                std::exception_ptr not_found;
#pragma omp critical (IdXMLFile_store_exception)
                not_found = std::make_exception_ptr(Exception::ElementNotFound(
                    __FILE__,
                    __LINE__,
                    OPENMS_PRETTY_FUNCTION,
                    "No accession " + protein_accession + " found in run '" + run_identifier +
                    "' for PSM " + p_hit.getSequence().toString() + "_" + StringUtils::toStr(p_hit.getCharge()) +
                    ". Please contact the maintainer of this tool e.g. on GitHub as this should not happen."));
                std::rethrow_exception(not_found);
              }
            }
          }

          if (!scratch.protein_accessions.empty())
          {
            out << " protein_refs=\"" << ListUtils::concatenate(scratch.protein_accessions, " ") << "\"";
          }

          out << " >\n";
          writeFragmentAnnotations_("UserParam", out, p_hit.getPeakAnnotations(), 4);
          write_user_params(out, p_hit, "\t\t\t\t<UserParam type=\"", scratch, UInt(-1), UInt(-1));

          out << "\t\t\t</PeptideHit>\n";
        }

        // do not write "spectrum_reference" or Constants::UserParam::SIGNIFICANCE_THRESHOLD since it is written as attribute already
        write_user_params(out, pep_id, "\t\t\t<UserParam type=\"", scratch, spectrum_reference_index, significance_threshold_index);
        out << "\t\t</PeptideIdentification>\n";
      };

      const Size block_size = 16;
      const SignedSize num_blocks = static_cast<SignedSize>((to_write.size() + block_size - 1) / block_size);
      const std::streamsize precision = os.precision();
      std::exception_ptr error;
      std::atomic<bool> failed(false);

#ifdef _OPENMP
      // at most one thread per block and at most 16: with more, the formatting outruns the serial write into the file and
      // the threads only wait
      const int team_size = static_cast<int>(std::min<SignedSize>({omp_get_max_threads(), 16, std::max<SignedSize>(num_blocks, 1)}));
#endif
#pragma omp parallel if (num_blocks > 1) num_threads(team_size)
      {
        Scratch scratch;
#pragma omp for ordered schedule(dynamic, 1)
        for (SignedSize b = 0; b < num_blocks; ++b)
        {
          const Size begin = static_cast<Size>(b) * block_size;
          const Size end = std::min(to_write.size(), begin + block_size);
          std::string text;
          std::exception_ptr block_error;
          if (!failed.load(std::memory_order_relaxed))
          {
            try
            {
              std::ostringstream block_os;
              block_os.precision(precision);
              for (Size k = begin; k < end; ++k)
              {
                write_peptide_identification(block_os, peptide_ids[to_write[k]], scratch);
              }
              text = block_os.str();
            }
            catch (...)
            {
              block_error = std::current_exception();
            }
          }
#pragma omp ordered
          {
            if (!failed.load(std::memory_order_relaxed))
            {
              if (block_error)
              {
                error = block_error;
                failed.store(true, std::memory_order_relaxed);
              }
              else
              {
                os.write(text.data(), static_cast<std::streamsize>(text.size()));
                setProgress(to_write[end - 1]);
              }
            }
          }
        }
      }
      if (error) std::rethrow_exception(error);

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

    for (const Size l : without_run)
    {
      warning(STORE,std::string("Omitting peptide identification because of missing ProteinIdentification with identifier '") + peptide_ids[l].getIdentifier() + "' while writing '" + filename + "'!");
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

  void IdXMLFile::onStartElement(const char16_t* qname, const Internal::XMLAttributes& attributes)
  {
    std::string tag = sm_.convert(qname);

    //START
    if (tag == "IdXML")
    {
      //check file version against schema version
      std::string file_version;
      prot_id_in_run_ = false;

      optionalAttributeAsString_(file_version, attributes, "version");
      if (file_version.empty())
      {
        file_version = "1.0";  //default version is 1.0
      }
      if (StringUtils::toDouble(file_version) > StringUtils::toDouble(version_))
      {
        warning(LOAD, "The XML file (" + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefined program behavior.");
      }

      //document id
      std::string document_id;
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

      std::string peak_unit;
      optionalAttributeAsString_(peak_unit, attributes, "peak_mass_tolerance_ppm");
      param_.fragment_mass_tolerance_ppm = peak_unit == "true" ? true : false;

      param_.precursor_mass_tolerance = attributeAsDouble_(attributes, "precursor_peak_tolerance");
      std::string precursor_unit;
      optionalAttributeAsString_(precursor_unit, attributes, "precursor_peak_tolerance_ppm");
      param_.precursor_mass_tolerance_ppm = precursor_unit == "true" ? true : false;

      //mass type
      std::string mass_type = attributeAsString_(attributes, "mass_type");
      if (mass_type == "monoisotopic")
      {
        param_.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
      }
      else if (mass_type == "average")
      {
        param_.mass_type = ProteinIdentification::PeakMassType::AVERAGE;
      }
      //enzyme
      std::string enzyme;
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
      std::string ref = attributeAsString_(attributes, "search_parameters_ref");
      if (!parameters_.contains(ref))
      {
        fatalError(LOAD,std::string("Invalid search parameters reference '") + ref + "'");
      }
      prot_id_.setSearchParameters(parameters_[ref]);

      //date
      prot_id_.setDateTime(DateTime::fromString(attributeAsString_(attributes, "date")));

      // set identifier (with UID to make downstream merging of prot_ids possible)
      // Note: technically, it would be preferable to prefix the UID for faster string comparison, but this results in random write-orderings during file store (breaks tests)
      prot_id_.setIdentifier(prot_id_.getSearchEngine() + '_' + attributeAsString_(attributes, "date") + '_' + StringUtils::toStr(UniqueIdGenerator::getUniqueId()));
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
      std::string accession = attributeAsString_(attributes, "accession");
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
      std::string tmp;
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
      std::string tmp3;
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
      pep_hit_.setSequence(AASequence::fromString(std::string(attributeAsString_(attributes, "sequence"))));

      //parse optional protein ids to determine accessions
      const char16_t* refs = attributes.value(sm_.convert("protein_refs").c_str());
      if (refs != nullptr)
      {
        std::string accession_string = sm_.convert(refs);
        StringUtils::trim(accession_string);
        std::vector<std::string> accessions;
        StringUtils::split(accession_string, ' ', accessions);
        if (!accession_string.empty() && accessions.empty())
        {
          accessions.push_back(accession_string);
        }

        for (std::vector<std::string>::const_iterator it = accessions.begin(); it != accessions.end(); ++it)
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
            fatalError(LOAD,std::string("Invalid protein reference '") + *it + "'");
          }
        }
      }

      //aa_before
      std::string tmp;
      optionalAttributeAsString_(tmp, attributes, "aa_before");

      if (!tmp.empty())
      {
        std::vector<std::string> parts;
        StringUtils::split(tmp, ' ', parts);
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
        std::vector<std::string> parts;
        StringUtils::split(tmp, ' ', parts);
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
        std::vector<std::string> parts;
        StringUtils::split(tmp, ' ', parts);
        if (peptide_evidences_.size() < parts.size())
        {
          peptide_evidences_.resize(parts.size());
        }

        for (Size i = 0; i != parts.size(); ++i)
        {
          peptide_evidences_[i].setStart(StringUtils::toInt32(parts[i]));
        }
      }

      //end
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "end");
      if (!tmp.empty())
      {
        std::vector<std::string> parts;
        StringUtils::split(tmp, ' ', parts);
        if (peptide_evidences_.size() < parts.size())
        {
          peptide_evidences_.resize(parts.size());
        }

        for (Size i = 0; i != parts.size(); ++i)
        {
          peptide_evidences_[i].setEnd(StringUtils::toInt32(parts[i]));
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

      std::string name = attributeAsString_(attributes, "name");
      std::string type = attributeAsString_(attributes, "type");

      // Handle specially encoded pepXML analysis results
      if (StringUtils::hasPrefix(name, "_ar_"))
      {
        // must be in PeptideHit (indicated by special _ar_ prefix)
        std::string sfx = StringUtils::substr(name, 4, name.size());
        std::string val_name = StringUtils::substr(sfx, sfx.find("_") + 1, sfx.size());
        if (StringUtils::hasPrefix(val_name, "subscore"))
        {
          std::string score_name = StringUtils::substr(val_name, val_name.find("_") + 1, val_name.size());
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
        std::string value = (std::string)attributeAsString_(attributes, "value");

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
        fatalError(LOAD,std::string("Invalid UserParam type '") + type + "' of parameter '" + name + "'");
      }
    }
  }

  void IdXMLFile::onEndElement(const char16_t* qname)
  {
    std::string tag = sm_.convert(qname);

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
        std::string spec = StringUtils::toStr(last_meta_->getMetaValue("EnzymeTermSpecificity"));
        if (spec != "unknown")
        {
          param_.enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(EnzymaticDigestion::getSpecificityByName(spec));
        }
      }
      last_meta_ = nullptr;
      ModificationDefinitionIO::registerFrom(param_); // before any peptide sequence is parsed
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
    groups, const std::string& group_name, const std::unordered_map<string, UInt>& accession_to_id,
    XMLHandler::ActionMode mode)
  {
    for (Size g = 0; g < groups.size(); ++g)
    {
      std::string name = group_name + "_" + StringUtils::toStr(g);
      if (meta.metaValueExists(name))
      {
        warning(mode,std::string("Metavalue '") + name + "' already exists. Overwriting...");
      }
      std::string accessions;
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
          accessions += "PH_" + StringUtils::toStr(pos->second);
        }
        else
        {
          fatalError(mode,std::string("Invalid protein reference '") + *acc_it + "'");
        }
      }
      std::string value =StringUtils::toStr(groups[g].probability) + "," + accessions;
      meta.setMetaValue(name, value);
    }
  }

  void IdXMLFile::getProteinGroups_(std::vector<ProteinIdentification::ProteinGroup>&
                                    groups, const std::string& group_name)
  {
    groups.clear();
    Size g_id = 0;
    std::string current_meta = group_name + "_" + StringUtils::toStr(g_id);
    StringList values;
    while (last_meta_->metaValueExists(current_meta)) // assumes groups have incremental g_IDs
    {
      // convert to proper ProteinGroup
      ProteinIdentification::ProteinGroup g;
      StringUtils::split(StringUtils::toStr(last_meta_->getMetaValue(current_meta)), ',', values);
      if (values.size() < 2)
      {
        fatalError(LOAD,std::string("Invalid UserParam for ProteinGroups (not enough values)'"));
      }
      g.probability = StringUtils::toDouble(values[0]);
      for (Size i_ind = 1; i_ind < values.size(); ++i_ind)
      {
        g.accessions.push_back(proteinid_to_accession_[values[i_ind]]);
      }
      groups.push_back(std::move(g));
      last_meta_->removeMetaValue(current_meta);
      current_meta = group_name + "_" + StringUtils::toStr(++g_id);
    }
  }

  std::ostream& IdXMLFile::createFlankingAAXMLString_(const std::vector<PeptideEvidence> & pes, std::ostream& os)
  {
    // Check if information on previous/following aa available. If not, we will not write it out
    bool has_aa_before_information(false);
    bool has_aa_after_information(false);
    std::string aa_string;

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
        os << " start=\"" << StringUtils::toStr(pes.begin()->getStart());
        for (std::vector<PeptideEvidence>::const_iterator it = pes.begin() + 1; it != pes.end(); ++it)
        {
          os << " " << StringUtils::toStr(it->getStart());
        }
        os << "\"";
      }

      if (has_aa_end_information)
      {
        os << " end=\"" << StringUtils::toStr(pes.begin()->getEnd());
        for (std::vector<PeptideEvidence>::const_iterator it = pes.begin() + 1; it != pes.end(); ++it)
        {
          os << " " << StringUtils::toStr(it->getEnd());
        }
        os << "\"";
      }
    }
    return os;
  }

  void IdXMLFile::writeFragmentAnnotations_(const std::string & tag_name, std::ostream & os,
                                            const std::vector<PeptideHit::PeakAnnotation>& annotations, UInt indent)
  {
    std::string val;
    PeptideHit::PeakAnnotation::writePeakAnnotationsString_(val, annotations);
    if (!val.empty())
    {
      os << std::string(indent, '\t') << "<" << writeXMLEscape(tag_name) << R"( type="string" name="fragment_annotation" value=")" << writeXMLEscape(val) << "\"/>" << "\n";
    }
  }

  void IdXMLFile::parseFragmentAnnotation_(const std::string& s, std::vector<PeptideHit::PeakAnnotation> & annotations)
  {
    if (s.empty()) { return; }
    StringList as;
    StringUtils::split_quoted(s, "|", as);

    // for each peak annotation: split string and fill fragment annotation entries
    StringList fields;
    for (const auto& pa : as)
    {
      StringUtils::split_quoted(pa, ",", fields);
      if (fields.size() != 4)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "Invalid fragment annotation. Four comma-separated fields required. std::string is: '" + pa + "'");
      }
      PeptideHit::PeakAnnotation fa;
      fa.mz = StringUtils::toDouble(fields[0]);
      fa.intensity = StringUtils::toDouble(fields[1]);
      fa.charge = StringUtils::toInt32(fields[2]);
      fa.annotation = StringUtils::unquote(fields[3]);
      annotations.push_back(fa);
    }
  }
} // namespace OpenMS
