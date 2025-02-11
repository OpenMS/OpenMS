// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Hendrik Weisser $
// $Authors: Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/PepXMLFile.h>

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include <fstream>

using namespace std;

namespace OpenMS
{

  String PepXMLFile::AminoAcidModification::toUnimodLikeString() const
  {
    String desc = "";
    if (massdiff_ >= 0)
    {
      desc += "+" + String(massdiff_);
    }
    else
    {
      desc += String(massdiff_);
    }

    if (!aminoacid_.empty() || !terminus_.empty())
    {
      desc += " (";
      if (!terminus_.empty())
      {
        if (is_protein_terminus_)
        {
          desc += "Protein ";
        }
        {
          String t = terminus_;
          desc += t.toUpper() + "-term";
        }
        if (!aminoacid_.empty())
        {
          desc += " ";
        }
      }
      if (!aminoacid_.empty())
      {
        String a = aminoacid_;
        desc += a.toUpper();
      }
      desc += ")";
    }
    return desc;
  }

  const String& PepXMLFile::AminoAcidModification::getDescription() const
  {
    return registered_mod_->getFullId();
  }

  bool PepXMLFile::AminoAcidModification::isVariable() const
  {
    return is_variable_;
  }

  double PepXMLFile::AminoAcidModification::getMass() const
  {
    return mass_;
  }

  double PepXMLFile::AminoAcidModification::getMassDiff() const
  {
    return massdiff_;
  }

  const String& PepXMLFile::AminoAcidModification::getTerminus() const
  {
    return terminus_;
  }

  const String& PepXMLFile::AminoAcidModification::getAminoAcid() const
  {
    return aminoacid_;
  }

  const ResidueModification* PepXMLFile::AminoAcidModification::getRegisteredMod() const
  {
    return registered_mod_;
  }

  const vector<String>& PepXMLFile::AminoAcidModification::getErrors() const
  {
    return errors_;
  }

  const ResidueModification*
  PepXMLFile::AminoAcidModification::lookupModInPreferredMods_(const vector<const ResidueModification*>& preferred_mods,
                                                       const String& aminoacid,
                                                       double massdiff,
                                                       const String& description,
                                                       const ResidueModification::TermSpecificity term_spec,
                                                       double tolerance)
  {
    for (const auto& pref_mod : preferred_mods)
    {
      if (description == pref_mod->getFullId())
      {
        return pref_mod;
      }
    }
    for (const auto& pref_mod : preferred_mods)
    {
      if ((aminoacid.empty() || aminoacid[0] == pref_mod->getOrigin()) &&
      (term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY || term_spec == pref_mod->getTermSpecificity()))
      {
        if (fabs(massdiff - pref_mod->getDiffMonoMass()) < tolerance)
        {
          return pref_mod;
        }
      }
    }
    return nullptr;
  }

  PepXMLFile::AminoAcidModification::AminoAcidModification(
      const String& aminoacid, const String& massdiff, const String& mass,
      String variable, const String& description, String terminus, const String& protein_terminus,
      const vector<const ResidueModification*>& preferred_fixed_mods,
      const vector<const ResidueModification*>& preferred_var_mods,
      double tolerance)
  {
    if (aminoacid.empty() && terminus.empty())
    {
      throw Exception::MissingInformation(__FILE__,
                                          __LINE__,
                                          OPENMS_PRETTY_FUNCTION,
                                          "Either terminus or amino acid origin or both needs to be set.");
    }

    aminoacid_ = aminoacid;
    massdiff_ = massdiff.toDouble();
    mass_ = mass.toDouble();
    is_variable_ = variable.toLower() == "y";
    description_ = description;
    registered_mod_ = nullptr;
    terminus_ = terminus.toLower();
    is_protein_terminus_ = false;
    term_spec_ = ResidueModification::NUMBER_OF_TERM_SPECIFICITY;

    if (terminus_ == "nc")
    {
      errors_.emplace_back("Warning: value 'nc' for aminoacid terminus not supported."
                    "The modification will be parsed as an unrestricted modification.");
    }

    if (aminoacid_.size() > 1)
    {
      errors_.emplace_back("Warning: Single modification specified for multiple amino acids. This is not supported."
                           "Please split them into one modification per amino acid. Proceeding with first AA...");
    }

    // BIG NOTE: According to the pepXML schema specification, protein terminus is either "c" or "n" if set.
    // BUT: Many tools will put "Y" or "N" there in conjunction with the terminus attribute.
    // Unfortunately there is an overlap for the letter "n". We will try to handle both based on upper/lower case:
    String protein_terminus_lower = protein_terminus;
    protein_terminus_lower = protein_terminus_lower.toLower();

    if (protein_terminus_lower == "y")
    {
      is_protein_terminus_ = true;
    }
    else if (protein_terminus_lower == "c")
    {
      is_protein_terminus_ = true;
      terminus_ = protein_terminus_lower; // protein_terminus takes precedence. I would assume they are the same
    }
    else if (protein_terminus == "n")
    {
      is_protein_terminus_ = true;
      terminus_ = protein_terminus;
    }
    else if (protein_terminus == "N")
    {
      is_protein_terminus_ = false;
    }

    // Now set our internal enum based on the inferred values
    if (terminus_ == "n")
    {
      if (is_protein_terminus_)
      {
        term_spec_ = ResidueModification::PROTEIN_N_TERM;
      }
      else
      {
        term_spec_ = ResidueModification::N_TERM;
      }
    }
    else if (terminus_ == "c")
    {
      if (is_protein_terminus_)
      {
        term_spec_ = ResidueModification::PROTEIN_C_TERM;
      }
      else
      {
        term_spec_ = ResidueModification::C_TERM;
      }
    }

    if (mass_ == massdiff_)
    {
      errors_.emplace_back("Warning: For modification with mass " + mass + ", mass == massdiff. This is wrong. "
       "Please report it to the maintainer of the tool that wrote the pepXML. OpenMS will try to calculate it manually, assuming massdiff is correct.");
      if (term_spec_ == ResidueModification::N_TERM || term_spec_ == ResidueModification::PROTEIN_N_TERM)
      {
        mass_ = massdiff_ + Residue::getInternalToNTerm().getMonoWeight();
      }
      else if (term_spec_ == ResidueModification::C_TERM || term_spec_ == ResidueModification::PROTEIN_C_TERM)
      {
        mass_ = massdiff_ + Residue::getInternalToCTerm().getMonoWeight();
      }
      else
      {
        mass_ = massdiff_ + ResidueDB::getInstance()->getResidue(aminoacid_)->getMonoWeight(Residue::ResidueType::Internal);
      }
    }

    if (isVariable())
    {
      registered_mod_ = lookupModInPreferredMods_(preferred_var_mods, aminoacid_, massdiff_,
                                                  description_, term_spec_, tolerance);
    }
    else
    {
      registered_mod_ = lookupModInPreferredMods_(preferred_fixed_mods, aminoacid_, massdiff_,
                                                  description_, term_spec_, tolerance);
    }
    //TODO push another warning if not found?

    if (registered_mod_ == nullptr)
    {
      // check if the modification is uniquely defined through its description (if given):
      if (!description.empty())
      {
        try
        {
          registered_mod_ = ModificationsDB::getInstance()->getModification(description, aminoacid, term_spec_);
        }
        catch (Exception::BaseException&)
        {
          errors_.emplace_back("Modification '" + description_ + "' of residue '" + aminoacid_ +
                               "' could not be matched. Trying by modification mass.");
        }
      }
      else
      {
        errors_.emplace_back("No modification description given. Trying to define by modification mass.");
      }
    }

    if (registered_mod_ == nullptr)
    {
      //TODO differentiate between integer and non-integer masses like in AASequence?
      // Are there non-integer masses in the header of PepXML?
      std::vector<const ResidueModification*> mods;
      // if terminus was not specified
      if (term_spec_ == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) // try least specific search first
      {
        ModificationsDB::getInstance()->searchModificationsByDiffMonoMassSorted(
            mods, massdiff_, mod_tol_, aminoacid_, ResidueModification::ANYWHERE);
      }
      if (mods.empty())
      {
        // for unknown terminus it looks for everything, otherwise just for the specific terminus
        // TODO we might also need to search for Protein-X-Term in case of X-Term
        //  since some tools seem to forget to annotate
        ModificationsDB::getInstance()->searchModificationsByDiffMonoMassSorted(
            mods, massdiff_, mod_tol_, aminoacid_, term_spec_);
      }
      if (!mods.empty())
      {
        registered_mod_ = mods[0];
        if (mods.size() > 1)
        {
          String mod_str = mods[0]->getFullId();
          for (const auto& m : mods)
          {
            mod_str += ", " + m->getFullId();
          }
          errors_.emplace_back("Modification '" + String(mass_) + "' is not uniquely defined by the given data. Using '" +
                               mods[0]->getFullId() + "' to represent any of '" + mod_str + "'.");
        }
      }
      // If we could not find a registered mod in our DB, create and register it. This will be used later for lookup in the sequences.
      else if (massdiff_ != 0)
      {
        // r will be nullptr if not found. The next line handles it.
        const Residue* r = ResidueDB::getInstance()->getResidue(aminoacid_[0]);
        //TODO check if it is better to create from mass or massdiff
        registered_mod_ = ResidueModification::createUnknownFromMassString(String(massdiff_),
                                                                                 massdiff_,
                                                                                 true,
                                                                                 term_spec_,
                                                                                 r);

        //Modification unknown, but trying to continue as we want to be able to read the rest despite
        // of the modifications but warning this will fail downstream
        errors_.emplace_back(
            "Modification '" + String(mass_) + "/delta " + String(massdiff_) +
            "' is unknown. Resuming with '" + registered_mod_->getFullId() +
            "', which could lead to failures using the data downstream.");
      }
    }
  }

  PepXMLFile::PepXMLFile() :
    XMLHandler("", "1.12"),
    XMLFile("/SCHEMAS/pepXML_v114.xsd", "1.14"),
    proteins_(nullptr),
    peptides_(nullptr),
    lookup_(nullptr),
    scan_map_(),
    analysis_summary_(false),
    keep_native_name_(false),
    search_score_summary_(false),
    preferred_fixed_modifications_({}),
    preferred_variable_modifications_({})
  {
    const ElementDB* db = ElementDB::getInstance();
    hydrogen_ = *db->getElement("Hydrogen");
  }

  const double PepXMLFile::mod_tol_ = 0.002;
  const double PepXMLFile::xtandem_artificial_mod_tol_ = 0.0005; // according to cpp in some old version of xtandem somehow very small fixed modification (electron mass?) gets annotated by X!Tandem. Don't add them as they interfere with other modifications.

  PepXMLFile::~PepXMLFile() = default;

  void PepXMLFile::store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, const String& mz_file, const String& mz_name, bool peptideprophet_analyzed, double rt_tolerance)
  {
    ofstream f(filename.c_str());
    if (!f)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    String search_engine_name;
    ProteinIdentification::SearchParameters search_params;
    if (!protein_ids.empty())
    {
      if (protein_ids.size() > 1)
      {
        warning(STORE, "More than one protein identification defined; only first search parameters are written into pepXML; search engine must be the same.");
      }
      search_params = protein_ids.begin()->getSearchParameters();
      if (protein_ids.begin()->getSearchEngine() == "XTandem")
      {
        search_engine_name = "X! Tandem";
      }
      else if (protein_ids.begin()->getSearchEngine() == "Mascot")
      {
        search_engine_name = "MASCOT";
      }
      else
      {
        search_engine_name = protein_ids.begin()->getSearchEngine();
        //Comet writes "Comet" in pep.xml, so this is ok
      }
    }

    f.precision(writtenDigits<double>(0.0));
    String raw_data;
    String base_name;
    SpectrumMetaDataLookup lookup;
    lookup.rt_tolerance = rt_tolerance;

    // The mz-File (if given)
    if (!mz_file.empty())
    {
      base_name = FileHandler::stripExtension(File::basename(mz_file));
      raw_data = FileTypes::typeToName(FileHandler::getTypeByFileName(mz_file));

      PeakMap experiment;
      FileHandler fh;
      fh.loadExperiment(mz_file, experiment, {}, ProgressLogger::NONE, false, false);
      lookup.readSpectra(experiment.getSpectra());
    }
    else
    {
      base_name = FileHandler::stripExtension(File::basename(filename));
      raw_data = "mzML";
    }
    // mz_name is input from IDFileConverter for 'base_name' attribute, only necessary if different from 'mz_file'.
    if (!mz_name.empty())
    {
      base_name = mz_name;
    }
    if (base_name.hasSubstring(".")) // spectrum query name is split by dot, otherwise correct charge can not be read.
    {
      replace(base_name.begin(), base_name.end(), '.', '_');
    }

    f << R"(<?xml version="1.0" encoding="UTF-8"?>)" << "\n";
    f << R"(<msms_pipeline_analysis date="2007-12-05T17:49:46" xmlns="http://regis-web.systemsbiology.net/pepXML" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd" summary_xml=".xml">)" << "\n";
    f << "<msms_run_summary base_name=\"" << base_name << R"(" raw_data_type="raw" raw_data=".)" << raw_data << "\" search_engine=\"" << search_engine_name << "\">" << "\n";
    String enzyme_name = search_params.digestion_enzyme.getName();
    f << "\t<sample_enzyme name=\"";
    f << enzyme_name.toLower() << "\">" << "\n";
    f << "\t\t<specificity cut=\"";
    if (!search_params.digestion_enzyme.getRegEx().empty())
    {
      vector<String> sub_regex;
      search_params.digestion_enzyme.getRegEx().split(")",sub_regex);
      boost::match_results<std::string::const_iterator> results;
      static const boost::regex e("(.*?)([A-Z]+)(.*?)");
      if (boost::regex_match(sub_regex[0], results, e))
      {
        f << results[2];
      }
      if (sub_regex[1].hasSubstring("!P"))
      {
        f << "\" no_cut=\"P";
      }
    }
    f << R"(" sense="C"/>)" << "\n";
    f << "\t</sample_enzyme>" << "\n";

    f << "\t<search_summary base_name=\"" << base_name;
    f << "\" search_engine=\"" << search_engine_name;
    f << "\" precursor_mass_type=\"";
    if (search_params.mass_type == ProteinIdentification::MONOISOTOPIC)
    {
      f << "monoisotopic";
    }
    else
    {
      f << "average";
    }
    f << "\" fragment_mass_type=\"";
    if (search_params.mass_type == ProteinIdentification::MONOISOTOPIC)
    {
      f << "monoisotopic";
    }
    else
    {
      f << "average";
    }
    f << R"(" out_data_type="" out_data="" search_id="1">)" << "\n";
    f << "\t\t<search_database local_path=\"" << search_params.db << R"(" type="AA"/>)" << "\n";


    // register modifications
    set<String> aa_mods;
    set<String> n_term_mods, c_term_mods;
    for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin();
         it != peptide_ids.end(); ++it)
    {
      if (!it->getHits().empty())
      {
        PeptideHit h = *it->getHits().begin();

        if (h.getSequence().isModified())
        {
          const AASequence& p = h.getSequence();
          if (p.hasNTerminalModification())
          {
            n_term_mods.insert(p.getNTerminalModification()->getFullId());
          }
          if (p.hasCTerminalModification())
          {
            c_term_mods.insert(p.getCTerminalModification()->getFullId());
          }

          for (Size i = 0; i != p.size(); ++i)
          {
            if (p[i].isModified())
            {
              aa_mods.insert(p[i].getModification()->getFullId());
            }
          }
        }
      }
    }

    // write modifications definitions
    // <aminoacid_modification aminoacid="C" massdiff="+58.01" mass="161.014664" variable="Y" binary="N" description="Carboxymethyl (C)"/>
    for (set<String>::const_iterator it = aa_mods.begin();
         it != aa_mods.end(); ++it)
    {
      const ResidueModification* mod = ModificationsDB::getInstance()->getModification(*it, "", ResidueModification::ANYWHERE);

      // compute mass of modified residue
      EmpiricalFormula ef = ResidueDB::getInstance()->getResidue(mod->getOrigin())->getFormula(Residue::Internal);
      ef += mod->getDiffFormula();

      f << "\t\t"
        << "<aminoacid_modification aminoacid=\"" << mod->getOrigin()
        << "\" massdiff=\"" << precisionWrapper(mod->getDiffMonoMass()) << "\" mass=\""
        << precisionWrapper(ef.getMonoWeight())
        << R"(" variable="Y" binary="N" description=")" << *it << "\"/>"
        << "\n";
    }

    for (set<String>::const_iterator it = n_term_mods.begin(); it != n_term_mods.end(); ++it)
    {
      const ResidueModification* mod = ModificationsDB::getInstance()->getModification(*it, "", ResidueModification::N_TERM);
      f << "\t\t"
        << R"(<terminal_modification terminus="n" massdiff=")"
        << precisionWrapper(mod->getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod->getMonoMass())
        << R"(" variable="Y" description=")" << *it
        << R"(" protein_terminus=""/>)" << "\n";
    }

    for (set<String>::const_iterator it = c_term_mods.begin(); it != c_term_mods.end(); ++it)
    {
      const ResidueModification* mod = ModificationsDB::getInstance()->getModification(*it, "", ResidueModification::C_TERM);
      f << "\t\t"
        << R"(<terminal_modification terminus="c" massdiff=")"
        << precisionWrapper(mod->getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod->getMonoMass())
        << R"(" variable="Y" description=")" << *it
        << R"(" protein_terminus=""/>)" << "\n";
    }

    f << "\t</search_summary>" << "\n";
    if (peptideprophet_analyzed)
    {
      f << "\t<analysis_timestamp analysis=\"peptideprophet\" time=\"2007-12-05T17:49:52\" id=\"1\"/>" << "\n";
    }

    // Scan index and scan number will be reconstructed if no spectrum lookup is possible to retrieve the values.
    // The scan index is generally zero-based and the scan number generally one-based.
    Int count(0);
    for (const PeptideIdentification& pep : peptide_ids)
    {
      if (pep.getHits().empty())
      {
        ++count;
        continue;
      }
      for (const PeptideHit& hit : pep.getHits())
      {
        PeptideHit h = hit;
        const AASequence& seq = h.getSequence();
        double precursor_neutral_mass = seq.getMonoWeight();

        int scan_index = count;
        int scan_nr = 0;

        if (lookup.empty())
        {
          if (pep.metaValueExists("RT_index")) // Setting metaValue "RT_index" in XTandemXMLFile in the case of X! Tandem.
          {
            scan_index = pep.getMetaValue("RT_index");
          }
          scan_nr = scan_index + 1;
        }
        else
        {
          if (pep.metaValueExists("spectrum_reference"))
          {
            //findByNativeID will fall back to RT lookup if none of the regexes registered in lookup can extract a meaningful ID or scan nr
            scan_index = lookup.findByNativeID(pep.getSpectrumReference());
          }
          else
          {
            scan_index = lookup.findByRT(pep.getRT());
          }

          SpectrumMetaDataLookup::SpectrumMetaData meta;
          lookup.getSpectrumMetaData(scan_index, meta);
          scan_nr = meta.scan_number;
        }
        // PeptideProphet requires this format for "spectrum" attribute (otherwise TPP parsing error)
        //  - see also the parser code if iProphet at http://sourceforge.net/p/sashimi/code/HEAD/tree/trunk/trans_proteomic_pipeline/src/Validation/InterProphet/InterProphetParser/InterProphetParser.cxx#l180
        //  strictly required attributes:
        //    - spectrum
        //    - assumed_charge
        //  optional attributes
        //    - retention_time_sec
        //    - swath_assay
        //    - experiment_label

        String spectrum_name = base_name + "." + scan_nr + "." + scan_nr + ".";
        if (pep.metaValueExists("pepxml_spectrum_name") && keep_native_name_)
        {
          spectrum_name = pep.getMetaValue("pepxml_spectrum_name");
        }

        f << "\t<spectrum_query spectrum=\"" << spectrum_name << h.getCharge() << "\""
          << " start_scan=\"" << scan_nr << "\""
          << " end_scan=\"" << scan_nr << "\""
          << " precursor_neutral_mass=\"" << precisionWrapper(precursor_neutral_mass) << "\""
          << " assumed_charge=\"" << h.getCharge() << "\" index=\"" << scan_index << "\"";

        if (pep.hasRT())
        {
          f << " retention_time_sec=\"" << pep.getRT() << "\" ";
        }

        if (!pep.getExperimentLabel().empty())
        {
          f << " experiment_label=\"" << pep.getExperimentLabel() << "\" ";
        }

        // "swath_assay" is an optional parameter used for SWATH-MS mostly and
        // may be set for a PeptideIdentification
        //   note that according to the parsing rules of TPP, this needs to be
        //   "xxx:yyy" where xxx is any string and yyy is probably an integer
        //   indicating the Swath window
        if (pep.metaValueExists("swath_assay"))
        {
          f << " swath_assay=\"" << pep.getMetaValue("swath_assay") << "\" ";
        }
        // "status" is an attribute that may be target or decoy
        if (pep.metaValueExists("status"))
        {
          f << " status=\"" << pep.getMetaValue("status") << "\" ";
        }

        f << ">\n";
        f << "\t<search_result>" << "\n";

        vector<PeptideEvidence> pes = h.getPeptideEvidences();

        // select first one if multiple are present as "leader"
        PeptideEvidence pe;
        if (!h.getPeptideEvidences().empty())
        {
          pe = pes[0];
        }

        f << "\t\t<search_hit hit_rank=\"1\" peptide=\""
          << seq.toUnmodifiedString() << "\" peptide_prev_aa=\""
          << pe.getAABefore() << "\" peptide_next_aa=\"" << pe.getAAAfter()
          << "\" protein=\"";

        f << pe.getProteinAccession();

        f << R"(" num_tot_proteins="1" num_matched_ions="0" tot_num_ions="0" calc_neutral_pep_mass=")" << precisionWrapper(precursor_neutral_mass)
          << R"(" massdiff="0.0" num_tol_term=")";
        Int num_tol_term = 1;
        if ((pe.getAABefore() == 'R' || pe.getAABefore() == 'K') && search_params.digestion_enzyme.getName() == "Trypsin")
        {
          num_tol_term = 2;
        }
        f << num_tol_term;
        f << R"(" num_missed_cleavages="0" is_rejected="0" protein_descr="Protein No. 1">)" << "\n";

        // multiple protein hits: <alternative_protein protein="sp|P0CZ86|GLS24_STRP3" num_tol_term="2" peptide_prev_aa="K" peptide_next_aa="-"/>
        if (pes.size() > 1)
        {
          for (Size k = 1; k != pes.size(); ++k)
          {
            f << "\t\t<alternative_protein protein=\"" << pes[k].getProteinAccession() << "\" num_tol_term=\"" << num_tol_term << "\"";

            if (pes[k].getAABefore() != PeptideEvidence::UNKNOWN_AA)
            {
              f << " peptide_prev_aa=\"" << pes[k].getAABefore() << "\"";
            }

            if (pes[k].getAAAfter() != PeptideEvidence::UNKNOWN_AA)
            {
              f << " peptide_next_aa=\"" << pes[k].getAAAfter() << "\"";
            }
            f << "/>" << "\n";
          }
        }
        if (seq.isModified())
        {
          f << "\t\t\t<modification_info modified_peptide=\""
            << seq.toBracketString() << "\"";

          if (seq.hasNTerminalModification())
          {
            const ResidueModification* mod = seq.getNTerminalModification();
            const double mod_nterm_mass = Residue::getInternalToNTerm().getMonoWeight() + mod->getDiffMonoMass();
            f << " mod_nterm_mass=\"" << precisionWrapper(mod_nterm_mass) << "\"";
          }

          if (seq.hasCTerminalModification())
          {
            const ResidueModification* mod = seq.getCTerminalModification();
            const double mod_cterm_mass = Residue::getInternalToCTerm().getMonoWeight() + mod->getDiffMonoMass();
            f << " mod_cterm_mass=\"" << precisionWrapper(mod_cterm_mass) << "\"";
          }

          f << ">" << "\n";

          for (Size i = 0; i != seq.size(); ++i)
          {
            if (seq[i].isModified())
            {
              const ResidueModification* mod = seq[i].getModification();
              // the modification position is 1-based
              f << "\t\t\t\t<mod_aminoacid_mass position=\"" << (i + 1)
                << "\" mass=\"" <<
                precisionWrapper(mod->getMonoMass() + seq[i].getMonoWeight(Residue::Internal)) << "\"/>" << "\n";
            }
          }

          f << "\t\t\t</modification_info>" << "\n";
        }

        // write out the (optional) search_score_summary that may be associated with peptide prophet results
        bool peptideprophet_written = false;
        if (!h.getAnalysisResults().empty())
        {
          // <analysis_result analysis="peptideprophet">
          //   <peptideprophet_result probability="0.0660" all_ntt_prob="(0.0000,0.0000,0.0660)">
          //     <search_score_summary>
          //       <parameter name="fval" value="0.7114"/>
          //       <parameter name="ntt" value="2"/>
          //       <parameter name="nmc" value="0"/>
          //       <parameter name="massd" value="-0.027"/>
          //       <parameter name="isomassd" value="0"/>
          //     </search_score_summary>
          //   </peptideprophet_result>
          // </analysis_result>

          for (const PeptideHit::PepXMLAnalysisResult& ar_it : h.getAnalysisResults())
          {
            f << "\t\t\t<analysis_result analysis=\"" << ar_it.score_type << "\">" << "\n";

            // get name of next tag
            String tagname = "peptideprophet_result";
            if (ar_it.score_type == "peptideprophet")
            {
              peptideprophet_written = true; // remember that we have now already written peptide prophet results
              tagname = "peptideprophet_result";
            }
            else if (ar_it.score_type == "interprophet")
            {
              tagname = "interprophet_result";
            }
            else
            {
              peptideprophet_written = true; // remember that we have now already written peptide prophet results
              warning(STORE, "Analysis type " + ar_it.score_type + " not supported, will use peptideprophet_result.");
            }

            f << "\t\t\t\t<" << tagname <<  " probability=\"" << ar_it.main_score;
            // TODO
            f << "\" all_ntt_prob=\"(" << ar_it.main_score << "," << ar_it.main_score
            << "," << ar_it.main_score << ")\">" << "\n";

            if (!ar_it.sub_scores.empty())
            {
              f << "\t\t\t\t\t<search_score_summary>" << "\n";
              for (std::map<String, double>::const_iterator subscore_it = ar_it.sub_scores.begin();
                  subscore_it != ar_it.sub_scores.end(); ++subscore_it)
              {
                f << "\t\t\t\t\t\t<parameter name=\""<< subscore_it->first << "\" value=\"" << subscore_it->second << "\"/>\n";
              }
              f << "\t\t\t\t\t</search_score_summary>" << "\n";
            }
            f << "\t\t\t\t</" << tagname << ">" << "\n";

            f << "\t\t\t</analysis_result>" << "\n";
          }
        }

        // deprecated way of writing out peptide prophet results (only if
        // requested explicitly and if we have not already written out the
        // peptide prophet results above through AnalysisResults
        if (peptideprophet_analyzed && !peptideprophet_written)
        {
          // if (!h.getAnalysisResults().empty()) { WARNING / }
          f << "\t\t\t<analysis_result analysis=\"peptideprophet\">" << "\n";
          f << "\t\t\t<peptideprophet_result probability=\"" << h.getScore()
            << "\" all_ntt_prob=\"(" << h.getScore() << "," << h.getScore()
            << "," << h.getScore() << ")\">" << "\n";

          f << "\t\t\t</peptideprophet_result>" << "\n";
          f << "\t\t\t</analysis_result>" << "\n";
        }
        else
        {
          bool haspep = pep.getScoreType() == "Posterior Error Probability" || pep.getScoreType() == "pep";
          bool percolator = false;
          if (search_engine_name == "X! Tandem")
          {
            // check if score type is XTandem or qvalue/fdr
            if (pep.getScoreType() == "XTandem")
            {
              f << "\t\t\t<search_score" << R"( name="hyperscore" value=")" << h.getScore() << "\"" << "/>\n";
              f << "\t\t\t<search_score" << R"( name="nextscore" value=")";
              if (h.metaValueExists("nextscore"))
              {
                f << h.getMetaValue("nextscore") << "\"" << "/>\n";
              }
              else
              {
                f << h.getScore() << "\"" << "/>\n";
              }
            }
            else if (h.metaValueExists("XTandem_score"))
            {
              f << "\t\t\t<search_score" << R"( name="hyperscore" value=")" << h.getMetaValue("XTandem_score") << "\"" << "/>\n";
              f << "\t\t\t<search_score" << R"( name="nextscore" value=")";
              if (h.metaValueExists("nextscore"))
              {
                f << h.getMetaValue("nextscore") << "\"" << "/>\n";
              }
              else
              {
                f << h.getMetaValue("XTandem_score") << "\"" << "/>\n";
              }
            }
            f << "\t\t\t<search_score" << R"( name="expect" value=")" << h.getMetaValue("E-Value") << "\"" << "/>\n";
          }
          else if (search_engine_name == "Comet")
          {
            f << "\t\t\t<search_score" << R"( name="xcorr" value=")" << h.getMetaValue("MS:1002252") << "\"" << "/>\n"; // name: Comet:xcorr
            f << "\t\t\t<search_score" << R"( name="deltacn" value=")" << h.getMetaValue("MS:1002253") << "\"" << "/>\n"; // name: Comet:deltacn
            f << "\t\t\t<search_score" << R"( name="deltacnstar" value=")" << h.getMetaValue("MS:1002254") << "\"" << "/>\n"; // name: Comet:deltacnstar
            f << "\t\t\t<search_score" << R"( name="spscore" value=")" << h.getMetaValue("MS:1002255") << "\"" << "/>\n"; // name: Comet:spscore
            f << "\t\t\t<search_score" << R"( name="sprank" value=")" << h.getMetaValue("MS:1002256") << "\"" << "/>\n"; // name: Comet:sprank
            f << "\t\t\t<search_score" << R"( name="expect" value=")" << h.getMetaValue("MS:1002257") << "\"" << "/>\n"; // name: Comet:expect
          }
          else if (search_engine_name == "MASCOT")
          {
            f << "\t\t\t<search_score" << R"( name="expect" value=")" << h.getMetaValue("EValue") << "\"" << "/>\n";
            f << "\t\t\t<search_score" << R"( name="ionscore" value=")" << h.getScore() << "\"" << "/>\n";
          }
          else if (search_engine_name == "OMSSA")
          {
            f << "\t\t\t<search_score" << R"( name="expect" value=")" << h.getScore() << "\"" << "/>\n";
          }
          else if (search_engine_name == "MSGFPlus")
          {
            f << "\t\t\t<search_score" << R"( name="expect" value=")" << h.getScore() << "\"" << "/>\n";
          }
          else if (search_engine_name == "Percolator")
          {
            double svm_score = 0.0;
            if (h.metaValueExists("MS:1001492"))
            {
              svm_score = static_cast<double>(h.getMetaValue("MS:1001492"));
              f << "\t\t\t<search_score" << R"( name="Percolator_score" value=")" << svm_score << "\"" << "/>\n";
            }
            else if (h.metaValueExists("Percolator_score"))
            {
              svm_score = static_cast<double>(h.getMetaValue("Percolator_score"));
              f << "\t\t\t<search_score" << R"( name="Percolator_score" value=")" << svm_score << "\"" << "/>\n";
            }

            double qval_score = 0.0;
            if (h.metaValueExists("MS:1001491"))
            {
              qval_score = static_cast<double>(h.getMetaValue("MS:1001491"));
              f << "\t\t\t<search_score" << R"( name="Percolator_qvalue" value=")" << qval_score << "\"" << "/>\n";
            }
            else if (h.metaValueExists("Percolator_qvalue"))
            {
              qval_score = static_cast<double>(h.getMetaValue("Percolator_qvalue"));
              f << "\t\t\t<search_score" << R"( name="Percolator_qvalue" value=")" << qval_score << "\"" << "/>\n";
            }

            double pep_score = 0.0;
            if (h.metaValueExists("MS:1001493"))
            {
              pep_score = static_cast<double>(h.getMetaValue("MS:1001493"));
            }
            else if (h.metaValueExists("Percolator_PEP"))
            {
              pep_score = static_cast<double>(h.getMetaValue("Percolator_PEP"));
            }
            else
            {
              if (!haspep)
              {
                // will be written later
                throw Exception::MissingInformation(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION,"Percolator PEP score missing for pepXML export of Percolator results.");
              }
            }
            f << "\t\t\t<search_score" << R"( name="Percolator_PEP" value=")" << pep_score << "\"" << "/>\n";

            f << "\t\t\t<analysis_result" << " analysis=\"peptideprophet\">\n";
            f << "\t\t\t\t<peptideprophet_result" << " probability=\"" << 1. - pep_score << "\"";
            f << " all_ntt_prob=\"(0.0000,0.0000," << 1. - pep_score << ")\"/>\n";
            f << "\t\t\t</analysis_result>" << "\n";
            percolator = true;
          } // Anything else
          else
          {
            f << "\t\t\t<search_score" << " name=\"" << pep.getScoreType() << "\" value=\"" << h.getScore() << "\"" << "/>\n";
          }
          // Any search engine with a PEP (e.g. also our IDPEP) except Percolator which has
          // written that part already
          if (haspep && !percolator)
          {
            f << "\t\t\t<search_score" << " name=\"" << pep.getScoreType() << "\" value=\"" << h.getScore() << "\"" << "/>\n";
            double probability = 1.0 - h.getScore();
            f << "\t\t\t<analysis_result" << " analysis=\"peptideprophet\">\n";
            f << "\t\t\t\t<peptideprophet_result" << " probability=\"" << probability << "\"";
            f << " all_ntt_prob=\"(0.0000,0.0000," << probability << ")\"/>\n";
            f << "\t\t\t</analysis_result>" << "\n";
          }
        }
        f << "\t\t</search_hit>" << "\n";
        f << "\t</search_result>" << "\n";
        f << "\t</spectrum_query>" << "\n";
      }
      ++count;
    }
    f << "</msms_run_summary>" << "\n";
    f << "</msms_pipeline_analysis>" << "\n";

    f.close();
  }

  void PepXMLFile::readRTMZCharge_(const xercesc::Attributes& attributes)
  {
    double mass = attributeAsDouble_(attributes, "precursor_neutral_mass");
    charge_ = attributeAsInt_(attributes, "assumed_charge");
    mz_ = (mass + hydrogen_mass_ * charge_) / charge_;
    rt_ = 0;
    // assume only one scan, i.e. ignore "end_scan":
    // "start and end_scan" are 1-based. "index" is 0-based
    scannr_ = attributeAsInt_(attributes, "start_scan");
    Size endscan = attributeAsInt_(attributes, "start_scan");
    if (scannr_ != endscan)
    {
      error(LOAD, "endscan not equal to startscan. Merged spectrum queries not supported. Parsing start scan nr. only.");
    }

    bool rt_present = optionalAttributeAsDouble_(rt_, attributes, "retention_time_sec");

    if (!rt_present) // get RT from experiment
    {
      if (lookup_ == nullptr || lookup_->empty())
      {
        // no lookup given, report non-fatal error
        error(LOAD, "Cannot get RT information - no spectra given");
        return;
      }

      Size index = (scannr_ != 0) ? lookup_->findByScanNumber(scannr_) :
        lookup_->findByReference(attributeAsString_(attributes, "spectrum"));
      SpectrumMetaDataLookup::SpectrumMetaData meta;
      lookup_->getSpectrumMetaData(index, meta);
      if (meta.ms_level == 2)
      {
        rt_ = meta.rt;
      }
      else
      {
        error(LOAD, "Cannot get RT information - scan mapping is incorrect");
      }
    }
  }

  void PepXMLFile::setParseUnknownScores(bool parse_unknown_scores)
  {
    this->parse_unknown_scores_ = parse_unknown_scores;
  }

  void PepXMLFile::load(const String& filename, vector<ProteinIdentification>&
                        proteins, vector<PeptideIdentification>& peptides,
                        const String& experiment_name
                        )
  {
    SpectrumMetaDataLookup lookup;
    load(filename, proteins, peptides, experiment_name, lookup);
  }

  void PepXMLFile::load(const String& filename, vector<ProteinIdentification>&
                        proteins, vector<PeptideIdentification>& peptides,
                        const String& experiment_name,
                        const SpectrumMetaDataLookup& lookup)
  {
    // initialize here, since "load" could be called several times:
    exp_name_ = "";
    prot_id_ = "";
    charge_ = 0;
    peptides.clear();
    peptides_ = &peptides;
    proteins.clear();
    proteins_ = &proteins;
    // assume mass type "average" (in case element "search_summary" is missing):
    hydrogen_mass_ = hydrogen_.getAverageWeight();

    file_ = filename; // filename for error messages in XMLHandler

    if (!experiment_name.empty())
    {
      exp_name_ = FileHandler::stripExtension(experiment_name);
      lookup_ = &lookup;
    }

    analysis_summary_ = false;
    wrong_experiment_ = false;
    // without experiment name, don't care about these two:
    seen_experiment_ = exp_name_.empty();
    checked_base_name_ = exp_name_.empty();

    parse_(filename, this);

    if (!seen_experiment_)
    {
      fatalError(LOAD, "Found no experiment with name '" + experiment_name + "'");
    }
    // clean up duplicate ProteinHits in each ProteinIdentification separately:
    // (can't use "sort" and "unique" because no "op<" defined for ProteinHit)
    for (ProteinIdentification& prot_it : proteins)
    {
      set<String> accessions;
      // modeled after "remove_if" in STL header "algorithm":
      vector<ProteinHit>::iterator first = prot_it.getHits().begin();
      vector<ProteinHit>::iterator result = first;
      for (; first != prot_it.getHits().end(); ++first)
      {
        String accession = first->getAccession();
        bool new_element = accessions.insert(accession).second;
        if (new_element) // don't remove
        {
          *result++ = *first;
        }
      }
      prot_it.getHits().erase(result, first);
    }

    // reset members
    exp_name_.clear();
    prot_id_.clear();
    date_.clear();
    proteins_ = nullptr;
    peptides_ = nullptr;
    lookup_ = nullptr;
    scan_map_.clear();
  }

  /*
    NOTE: numbering schemes for multiple searches
    ---------------------------------------------

    pepXML can contain search results for multiple files and for multiple
    searches of those files. We have seen two different variants (both seem to
    be valid pepXML):

    1. One "msms_run_summary" per searched file, containing multiple
    "search_summary" elements (for each search); counting starts at "1" in each
    "msms_run_summary"; produced e.g. by the TPP:

    <msms_run_summary basename="File1">
      <search_summary search_engine="A" search_id="1">
      <search_summary search_engine="B" search_id="2">
    ...
    <msms_run_summary basename="File2">
      <search_summary search_engine="A" search_id="1">
      <search_summary search_engine="B" search_id="2">
    ...

    2. One "msms_run_summary" per search of a file, containing one
    "search_summary" element; counting is sequential across "msms_run_summary"
    sections; produced e.g. by ProteomeDiscoverer:

    <msms_run_summary basename="File1">
      <search_summary search_engine="A" search_id="1">
    ...
    <msms_run_summary basename="File1">
      <search_summary search_engine="B" search_id="2">
    ...

    The "search_id" numbers are necessary to associate search hits with the
    correct search runs. In the parser, we keep track of the search runs per
    "msms_run_summary" in the "current_proteins_" vector. Importantly, this
    means that for variant 2 of the numbering, the "search_id" number may be
    higher than the number of elements in the vector! However, in this case we
    know that the last - and only - element should be the correct one.
  */

  void PepXMLFile::startElement(const XMLCh* const /*uri*/,
                                const XMLCh* const /*local_name*/,
                                const XMLCh* const qname,
                                const xercesc::Attributes& attributes)
  {
    String element = sm_.convert(qname);

    // cout << "Start: " << element << "\n";

    if (element == "msms_run_summary") // parent: "msms_pipeline_analysis"
    {
      String ms_run_path;
      if (!exp_name_.empty())
      {
        String base_name = attributeAsString_(attributes, "base_name");
        if (!base_name.empty())
        {
          wrong_experiment_ = !base_name.hasSuffix(exp_name_);
          seen_experiment_ = seen_experiment_ || !wrong_experiment_;
          checked_base_name_ = true;
        }
        else // really shouldn't happen, but does for Mascot export to pepXML
        {
          error(LOAD, "'base_name' attribute of 'msms_run_summary' element is empty");
          // continue for now, but check later in 'search_summary':
          wrong_experiment_ = false;
          checked_base_name_ = false;
        }
        String raw_data = attributeAsString_(attributes, "raw_data");
        if (!base_name.empty() && !raw_data.empty())
        {
          ms_run_path = base_name + "." + raw_data;
        }
      }
      if (wrong_experiment_) return;

      // create a ProteinIdentification in case "search_summary" is missing:
      ProteinIdentification protein;
      protein.setDateTime(date_);
      prot_id_ = "unknown_" + date_.getDate();
      enzyme_ = "unknown_enzyme";
      // "prot_id_" will be overwritten if elem. "search_summary" is present
      protein.setIdentifier(prot_id_);
      proteins_->push_back(protein);
      if (!ms_run_path.empty())
      {
        protein.setPrimaryMSRunPath(StringList(1, ms_run_path));
      }
      current_proteins_.clear();
      current_proteins_.push_back(--proteins_->end());
    }
    else if (element == "analysis_summary") // parent: "msms_pipeline_analysis"
    {
      // this element can contain "search summary" elements, which we only
      // expect as subelements of "msms_run_summary", so skip the whole thing
      analysis_summary_ = true;
    }
    else if (wrong_experiment_ || analysis_summary_)
    {
      // do nothing here (this case exists to prevent parsing of elements for
      // experiments we're not interested in or for analysis summaries)
    }
    // now, elements occurring more frequently are generally closer to the top
    else if (element == "search_score") // parent: "search_hit"
    {
      String name = attributeAsString_(attributes, "name");
      double value;

      // TODO: deal with different scores
      if (name == "expect") // X!Tandem, Mascot or MSFragger E-value
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType(name);
        current_peptide_.setHigherScoreBetter(false);
        if (search_engine_ == "Comet")
        {
          peptide_hit_.setMetaValue("MS:1002257", value); // name: Comet:expectation value
        }
        else if (search_engine_ == "X! Tandem" || search_engine_ == "MSFragger") // TODO check if there are separate CVs?
        {
          peptide_hit_.setMetaValue("MS:1001330", value); // name: X\!Tandem:expect
        }
        else if (search_engine_ == "Mascot")
        {
          peptide_hit_.setMetaValue("MS:1001172", value); // name: Mascot:expectation value
        }
        //TODO: there is no (generic) umbrella term for expect val in the CV right now
      }
      else if (name == "mvh") // MyriMatch score
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType(name);
        current_peptide_.setHigherScoreBetter(true);
      }
      // if (name == "hyperscore")
      // { // X!Tandem score
      //   value = attributeAsDouble_(attributes, "value");
      //   peptide_hit_.setScore(value);
      //   current_peptide_.setScoreType(name); // add "X!Tandem" to name?
      //   current_peptide_.setHigherScoreBetter(true);
      // }
      else if (name == "xcorr") // Sequest score
      { // and use the mvh
        value = attributeAsDouble_(attributes, "value");
        if (search_engine_ != "MyriMatch") //MyriMatch has also an xcorr, but we want to ignore it
        {
          peptide_hit_.setScore(value);
          current_peptide_.setScoreType(name); // add "Sequest" to name?
          current_peptide_.setHigherScoreBetter(true);
        }
        if (search_engine_ == "Comet")
        {
          peptide_hit_.setMetaValue("MS:1002252", value); // name: Comet:xcorr
          peptide_hit_.setMetaValue("COMET:xcorr", value); // name: COMET:xcorr
        }
        else
        {
          peptide_hit_.setMetaValue("MS:1001155", value); // name: SEQUEST:xcorr
        }
        //TODO: no other xcorr or generic xcorr in the CV right now, use SEQUEST:xcorr
      }
      else if (name == "fval") // SpectraST score
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType(name);
        current_peptide_.setHigherScoreBetter(true);
        peptide_hit_.setMetaValue("MS:1001419", value); // def: "SpectraST spectrum score.
      }
      else if (name == "hyperscore")
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setMetaValue("hyperscore", value);
      }
      else if (name == "nextscore")
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setMetaValue("nextscore", value);
      }
      else if (search_engine_ == "Comet")
      {
        if (name == "deltacn")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("MS:1002253", value); // name: Comet:deltacn
          peptide_hit_.setMetaValue("COMET:deltaCn", value);
        }
        else if (name == "spscore")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("MS:1002255", value); // name: Comet:spscore
          peptide_hit_.setMetaValue("COMET:spscore", value); // name: Comet:spscore
        }
        else if (name == "sprank")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("MS:1002256", value); // name: Comet:sprank
          peptide_hit_.setMetaValue("COMET:sprank", value); // name: Comet:sprank
        }
        else if (name == "deltacnstar")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("MS:1002254", value); // name: Comet:deltacnstar
          peptide_hit_.setMetaValue("COMET:deltacnstar", value); // name: Comet:deltacnstar
        }
        else if (name == "lnrSp")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("Comet:lnrSp", value); // name: Comet:lnrSp
          peptide_hit_.setMetaValue("COMET:lnRankSP", value); // name: COMET:lnRankSP
        }              
        else if (name == "deltLCn")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("COMET:deltaLCn", value); // name: Comet:deltLCn
        }
        else if (name == "lnExpect")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("COMET:lnExpect", value); // name: Comet:lnExpect          
        }
        else if (name == "IonFrac")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("Comet:IonFrac", value); // name: Comet:IonFrac
          peptide_hit_.setMetaValue("COMET:IonFrac", value); // matched_ions / total_ions
        }
        else if (name == "lnNumSP")
        {
          value = attributeAsDouble_(attributes, "value");
          peptide_hit_.setMetaValue("COMET:lnNumSP", value); // name: Comet:lnNumSP
        }        
      }
      else if (parse_unknown_scores_)
      {
        String strvval = attributeAsString_(attributes, "value");
        if(name.hasSuffix("_ions")) //TODO create a dictionary of which known scores are which type?
        {
          peptide_hit_.setMetaValue(name, attributeAsInt_(attributes, "value")); // e.g. name: Comet:matched_b1_ions
        }
        else
        {
          try
          {
            peptide_hit_.setMetaValue(name, attributeAsDouble_(attributes, "value")); // Any other generic score
          }
          catch (Exception::ConversionError&)
          {
            //TODO warn about non-numeric score? Or even do not catch the conversion error?
            peptide_hit_.setMetaValue(name, attributeAsString_(attributes, "value")); // Any other generic score (fallback String)
          }
          
        }
      }
    }
    else if (element == "search_hit") // parent: "search_result"
    { // creates a new PeptideHit
      current_sequence_ = attributeAsString_(attributes, "peptide");
      current_modifications_.clear();
      PeptideEvidence pe;
      peptide_hit_ = PeptideHit();
      peptide_hit_.setRank(attributeAsInt_(attributes, "hit_rank"));
      peptide_hit_.setCharge(charge_); // from parent "spectrum_query" tag
      String prev_aa, next_aa;
      if (optionalAttributeAsString_(prev_aa, attributes, "peptide_prev_aa"))
      {
        pe.setAABefore(prev_aa[0]);
      }

      if (optionalAttributeAsString_(next_aa, attributes, "peptide_next_aa"))
      {
        pe.setAAAfter(next_aa[0]);
      }
      if (search_engine_ == "Comet")
      {
        String value;
        if (optionalAttributeAsString_(value, attributes, "num_matched_ions"))
        {
          peptide_hit_.setMetaValue("MS:1002258", value); // name: Comet:matched ions

        }
        if (optionalAttributeAsString_(value, attributes, "tot_num_ions"))
        {
          peptide_hit_.setMetaValue("MS:1002259", value); // name: Comet:total ions
        }
        if (optionalAttributeAsString_(value, attributes, "num_matched_peptides"))
        {
          peptide_hit_.setMetaValue("num_matched_peptides", value);
        }
      }

      double pepmassdiff = attributeAsDouble_(attributes, "massdiff");
      peptide_hit_.setMetaValue(Constants::UserParam::ISOTOPE_ERROR, std::lrint(pepmassdiff/Constants::C13C12_MASSDIFF_U));

      String protein = attributeAsString_(attributes, "protein");
      protein.trim();
      pe.setProteinAccession(protein);

      ProteinHit hit;
      hit.setAccession(protein);

      if (has_decoys_)
      {
        String curr_status("");
        bool current_prot_is_decoy = protein.hasPrefix(decoy_prefix_);
        if (peptide_hit_.metaValueExists("target_decoy"))
        {
          curr_status = peptide_hit_.getMetaValue("target_decoy");
        }
        if (curr_status.empty())
        {
          peptide_hit_.setMetaValue("target_decoy", current_prot_is_decoy ? "decoy" : "target");
        }
        else if (curr_status == "target" && current_prot_is_decoy)
        {
          peptide_hit_.setMetaValue("target_decoy", "target+decoy");
        }
        else if (curr_status == "decoy" && !current_prot_is_decoy)
        {
          peptide_hit_.setMetaValue("target_decoy", "target+decoy");
        }

        hit.setMetaValue("target_decoy", current_prot_is_decoy ? "decoy" : "target");
      }
      peptide_hit_.addPeptideEvidence(pe);

      // depending on the numbering scheme used in the pepXML, "search_id_"
      // may appear to be "out of bounds" - see NOTE above:
      current_proteins_[min(UInt(current_proteins_.size()), search_id_) - 1]->insertHit(hit);
    }
    else if (element == "search_result") // parent: "spectrum_query"
    {
      // creates a new PeptideIdentification
      current_peptide_ = PeptideIdentification();
      current_peptide_.setRT(rt_);
      current_peptide_.setMZ(mz_);
      current_peptide_.setBaseName(current_base_name_);

      search_id_ = 1; // default if attr. is missing (ref. to "search_summary")
      optionalAttributeAsUInt_(search_id_, attributes, "search_id");
      // depending on the numbering scheme used in the pepXML, "search_id_"
      // may appear to be "out of bounds" - see NOTE above:
      String identifier = current_proteins_[min(UInt(current_proteins_.size()), search_id_) - 1]->getIdentifier();
      current_peptide_.setIdentifier(identifier);

      // set optional attributes
      if (!native_spectrum_name_.empty() && keep_native_name_)
      {
        current_peptide_.setMetaValue("pepxml_spectrum_name", native_spectrum_name_);
      }
      //TODO: we really need something uniform here, like scan number - and not in metainfointerface
      if (SpectrumLookup::isNativeID(native_spectrum_name_))
      {
        current_peptide_.setSpectrumReference( native_spectrum_name_);
      }
      else if (scannr_ != 0)
      {
        current_peptide_.setSpectrumReference( String("scan=") + String(scannr_));
      }
      //TODO else error?
      
        
      if (!experiment_label_.empty())
      {
        current_peptide_.setExperimentLabel(experiment_label_);
      }
      if (!swath_assay_.empty())
      {
        current_peptide_.setMetaValue("swath_assay", swath_assay_);
      }
      if (!status_.empty())
      {
        current_peptide_.setMetaValue("status", status_);
      }
    }
    else if (element == "spectrum_query") // parent: "msms_run_summary"
    {
      // sample:
      // <spectrum_query spectrum="foobar.02552.02552.2" start_scan="2552"
      // end_scan="2552" precursor_neutral_mass="1168.6176" assumed_charge="2"
      // index="10" retention_time_sec="488.652" experiment_label="urine"
      // swath_assay="EIVLTQSPGTL2:9" status="target">

      readRTMZCharge_(attributes); // sets "rt_", "mz_", "charge_"

      // retrieve optional attributes
      native_spectrum_name_ = "";
      experiment_label_ = "";
      swath_assay_ = "";
      status_ = "";
      optionalAttributeAsString_(native_spectrum_name_, attributes, "spectrum"); //TODO store separately? Do we ever need the pepXML internal ref?
      optionalAttributeAsString_(native_spectrum_name_, attributes, "spectrumNativeID"); //some engines write that optional attribute - is preferred to spectrum
      optionalAttributeAsString_(native_spectrum_name_, attributes, "native_id"); //MSFragger specific native ID
      optionalAttributeAsString_(experiment_label_, attributes, "experiment_label");
      optionalAttributeAsString_(swath_assay_, attributes, "swath_assay");
      optionalAttributeAsString_(status_, attributes, "status");
    }
    else if (element == "analysis_result") // parent: "search_hit"
    {
      current_analysis_result_ = PeptideHit::PepXMLAnalysisResult();
      current_analysis_result_.score_type = attributeAsString_(attributes, "analysis");
    }
    else if (element == "search_score_summary")
    {
      search_score_summary_ = true;
    }
    else if (element == "parameter") // parent: "search_score_summary"
    {
      // If we are within a search_score_summary, add the read in values to the current AnalysisResult
      if (search_score_summary_)
      {
        String name = attributeAsString_(attributes, "name");
        double value = attributeAsDouble_(attributes, "value");
        current_analysis_result_.sub_scores[name] = value;
      }
      else if (search_summary_)
      {
        String name = attributeAsString_(attributes, "name");
        if (name == "fragment_bin_tol")
        {
          double value = attributeAsDouble_(attributes, "value");
          params_.fragment_mass_tolerance = value/2.0;
          params_.fragment_mass_tolerance_ppm = false;
        }
        else if (name == "peptide_mass_tolerance")
        {
          double value = attributeAsDouble_(attributes, "value");
          params_.precursor_mass_tolerance = value;
        }
        // this is quite comet specific, but so is parameter name peptide_mass_units, see comet configuration file documentation
        else if (name == "peptide_mass_units")
        {
          int value = attributeAsInt_(attributes, "value");
          switch (value) {
          case 0: // comet value 0 type amu
              params_.precursor_mass_tolerance_ppm = false;
              break;
          case 1: // comet value 1 type mmu
              params_.precursor_mass_tolerance_ppm = false;
              break;
          case 2: // comet value 1 type ppm
              params_.precursor_mass_tolerance_ppm = true;
              break;
          default:
              break;
          }
        }
        else if (name == "decoy_search")
        {
          has_decoys_ = attributeAsInt_(attributes, "value");
        }
        else if (name == "decoy_prefix")
        {
          decoy_prefix_ = attributeAsString_(attributes, "value");
        }
      }
      else
      {
        // currently not handled
      }
    }
    else if (element == "peptideprophet_result") // parent: "analysis_result" (in "search_hit")
    {
      // PeptideProphet probability overwrites original search score
      // maybe TODO: deal with meta data associated with PeptideProphet search
      double value = attributeAsDouble_(attributes, "probability");
      if (current_peptide_.getScoreType() != "InterProphet probability")
      {
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType("PeptideProphet probability");
        current_peptide_.setHigherScoreBetter(true);
      }
      current_analysis_result_.main_score = value;
      current_analysis_result_.higher_is_better = true;
    }
    else if (element == "interprophet_result") // parent: "analysis_result" (in "search_hit")
    {
      // InterProphet probability overwrites PeptideProphet probability and
      // original search score
      double value = attributeAsDouble_(attributes, "probability");
      peptide_hit_.setScore(value);
      current_peptide_.setScoreType("InterProphet probability");
      current_peptide_.setHigherScoreBetter(true);
      current_analysis_result_.main_score = value;
      current_analysis_result_.higher_is_better = true;
    }
    else if (element == "modification_info") // parent: "search_hit" (in "search result")
    {
      // Has N-Term Modification
      double mod_nterm_mass;
      if (optionalAttributeAsDouble_(mod_nterm_mass, attributes, "mod_nterm_mass")) // this specifies a terminal modification
      {
        // look up the modification in the search_summary by mass
        bool found = false;
        for (const AminoAcidModification& it : variable_modifications_)
        {
          if ((fabs(mod_nterm_mass - it.getMass()) < mod_tol_) && it.getTerminus() == "n")
          {
            current_modifications_.emplace_back(it.getRegisteredMod(), Size(-1)); // position not needed for terminus
            found = true;
            break; // only one modification should match, so we can stop the loop here
          }
        }
        if (!found)
        {
          for (const AminoAcidModification& it : fixed_modifications_)
          {
            if ((fabs(mod_nterm_mass - it.getMass()) < mod_tol_) && it.getTerminus() == "n")
            {
              current_modifications_.emplace_back(it.getRegisteredMod(), Size(-1)); // position not needed for terminus
              found = true;
              break; // only one modification should match, so we can stop the loop here
            }
          }
        }

        if (!found)
        {
          // It was not registered in the pepXML header. Search it in DB by massdiff = mod_nterm_mass - orig_nterm_mass.
          std::vector<const ResidueModification*> mods;
          try
          {
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, mod_nterm_mass - Residue::getInternalToNTerm().getMonoWeight(), mod_tol_, "", ResidueModification::N_TERM);
          }
          catch(...)
          {
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, mod_nterm_mass - Residue::getInternalToNTerm().getMonoWeight(), mod_tol_, "", ResidueModification::PROTEIN_N_TERM);
          }

          if (!mods.empty())
          {
            current_modifications_.emplace_back(mods[0], Size(-1)); // -1, because position does not matter
          }
          else
          {
            error(LOAD, String("Cannot find N-terminal modification with mass " + String(mod_nterm_mass) + "."));
          }
          //TODO we should create and register a mod here
        }
      }

      // Has C-Term Modification
      double mod_cterm_mass;
      if (optionalAttributeAsDouble_(mod_cterm_mass, attributes, "mod_cterm_mass")) // this specifies a terminal modification
      {
        // look up the modification in the search_summary by mass
        bool found = false;
        for (const AminoAcidModification& amino : variable_modifications_)
        {
          if ((fabs(mod_cterm_mass - amino.getMass()) < mod_tol_) && amino.getTerminus() == "c")
          {
            current_modifications_.emplace_back(amino.getRegisteredMod(), Size(-1)); // position not needed for terminus
            found = true;
            break; // only one modification should match, so we can stop the loop here
          }
        }
        if (!found)
        {
          for (const AminoAcidModification& amino : fixed_modifications_)
          {
            if ((fabs(mod_cterm_mass - amino.getMass()) < mod_tol_) && amino.getTerminus() == "c")
            {
              current_modifications_.emplace_back(amino.getRegisteredMod(), Size(-1)); // position not needed for terminus
              found = true;
              break; // only one modification should match, so we can stop the loop here
            }
          }
        }

        if (!found)
        {
          // It was not registered in the pepXML header. Search it in DB.
          std::vector<const ResidueModification*> mods;
          try
          {
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, mod_cterm_mass - Residue::getInternalToCTerm().getMonoWeight(), mod_tol_, "", ResidueModification::C_TERM);
          }
          catch(...)
          {
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, mod_cterm_mass - Residue::getInternalToCTerm().getMonoWeight(), mod_tol_, "", ResidueModification::PROTEIN_C_TERM);
          }

          if (!mods.empty())
          {
            current_modifications_.emplace_back(mods[0], Size(-1)); // -1, because position does not matter
          }
          else
          {
            error(LOAD, String("Cannot find N-terminal modification with mass " + String(mod_nterm_mass) + "."));
          }
          //TODO we should create and register a mod here
        }
      }
    }
    else if (element == "alternative_protein") // parent: "search_hit"
    {
      String protein = attributeAsString_(attributes, "protein");
      PeptideEvidence pe;
      pe.setProteinAccession(protein);

      ProteinHit hit;
      hit.setAccession(protein);

      if (has_decoys_)
      {
        String curr_status("");
        bool current_prot_is_decoy = protein.hasPrefix(decoy_prefix_);
        if (peptide_hit_.metaValueExists("target_decoy"))
        {
          curr_status = peptide_hit_.getMetaValue("target_decoy");
        }
        if (curr_status.empty())
        {
          peptide_hit_.setMetaValue("target_decoy", current_prot_is_decoy ? "decoy" : "target");
        }
        else if (curr_status == "target" && current_prot_is_decoy)
        {
          peptide_hit_.setMetaValue("target_decoy", "target+decoy");
        }
        else if (curr_status == "decoy" && !current_prot_is_decoy)
        {
          peptide_hit_.setMetaValue("target_decoy", "target+decoy");
        }
        hit.setMetaValue("target_decoy", current_prot_is_decoy ? "decoy" : "target");
      }
      peptide_hit_.addPeptideEvidence(pe);

      // depending on the numbering scheme used in the pepXML, "search_id_"
      // may appear to be "out of bounds" - see NOTE above:
      current_proteins_[min(UInt(current_proteins_.size()), search_id_) - 1]->
      insertHit(hit);
    }
    else if (element == "mod_aminoacid_mass") // parent: "modification_info" (in "search_hit")
    {
      // this element should only be used for internal AA mods OR Terminal mods at a specific
      // amino acid (pepXML limitation)
      double modification_mass = attributeAsDouble_(attributes, "mass");
      Size modification_position = attributeAsInt_(attributes, "position");

      // the modification position is 1-based
      String origin = String(current_sequence_[modification_position - 1]);
      String temp_description = "";

      //TODO can we infer fixed/variable from the static/variable (diffmass) attributes in pepXML?
      // Only in some cases probably, since it is an optional attribute
      bool found = lookupAddFromHeader_(modification_mass, modification_position - 1, fixed_modifications_);

      if (!found)
      {
        found = lookupAddFromHeader_(modification_mass, modification_position - 1, variable_modifications_);
      }

      if (!found)
      {
        // try by PSI mod ID if present
        String psimod_id;
        optionalAttributeAsString_(psimod_id, attributes, "id");
        if (!psimod_id.empty())
        {
          try
          {
            current_modifications_.emplace_back(ModificationsDB::getInstance()->getModification(psimod_id, origin), modification_position - 1);
            found = true;
          }
          catch (...) {}
        }

        if (!found)
        {
          // Lookup in our DB
          //TODO also here, maybe the static/variable attribute is better for diffmass if present?
          double diffmass = modification_mass - ResidueDB::getInstance()->getResidue(origin)->getMonoWeight(Residue::Internal);
          vector<const ResidueModification*> mods;
          // try least specific search first:
          ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, diffmass, mod_tol_, origin, ResidueModification::ANYWHERE);
          if (mods.empty())
          {
            if (modification_position == 1)
            {
              ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, diffmass, mod_tol_, origin, ResidueModification::N_TERM);
              if (mods.empty()) ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, diffmass, mod_tol_, origin, ResidueModification::PROTEIN_N_TERM);
            }
            else if (modification_position == current_sequence_.length())
            {
              ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, diffmass, mod_tol_, origin, ResidueModification::C_TERM);
              if (mods.empty()) ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, diffmass, mod_tol_, origin, ResidueModification::PROTEIN_C_TERM);
            }
          }
          if (!mods.empty())
          {
            if (mods.size() > 1)
            {
              warning(LOAD, String("Modification '") + String(modification_mass) + "' of residue " + String(origin) + " at position "
              + String(modification_position) + " in '" + current_sequence_ + "' not registered in pepXML header nor uniquely defined in DB." +
              " Using " + mods[0]->getFullId());
              current_modifications_.emplace_back(mods[0], modification_position - 1);
            }
          }
          else
          {
            // still nothing found, register as unknown as last resort
            current_modifications_.emplace_back(
                ResidueModification::createUnknownFromMassString(
                    String(modification_mass),
                    modification_mass,
                    false,
                    ResidueModification::ANYWHERE, // since it is at an amino acid it is probably NOT a terminal mod
                    ResidueDB::getInstance()->getResidue(current_sequence_[modification_position - 1])),
                modification_position - 1
            );
          }
        }
      }
    }
    else if (element == "aminoacid_modification" || element == "terminal_modification") // parent: "search_summary"
    {
      String description;
      optionalAttributeAsString_(description, attributes, "description");
      String massdiff = attributeAsString_(attributes, "massdiff");
      String mass = attributeAsDouble_(attributes, "mass");
      String protein_terminus_entry;
      String aminoacid;
      String terminus;
      String is_variable = attributeAsString_(attributes, "variable");

      if (element == "aminoacid_modification")
      {
        aminoacid = attributeAsString_(attributes, "aminoacid");
        optionalAttributeAsString_(terminus, attributes, "peptide_terminus");
        optionalAttributeAsString_(protein_terminus_entry, attributes, "protein_terminus");
        // can't set term_spec to "ANYWHERE", because terminal mods. may not be registered as such
        // (this is an issue especially with Comet)!
      }
      else //terminal_modification
      {
        // Somehow very small fixed modifications (electron mass?) get annotated by X!Tandem. Don't add them as they interfere with other mods.
        if (fabs(massdiff.toDouble()) < xtandem_artificial_mod_tol_)
        {
          return;
        }

        optionalAttributeAsString_(aminoacid, attributes, "aminoacid");
        terminus = attributeAsString_(attributes, "terminus");
        // WARNING: Many search engines annotate this wrong! This is supposed to be empty or "n" or "c".
        // But many SEs annotate it as "Y" and "N" for Yes and No. We handle this in the AAMod ctor
        protein_terminus_entry = attributeAsString_(attributes, "protein_terminus");
      }

      AminoAcidModification aa_mod{
        aminoacid, massdiff, mass, is_variable, description, terminus, protein_terminus_entry,
        preferred_fixed_modifications_, preferred_variable_modifications_, mod_tol_
      };

      const vector<String>& errs = aa_mod.getErrors();
      if (!errs.empty())
      {
        error(LOAD, "Errors during parsing of aminoacid/terminal modification element:");
        for (const auto& e : errs)
        {
          error(LOAD, e);
        }
      }

      if (aa_mod.getRegisteredMod() != nullptr)
      {
        if (aa_mod.isVariable())
        {
          variable_modifications_.push_back(aa_mod);
          params_.variable_modifications.push_back(aa_mod.getDescription());
        }
        else
        {
          fixed_modifications_.push_back(aa_mod);
          params_.fixed_modifications.push_back(aa_mod.getDescription());
        }
      }
    }
    else if (element == "search_summary") // parent: "msms_run_summary"
    { // creates a new ProteinIdentification (usually)
      search_summary_ = true;
      current_base_name_ = "";
      optionalAttributeAsString_(current_base_name_, attributes, "base_name");
      if (!checked_base_name_) // work-around for files exported by Mascot
      {
        if (current_base_name_.hasSuffix(exp_name_))
        {
          seen_experiment_ = true;
        }
        else // wrong experiment after all - roll back changes that were made
        {
          proteins_->pop_back();
          current_proteins_.clear();
          wrong_experiment_ = true;
          return;
        }
      }

      fixed_modifications_.clear();
      variable_modifications_.clear();
      params_ = ProteinIdentification::SearchParameters();
      params_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_));
      String mass_type = attributeAsString_(attributes, "precursor_mass_type");
      if (mass_type == "monoisotopic")
      {
        hydrogen_mass_ = hydrogen_.getMonoWeight();
      }
      else
      {
        hydrogen_mass_ = hydrogen_.getAverageWeight();
        if (mass_type != "average")
        {
          error(LOAD, "'precursor_mass_type' attribute of 'search_summary' tag should be 'monoisotopic' or 'average', not '" + mass_type + "' (assuming 'average')");
        }
      }
      // assuming "SearchParameters::mass_type" refers to the fragment mass
      mass_type = attributeAsString_(attributes, "fragment_mass_type");
      if (mass_type == "monoisotopic")
      {
        params_.mass_type = ProteinIdentification::MONOISOTOPIC;
      }
      else
      {
        if (mass_type == "average")
        {
          params_.mass_type = ProteinIdentification::AVERAGE;
        }
        else
        {
          error(LOAD, "'fragment_mass_type' attribute of 'search_summary' tag should be 'monoisotopic' or 'average', not '" + mass_type + "'");
        }
      }

      search_engine_ = attributeAsString_(attributes, "search_engine");

      //TODO why do we not store versions?
      if (search_engine_ == "X! Tandem")
      {
        String ver;
        optionalAttributeAsString_(ver, attributes, "search_engine_version");
        if (ver.hasPrefix("MSFragger"))
        {
          search_engine_ = "MSFragger";
        }
      }

      // generate a unique identifier for every search engine run.
      prot_id_ = search_engine_ + "_" + date_.getDate() + "_" + date_.getTime();

      search_id_ = 1;
      optionalAttributeAsUInt_(search_id_, attributes, "search_id");
      vector<ProteinIdentification>::iterator prot_it;
      if (search_id_ <= proteins_->size()) // ProteinIdent. was already created for "msms_run_summary" -> add to it
      {
        prot_it = current_proteins_.back();
      }
      else // create a new ProteinIdentification
      {
        ProteinIdentification protein;
        protein.setDateTime(date_);
        proteins_->push_back(protein);
        prot_it = --proteins_->end();
        prot_id_ = prot_id_ + "_" + search_id_; // make sure the ID is unique
      }
      prot_it->setSearchEngine(search_engine_);
      prot_it->setIdentifier(prot_id_);
    }
    else if (element == "sample_enzyme") // parent: "msms_run_summary"
    { // special case: search parameter that occurs *before* "search_summary"!
      enzyme_ = attributeAsString_(attributes, "name");

      if (enzyme_ == "stricttrypsin") enzyme_ = "Trypsin/P"; // MSFragger synonyme

      if (ProteaseDB::getInstance()->hasEnzyme(enzyme_.toLower()))
      {
        params_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_));
      }
    }
    else if (element == "specificity" && params_.digestion_enzyme.getName() == "unknown_enzyme") // parent: "sample_enzyme"
    { // special case: search parameter that occurs *before* "search_summary"!
      String cut_before = attributeAsString_(attributes, "cut");
      String no_cut_after = attributeAsString_(attributes, "no_cut");
      String sense = attributeAsString_(attributes, "sense");
      params_.digestion_enzyme = DigestionEnzymeProtein(DigestionEnzyme(
          "user-defined," + enzyme_ + "," + cut_before + "," + no_cut_after + "," + sense,
          cut_before, no_cut_after, sense));
    }
    else if (element == "enzymatic_search_constraint") // parent: "search_summary"
    {
      //TODO we should not overwrite the enzyme here! Luckily in most files it is the same
      // enzyme as in sample_enzyme or something useless like "default".
      ///<enzymatic_search_constraint enzyme="nonspecific" max_num_internal_cleavages="1" min_number_termini="2"/>
      enzyme_ = attributeAsString_(attributes, "enzyme");    
      if (enzyme_ == "stricttrypsin") enzyme_ = "Trypsin/P"; // MSFragger synonyme

      if (ProteaseDB::getInstance()->hasEnzyme(enzyme_))
      {
        DigestionEnzymeProtein enzyme_to_set = *(ProteaseDB::getInstance()->getEnzyme(enzyme_.toLower()));
        if (params_.digestion_enzyme.getName() != "unknown_enzyme" && enzyme_to_set != params_.digestion_enzyme)
        {
          error(LOAD, "More than one enzyme found. This is currently not supported. Proceeding with last encountered only.");
        }
        params_.digestion_enzyme = enzyme_to_set;
      }

      // Always add the additional infos (e.g. for enzymatic_search_constraint enzyme="default" like in MSFragger)
      //TODO maybe check for keyword names like "default"?
      //TODO the question is, what to do with the following infos if the enzymes do not match?
      // We always parse and set it for now, to be
      // a) backwards-compatible
      // b) it's the only info we have
      int mc = attributeAsInt_(attributes, "max_num_internal_cleavages");
      params_.missed_cleavages = mc;

      int min_num_termini = attributeAsInt_(attributes, "min_number_termini");
      if (min_num_termini >= 0 && min_num_termini <= 2)
      {
        params_.enzyme_term_specificity = EnzymaticDigestion::Specificity(min_num_termini);
      }
    }
    else if (element == "search_database") // parent: "search_summary"
    {
      params_.db = attributeAsString_(attributes, "local_path");
      if (params_.db.empty())
      {
        optionalAttributeAsString_(params_.db, attributes, "database_name");
      }
    }
    else if (element == "msms_pipeline_analysis") // root
    {
      String date = attributeAsString_(attributes, "date");
      // fix for corrupted xs:dateTime format:
      if ((date[4] == ':') && (date[7] == ':') && (date[10] == ':'))
      {
        error(LOAD, "Format of attribute 'date' in tag 'msms_pipeline_analysis' does not comply with standard 'xs:dateTime'");
        date[4] = '-';
        date[7] = '-';
        date[10] = 'T';
      }
      date_ = asDateTime_(date);
    }
  }

  bool PepXMLFile::lookupAddFromHeader_(double modification_mass,
                                        Size modification_position,
                                        vector<AminoAcidModification> const& header_mods)
  {
    bool found = false;
    for (const auto& m : header_mods)
    {
      //TODO do we really want to allow another tolerance to the masses in the header?
      if (fabs(modification_mass - m.getMass()) < mod_tol_)
      {
        if (m.getAminoAcid().hasSubstring(current_sequence_[modification_position]))
        {
          current_modifications_
              .emplace_back(m.getRegisteredMod(), modification_position); // position not needed for terminus
          found = true;
          break; // only one modification should match, so we can stop the loop here
        }
      }
    }
    return found;
  }

  void PepXMLFile::endElement(const XMLCh* const /*uri*/,
                              const XMLCh* const /*local_name*/,
                              const XMLCh* const qname)
  {
    String element = sm_.convert(qname);

    // cout << "End: " << element << "\n";

    if (element == "analysis_summary")
    {
      analysis_summary_ = false;
    }
    else if (element == "search_score_summary")
    {
      search_score_summary_ = false;
    }
    else if (element == "analysis_result") // parent: "search_hit"
    {
      peptide_hit_.addAnalysisResults(current_analysis_result_);
    }
    else if (wrong_experiment_ || analysis_summary_)
    {
      // do nothing here (skip all elements that belong to the wrong experiment
      // or to an analysis summary)
    }
    else if (element == "spectrum_query")
    {
      // clear optional attributes
      native_spectrum_name_ = "";
      experiment_label_ = "";
      swath_assay_ = "";
      status_ = "";
    }
    else if (element == "search_hit")
    {
      AASequence temp_aa_sequence = AASequence::fromString(current_sequence_);

      //Note: using our AASequence::fromString on the modified_sequence of
      // the modification_info element is probably not possible since modifications may have special
      // symbols that we would need to save and lookup additionally
      //TODO one might argue that mods from XML attributes are better specified and should be preferred

      // Modifications explicitly annotated at the current search_hit take preference over
      // implicit fixed mods
      //TODO use peptide_hit_.getPeptideEvidences().back().getAAAfter()/Before() to see if Protein term is applicable at all
      // But: Leave it out for now since I bet there is software out there, not adhering to the specification
      for (const auto& mod : current_modifications_)
      {
        if (mod.first->getTermSpecificity() == ResidueModification::N_TERM ||
            mod.first->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
        {
          if (!temp_aa_sequence.hasNTerminalModification())
          {
            temp_aa_sequence.setNTerminalModification(mod.first);
          }
          else
          {
            warning(LOAD, "Multiple N-term mods specified for search_hit with sequence " + current_sequence_
              + " proceeding with first.");
          }
        }
        else if (mod.first->getTermSpecificity() == ResidueModification::C_TERM ||
                 mod.first->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
        {
          if (!temp_aa_sequence.hasCTerminalModification())
          {
            temp_aa_sequence.setCTerminalModification(mod.first);
          }
          else
          {
            warning(LOAD, "Multiple C-term mods specified for search_hit with sequence " + current_sequence_
              + " proceeding with first.");
          }
        }
        else // internal
        {
          if (!temp_aa_sequence[mod.second].isModified())
          {
            temp_aa_sequence.setModification(mod.second, mod.first->getFullId());
          }
          else
          {
            warning(LOAD, "Multiple mods for position " + String(mod.second)
              + " specified for search_hit with sequence " + current_sequence_ + " proceeding with first.");
          }
        }
      }

      // Now apply implicit fixed modifications at positions where there is no modification yet.
      for (const auto& mod : fixed_modifications_)
      {
        if (mod.getRegisteredMod()->getTermSpecificity() == ResidueModification::N_TERM ||
            mod.getRegisteredMod()->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
        {
          if (!temp_aa_sequence.hasNTerminalModification())
          {
            temp_aa_sequence.setNTerminalModification(mod.getRegisteredMod());
          }
          else
          {
            warning(LOAD, "Trying to add a fixed N-term modification from the search_summary to an already"
                          " annotated and modified N-terminus of " + current_sequence_
                          + " ... skipping.");
          }
        }
        else if (mod.getRegisteredMod()->getTermSpecificity() == ResidueModification::C_TERM ||
            mod.getRegisteredMod()->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
        {
          if (!temp_aa_sequence.hasCTerminalModification())
          {
            temp_aa_sequence.setCTerminalModification(mod.getRegisteredMod());
          }
          else
          {
            warning(LOAD, "Trying to add a fixed C-term modification from the search_summary to an already"
                          " annotated and modified N-terminus of " + current_sequence_
                          + " ... skipping.");
          }
        }
        else // go through the sequence and look for unmodified residues that match:
        {
          for (Size s = 0; s < temp_aa_sequence.size(); s++)
          {
            Residue const* r = &temp_aa_sequence[s];
            if (!r->isModified())
            {
              if (mod.getAminoAcid().hasSubstring(r->getOneLetterCode()))
              {
                // Note: this calls setModification_ on a new Residue which changes its
                // weight to the weight of the modification (set above)
                temp_aa_sequence.setModification(s,
                                                 ResidueDB::getInstance()->getModifiedResidue(r, mod.getRegisteredMod()->getFullId()));
              }
            }
            //TODO else warn as well? Skip for now.
          }
        }
      }

      peptide_hit_.setSequence(temp_aa_sequence);
      current_peptide_.insertHit(peptide_hit_);
    }
    else if (element == "search_result")
    {
      peptides_->push_back(current_peptide_);
    }
    else if (element == "search_summary")
    {
      // In idXML we only store search engine and date as identifier, but to distinguish two identification runs these values must be unique.
      // As a workaround to support multiple runs, we make the date unique by adding one second for every additional identification run.
      UInt hour, minute, second;
      date_.getTime(hour, minute, second);
      hour = (hour + (minute + (second + 1) / 60) / 60) % 24;
      minute = (minute + (second + 1) / 60) % 60;
      second = (second + 1) % 60;
      date_.setTime(hour, minute, second);

      current_proteins_.back()->setSearchParameters(params_);
      search_summary_ = false;
    }
  }

  void PepXMLFile::setPreferredFixedModifications(const std::vector<const ResidueModification*>& mods)
  {
    preferred_fixed_modifications_ = mods;
  }

  void PepXMLFile::setPreferredVariableModifications(const std::vector<const ResidueModification*>& mods)
  {
    preferred_variable_modifications_ = mods;
  }

  /*std::vector<int> PepXMLFile::getIsotopeErrorsFromIntSetting_(int intSetting)
  {
    static std::vector<std::vector<int>> cometIsoLists_{{},{1},{1,2},{1,2,3},{-8,-4,4,8},{-1,1,2,3}};
    return cometIsoLists_[intSetting];
  }
  */

} // namespace OpenMS
