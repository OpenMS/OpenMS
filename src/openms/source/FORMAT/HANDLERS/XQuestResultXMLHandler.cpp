// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XQuestResultXMLHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/CrossLinksDB.h>


#include <cassert>

using namespace std;
using namespace xercesc;

namespace OpenMS::Internal
{

    // Initialize static const members
    std::map< Size, String > XQuestResultXMLHandler::enzymes
    {
      std::make_pair(1, "trypsin"), std::make_pair(2, "chymotrypsin"), std::make_pair(3, "unknown_enzyme"),
      std::make_pair(9, "unknown_enzyme"), std::make_pair(10, "unknown_enzyme"), std::make_pair(14, "unknown_enzyme"),
      std::make_pair(15, "unknown_enzyme"), std::make_pair(16, "unknown_enzyme"), std::make_pair(17, "unknown_enzyme"),
      std::make_pair(18, "unknown_enzyme"), std::make_pair(20, "unknown_enzyme")
    };

    std::map< String, UInt> XQuestResultXMLHandler::months
    {
      std::make_pair("Jan", 1), std::make_pair("Feb", 2), std::make_pair("Mar", 3), std::make_pair("Apr", 4),
      std::make_pair("May", 5), std::make_pair("Jun", 6), std::make_pair("Jul", 7), std::make_pair("Aug", 8),
      std::make_pair("Sep", 9), std::make_pair("Oct", 10), std::make_pair("Nov", 11), std::make_pair("Dec", 12)
    };

    // reader
    XQuestResultXMLHandler::XQuestResultXMLHandler(const String &filename,
                                                   std::vector< PeptideIdentification > & pep_ids,
                                                   std::vector< ProteinIdentification > & prot_ids
                                                  ) :
      XMLHandler(filename, "1.0"),
      pep_ids_(&pep_ids),
      prot_ids_(&prot_ids),
      n_hits_(0),
      min_score_(0),
      max_score_(0)
    {
      // Initialize the one and only protein identification
      this->prot_ids_->clear();
      ProteinIdentification prot_id;
      prot_id.setSearchEngine("xQuest");
      prot_id.setSearchEngineVersion(VersionInfo::getVersion());
      prot_id.setMetaValue("SpectrumIdentificationProtocol", DataValue("MS:1002494")); // crosslinking search = MS:1002494
      this->prot_ids_->push_back(prot_id);

      // Fetch the enzymes database
      this->enzymes_db_ = ProteaseDB::getInstance();

      // TODO Produce some warnings that are associated with the reading of xQuest result files
      // OPENMS_LOG_WARN << "WARNING: Fixed modifications are not available in the xQuest input file and will thus be not present in the loaded data!\n" << std::endl;
    }

    // writer
    XQuestResultXMLHandler::XQuestResultXMLHandler(const std::vector<ProteinIdentification>& pro_id,
                                                   const std::vector<PeptideIdentification>& pep_id,
                                                   const String& filename,
                                                   const String& version
                                                 ) :
      XMLHandler(filename, version),
      pep_ids_(nullptr),
      prot_ids_(nullptr),
      cpro_id_(&pro_id),
      cpep_id_(&pep_id)
    {
    }

    XQuestResultXMLHandler::~XQuestResultXMLHandler()
    = default;

    void XQuestResultXMLHandler::extractDateTime_(const String & xquest_datetime_string, DateTime & date_time) const
    {
      StringList xquest_datetime_string_split;
      StringUtils::split(xquest_datetime_string,' ', xquest_datetime_string_split);
      if (this->is_openpepxl_)
      {
        // Example: 2017-03-17 23:04:50
        date_time.setDate(xquest_datetime_string_split[0]);
        date_time.setTime(xquest_datetime_string_split[1]);
      }
      else
      {
        // Example: Fri Dec 18 12:28:42 2015
        // Example: Fri Dec  8 12:28:42 2015
        // for single digit days, there are two spaces and the day is in the 4th string instead of the 3rd
        // that also moves the time and year one slot further
        UInt correction = 0;
        String day_string = xquest_datetime_string_split[2];
        if (day_string.empty())
        {
          correction = 1;
        }
        UInt day = xquest_datetime_string_split[2+correction].toInt();
        UInt year = xquest_datetime_string_split[4+correction].toInt();
        UInt month = XQuestResultXMLHandler::months[xquest_datetime_string_split[1]];
        date_time.setDate(month, day, year);
        date_time.setTime(xquest_datetime_string_split[3+correction]);
      }
    }

    // Extracts the position of the Cross-Link for intralinks and crosslinks
    void XQuestResultXMLHandler::getLinkPosition_(const xercesc::Attributes & attributes, std::pair<SignedSize, SignedSize> & pair)
    {
      String xlink_position = this->attributeAsString_(attributes, "xlinkposition");
      StringList xlink_position_split;
      StringUtils::split(xlink_position, "," ,xlink_position_split);

      pair.first = xlink_position_split[0].toInt();
      pair.second = xlink_position_split.size() == 2 ? xlink_position_split[1].toInt() : 0;
    }

    void XQuestResultXMLHandler::setPeptideEvidence_(const String & prot_string, PeptideHit & pep_hit)
    {
      StringList prot_list;
      StringUtils::split(prot_string, ",", prot_list);
      vector< PeptideEvidence > evidences;
      evidences.reserve(prot_list.size());

      for (StringList::const_iterator prot_list_it = prot_list.begin();
           prot_list_it != prot_list.end(); ++prot_list_it)
      {
        PeptideEvidence pep_ev;
        String accession = *prot_list_it;

        if (this->accessions_.find(accession) == this->accessions_.end())
        {
          this->accessions_.insert(accession);

          ProteinHit prot_hit;
          prot_hit.setAccession(accession);
          prot_hit.setMetaValue("target_decoy", accession.hasSubstring(this->decoy_string_) ? "decoy" : "target");

          (*this->prot_ids_)[0].getHits().push_back(prot_hit);
        }

        pep_ev.setProteinAccession(accession);
        pep_ev.setStart(PeptideEvidence::UNKNOWN_POSITION); // These information are not available in the xQuest result file
        pep_ev.setEnd(PeptideEvidence::UNKNOWN_POSITION);
        pep_ev.setAABefore(PeptideEvidence::UNKNOWN_AA);
        pep_ev.setAAAfter(PeptideEvidence::UNKNOWN_AA);

        evidences.push_back(pep_ev);
      }
      pep_hit.setPeptideEvidences(evidences);
    }

    // Assign all attributes in the peptide_id_attributes map to the MetaInfoInterface object
    void XQuestResultXMLHandler::addMetaValues_(MetaInfoInterface & meta_info_interface)
    {
      for (std::map<String, DataValue>::const_iterator it = this->peptide_id_meta_values_.begin();
           it != this->peptide_id_meta_values_.end(); ++it)
      {
        std::pair<String, DataValue> item = *it;
        meta_info_interface.setMetaValue(item.first, item.second);
      }
    }

    double XQuestResultXMLHandler::getMinScore() const
    {
      return this->min_score_;
    }

    double XQuestResultXMLHandler::getMaxScore() const
    {
      return this->max_score_;
    }

    UInt XQuestResultXMLHandler::getNumberOfHits() const
    {
      return this->n_hits_;
    }

    void XQuestResultXMLHandler::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname)
    {
      String tag = StringManager::convert(qname);
      if (tag == "xquest_results")
      {
        if (!this->is_openpepxl_)
        {
          ProteinIdentification::SearchParameters search_params((*this->prot_ids_)[0].getSearchParameters());
          search_params.charges = ListUtils::concatenate(this->charges_, ",");

          // min and max searched precursor charge not written out in xQuest
          // determination by charges in found results is not as clean, but is the best we can do
          search_params.setMetaValue("precursor:min_charge", this->min_precursor_charge_);
          search_params.setMetaValue("precursor:max_charge", this->max_precursor_charge_);

          (*this->prot_ids_)[0].setSearchParameters(search_params);
        }
      }
    }

    void XQuestResultXMLHandler::startElement(const XMLCh * const, const XMLCh * const, const XMLCh * const qname, const Attributes &attributes)
    {
      String tag = StringManager::convert(qname);
      // Extract meta information from the xquest_results tag
      if (tag == "xquest_results")
      {
        //cout << "Parse xQuest search settings" << endl;
        // Decide whether this Block is original xQuest or OpenPepXL
        String xquest_version = this->attributeAsString_(attributes, "xquest_version");
        this->is_openpepxl_ = xquest_version.hasSubstring("OpenPepXL");

        // Date and Time of Search
        DateTime date_time;
        //cout << "Parse Date" << endl;
        this->extractDateTime_(this->attributeAsString_(attributes, "date"), date_time);
        (*this->prot_ids_)[0].setDateTime(date_time);

        // Set the search parameters
        ProteinIdentification::SearchParameters search_params;

        //cout << "Parse Enzyme" << endl;
        // General
        if (this->is_openpepxl_) // Enzyme via name
        {
          search_params.digestion_enzyme = dynamic_cast<const DigestionEnzymeProtein&>(*this->enzymes_db_->getEnzyme(this->attributeAsString_(attributes, "enzyme_name")));
        }
        else // Enzyme via enzyme number in xQuest
        {
          search_params.digestion_enzyme = dynamic_cast<const DigestionEnzymeProtein&>(*this->enzymes_db_->getEnzyme(XQuestResultXMLHandler::enzymes[this->attributeAsInt_(attributes, "enzyme_num")]));
        }

        search_params.missed_cleavages = this->attributeAsInt_(attributes, "missed_cleavages");
        search_params.db = this->attributeAsString_(attributes, "database");
        search_params.precursor_mass_tolerance = this->attributeAsDouble_(attributes, "ms1tolerance");
        String tolerancemeasure_ms1 = this->attributeAsString_(attributes, this->is_openpepxl_ ? "tolerancemeasure_ms1" : "tolerancemeasure");
        search_params.precursor_mass_tolerance_ppm = tolerancemeasure_ms1 == "ppm";
        search_params.fragment_mass_tolerance = this->attributeAsDouble_(attributes, "ms2tolerance");
        String tolerancemeasure_ms2 = this->attributeAsString_(attributes, "tolerancemeasure_ms2");
        search_params.fragment_mass_tolerance_ppm = tolerancemeasure_ms2 != "Da";

        //cout << "Parse Mods" << endl;
        // variable Modifications
        vector< String > variable_mod_list;
        vector< String > variable_mod_split;
        String var_mod_string;
        if (this->optionalAttributeAsString_(var_mod_string, attributes, "variable_mod") && !var_mod_string.empty())
        {
          StringUtils::split(var_mod_string, ",", variable_mod_split);
          if (variable_mod_split[0].size() == 1 && variable_mod_split[0] != "0") // xQuest style mods=  "one-letter-code,mass"
          {
            double mod_mass = double(DataValue(variable_mod_split[1]));
            std::vector<String> mods;
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, mod_mass, 0.01, variable_mod_split[0]);
            if (!mods.empty())
            {
              variable_mod_list.push_back(mods[0]);
            }
          }
          search_params.variable_modifications = variable_mod_list;
        }
        // fixed Modifications
        String fixed_mod_string;
        StringList fixed_mod_list;
        if (this->optionalAttributeAsString_(fixed_mod_string, attributes, "fixed_mod") && !fixed_mod_string.empty())
        {
          fixed_mod_list = ListUtils::create<String>(fixed_mod_string);
          search_params.fixed_modifications = fixed_mod_list;
        }

        //cout << "Parse decoy stuff" << endl;
        String decoy_prefix;
        // if this info is not available, we can assume the decoy string is a prefix, since that is the standard way
        if (!this->optionalAttributeAsString_(decoy_prefix, attributes, "decoy_prefix"))
        {
          decoy_prefix = "1";
        }

        // change the default decoy string, if the parameter is given
        String current_decoy_string;
        if (this->optionalAttributeAsString_(current_decoy_string, attributes, "decoy_string") && !current_decoy_string.empty())
        {
          this->decoy_string_ = current_decoy_string;
        }

        // do some stringstream magic to turn "1" or "0" strings into booleans
        bool decoy_prefix_bool;
        std::istringstream is(decoy_prefix);
        is >> decoy_prefix_bool;

        // Meta Values
        String database_dc;
        if (this->optionalAttributeAsString_(database_dc, attributes, "database_dc") && !database_dc.empty())
        {
          search_params.setMetaValue("input_decoys", DataValue(database_dc));
        }
        search_params.setMetaValue("decoy_prefix", DataValue(decoy_prefix_bool));
        search_params.setMetaValue("decoy_string", DataValue(this->decoy_string_));

        search_params.setMetaValue("fragment:mass_tolerance_xlinks", DataValue(this->attributeAsDouble_(attributes, "xlink_ms2tolerance")));
        String monolink_masses_string_raw;
        StringList monolink_masses_string;

        if (this->optionalAttributeAsString_(monolink_masses_string_raw, attributes, "monolinkmw") && !monolink_masses_string_raw.empty())
        {
          monolink_masses_string = ListUtils::create<String>(monolink_masses_string_raw);
        }

        if (!monolink_masses_string.empty())
        {
          DoubleList monolink_masses;
          for (String monolink_string : monolink_masses_string)
          {
            monolink_masses.push_back(monolink_string.trim().toDouble());
          }
          search_params.setMetaValue("cross_link:mass_monolink", monolink_masses);
        }

        // xQuest uses the non-standard character "\u2212" for the minus in negative numbers. This can happen for zero-length cross-linkers.
        // Replace it with a proper "-" (minus), if there is one, to be able to convert it to a negative double.
        String xlinkermw = this->attributeAsString_(attributes, "xlinkermw");
        xlinkermw.substitute("\u2212", "-");

        search_params.setMetaValue("cross_link:mass_mass", DataValue(xlinkermw.toDouble()));
        this->cross_linker_name_ = this->attributeAsString_(attributes, "crosslinkername");
        search_params.setMetaValue("cross_link:name", DataValue(this->cross_linker_name_));
        String iso_shift = this->attributeAsString_(attributes, "cp_isotopediff");
        if (!iso_shift.empty())
        {
          search_params.setMetaValue("cross_link:mass_isoshift", iso_shift.toDouble());
        }

        bool ntermxlinkable;
        std::istringstream is_nterm(this->attributeAsString_(attributes, "ntermxlinkable"));
        is_nterm >> ntermxlinkable;

        //cout << "Parse AArequired" << endl;
        String aarequired;
        // older xQuest versions only allowed homobifunctional cross-linkers
        if (this->optionalAttributeAsString_(aarequired, attributes, "AArequired") && !aarequired.empty())
        {
          if (ntermxlinkable)
          {
            aarequired += ",N-term";
          }
          search_params.setMetaValue("cross_link:residue1", ListUtils::create<String>(aarequired));
          search_params.setMetaValue("cross_link:residue2", ListUtils::create<String>(aarequired));
        }
        else
        {
          String aarequired1 = this->attributeAsString_(attributes, "AArequired1");
          String aarequired2 = this->attributeAsString_(attributes, "AArequired2");
          if (ntermxlinkable)
          {
            if ( !(aarequired1.hasSubstring("N-term") || aarequired2.hasSubstring("N-term")) )
            {
              aarequired1 += ",N-term";
              aarequired2 += ",N-term";
            }
          }
          search_params.setMetaValue("cross_link:residue1", ListUtils::create<String>(aarequired1));
          search_params.setMetaValue("cross_link:residue2", ListUtils::create<String>(aarequired2));
        }

        if (this->is_openpepxl_)
        {
          //cout << "Parse OPXL specific settings" << endl;
          String searched_charges = this->attributeAsString_(attributes, "charges");
          search_params.charges = searched_charges;
          IntList charge_ints = ListUtils::create<Int>(searched_charges);
          std::sort(charge_ints.begin(), charge_ints.end());
          Int min_charge = charge_ints[0];
          Int max_charge = charge_ints.back();
          search_params.setMetaValue("precursor:min_charge", min_charge);
          search_params.setMetaValue("precursor:max_charge", max_charge);

          StringList ms_run = ListUtils::create<String>(this->attributeAsString_(attributes, "run_path"));
        }

        (*this->prot_ids_)[0].setSearchParameters(search_params);
      }
      else if (tag == "spectrum_search")
      {
        // Examples of lines to be parsed with this code
        // <spectrum_search spectrum="GUA1354-S15-A-LRRK2_DSG_A4.light.2616_GUA1354-S15-A-LRRK2_DSG_A4.heavy.2481" mz_precursor="590.556396484375" scantype="light_heavy" charge_precursor="4" Mr_precursor="2358.19648007042" rtsecscans="2231.988:2194.8258"                mzscans="590.556396484375:592.065673828125" >
        // <spectrum_search spectrum="GUA1354-S15-A-LRRK2_DSG_A4.light.1327_GUA1354-S15-A-LRRK2_DSG_A4.heavy.1327" mz_precursor="1008.83288574219" scantype="light"       charge_precursor="3" Mr_precursor="3023.47682782626" rtsecscans="2796.68020000002:2796.68020000002" mzscans="1008.83288574219:1008.83288574219" >
        // <spectrum_search Mr_precursor="1465.880913324" addedMass="0" apriori_pmatch_common="0.0311" apriori_pmatch_xlink="0.0658" charge_precursor="3" ionintensity_stdev="5.73" iontag_ncandidates="240" mean_ionintensity="2.28" mz_precursor="489.63479614" mzscans="489.63479614:493.6600647" ncommonions="71" nxlinkions="102" rtsecscans="2491:2477" scantype="light_heavy" spectrum="aleitner_M1012_006.c.02942.02942.3_aleitner_M1012_006.c.02913.02913.3">
        //cout << "Parse Spectrum" << endl;
        // Update retention time of light
        StringList rt_split;
        StringUtils::split(this->attributeAsString_(attributes, "rtsecscans"), ":", rt_split);
        this->rt_light_ = rt_split[0].toDouble();
        this->rt_heavy_ = rt_split[1].toDouble();

        String mz_scans = this->attributeAsString_(attributes, "mzscans");
        if (!mz_scans.empty())
        {
          StringList mz_split;
          StringUtils::split(this->attributeAsString_(attributes, "mzscans"), ":", mz_split);
          this->mz_light_ = mz_split[0].toDouble();
          this->mz_heavy_ = mz_split[1].toDouble();
        }
        else
        {
          double mz_precursor = this->attributeAsString_(attributes, "mz_precursor").toDouble();
          this->mz_light_ = mz_precursor;
          this->mz_heavy_ = mz_precursor;
        }

        // Update min and max precursor charge
        UInt charge_precursor = this->attributeAsInt_(attributes, "charge_precursor");
        if (!this->is_openpepxl_)
        {
          //cout << "Parse xQuest Spectrum" << endl;
          if (charge_precursor < this->min_precursor_charge_)
          {
            this->min_precursor_charge_ = charge_precursor;
          }
          if (charge_precursor > this->max_precursor_charge_)
          {
            this->max_precursor_charge_ = charge_precursor;
          }
          this->charges_.insert(charge_precursor);

          // read input filename (will not contain file type this way), example:
          // "C_Lee_141014_CRM_dialysis_NCE20_1.06904.06904.3_C_Lee_141014_CRM_dialysis_NCE20_1.06904.06904.3"
          String spectrum = this->attributeAsString_(attributes, "spectrum");

          // split into light and heavy (or light and light if unlabeled)
          StringList split_spectra = XQuestResultXMLHandler::splitByMiddle(spectrum, '_');
          vector<String> split_spectrum_light;
          vector<String> split_spectrum_heavy;
          // split away the spectrum indices from the file name
          StringUtils::split(split_spectra[0], ".", split_spectrum_light);
          StringUtils::split(split_spectra[1], ".", split_spectrum_heavy);


          String file_name = split_spectrum_light[0];
          if (std::find(this->ms_run_path_.begin(), this->ms_run_path_.end(), file_name) == this->ms_run_path_.end())
          {
            this->ms_run_path_.push_back(file_name);
          }
          this->spectrum_input_file_ = file_name;

          // read spectrum indices
          if (split_spectrum_light[split_spectrum_light.size()-1].size() > 1)
          {
            //cout << "Parse Spectrum index version 1" << endl;
            //cout << endl << split_spectrum_light[split_spectrum_light.size()-1] << endl;
            //cout << endl << split_spectrum_heavy[split_spectrum_heavy.size()-1] << endl;
            this->spectrum_index_light_ = split_spectrum_light[split_spectrum_light.size()-1].toInt();
            this->spectrum_index_heavy_ = split_spectrum_heavy[split_spectrum_heavy.size()-1].toInt();
          }
          else
          {
            //cout << "Parse Spectrum index version 2" << endl;
            //cout << endl << split_spectrum_light[split_spectrum_light.size()-2] << endl;
            //cout << endl << split_spectrum_heavy[split_spectrum_heavy.size()-2] << endl;
            this->spectrum_index_light_ = split_spectrum_light[split_spectrum_light.size()-2].toInt();
            this->spectrum_index_heavy_ = split_spectrum_heavy[split_spectrum_heavy.size()-2].toInt();
          }
        }
        else
        {
          this->spectrum_index_light_ = this->attributeAsInt_(attributes, "scan_index_light");
          this->spectrum_index_heavy_ = this->attributeAsInt_(attributes, "scan_index_heavy");

          ProteinIdentification::SearchParameters search_params((*this->prot_ids_)[0].getSearchParameters());
          if (!search_params.metaValueExists("input_mzML"))
          {
            String spectrum = this->attributeAsString_(attributes, "spectrum");
            vector<String> split_spectrum;
            StringUtils::split(spectrum, ".", split_spectrum);
            String file_name = split_spectrum[0];
            search_params.setMetaValue("input_mzML", file_name + String(".mzML"));
            (*this->prot_ids_)[0].setSearchParameters(search_params);
          }
        }
      }
      else if (tag == "search_hit")
      {
        // Examples of lines to be parsed with this code

        // <search_hit search_hit_rank="1" id="DNSTMGYMAAKK-RDVEKFLSK-a11-b5" type="xlink" structure="DNSTMGYMAAKK-RDVEKFLSK" seq1="DNSTM(Oxidation)GYM(Oxidation)AAKK" seq2="RDVEKFLSK" prot1="tr|Q8TBA7|Q8TBA7_HUMAN" prot2="sp|Q5S007-v1|LRRK2_HUMAN" topology="a11-b5" xlinkposition="11,5" Mr="2564.2250873787" mz="855.748972259671" charge="3" xlinkermass="96.0211294" measured_mass="2564.22762128328"
        // error="0.000844634859959115" error_rel="0.987012415251626" xlinkions_matched="6" backboneions_matched="1" xcorrx="0.312314444528579" xcorrb="-0.0506118717404067" match_odds="0.794234705691207" prescore="0.0369274467229843" num_of_matched_ions_alpha="3" num_of_matched_ions_beta="4" num_of_matched_common_ions_alpha="1" num_of_matched_common_ions_beta="0" num_of_matched_xlink_ions_alpha="2" num_of_matched_xlink_ions_beta="4"
        // TIC="0.0292408974147396" wTIC="0.026377408862402" intsum="0.397526955232024" HyperCommon="0.743940400979002" HyperXlink="34.1231158133129" HyperAlpha="16.0630790689233" HyperBeta="6.84199589723582" HyperBoth="31.1180197102582" selected="false" target_decoy="target" protein_references="unique" annotated_spec="" score="2.32103769126514" >

        // <search_hit search_hit_rank="3" id="MGIKTSEGTPGFRAPEVAR-HKMSYSGR-a4-b2" type="xlink" structure="MGIKTSEGTPGFRAPEVAR-HKMSYSGR" seq1="M(Oxidation)GIKTSEGTPGFRAPEVAR" seq2="HKMSYSGR" prot1="sp|Q5S007-v1|LRRK2_HUMAN" prot2="sp|Q5S007-v1|LRRK2_HUMAN" topology="a4-b2" xlinkposition="4,2" Mr="3079.4967874314" mz="770.881473324621" charge="4" xlinkermass="96.0211294" measured_mass="3079.49506405479"
        // error="-0.000430844152219834" error_rel="-0.558898049996855" xlinkions_matched="14" backboneions_matched="6" xcorrx="0.198434093695336" xcorrb="0.00514737154810852" match_odds="1.45901170826174" prescore="0.0599999986588955" num_of_matched_ions_alpha="15" num_of_matched_ions_beta="5" num_of_matched_common_ions_alpha="5" num_of_matched_common_ions_beta="1" num_of_matched_xlink_ions_alpha="10" num_of_matched_xlink_ions_beta="4"
        // TIC="0.0562770907575218" wTIC="0.0370273112047904" intsum="0.818966233637184" HyperCommon="6.80908719125821" HyperXlink="33.1079286508253" HyperAlpha="15.5319805998036" HyperBeta="1.62767939400878" HyperBoth="23.997840801109" selected="false" target_decoy="target" protein_references="unique" annotated_spec="" score="2.69829871110556" >

        // <search_hit Mr="2145.18339" TIC="0.08237" TIC_alpha="0.03287" TIC_beta="0.04951" annotated_spec="" apriori_match_probs="0.99970" apriori_match_probs_log="-0.00013" backboneions_matched="" charge="3" error="1.6" error_rel="-1.6" id="KSKTLQYFA-KQYSAKAK-a1-b1" intsum="91.91980" match_error_mean="-8.04546309837745" match_error_stdev="278.931294616457" match_odds="2.85579" match_odds_alphacommon="1.77210" match_odds_alphaxlink="1.98118"
        // match_odds_betacommon="2.35354" match_odds_betaxlink="5.31633" measured_mass="2145.1800" mz="716.06781" num_of_matched_common_ions_alpha="1" num_of_matched_common_ions_beta="1" num_of_matched_ions_alpha="3" num_of_matched_ions_beta="5" num_of_matched_xlink_ions_alpha="2" num_of_matched_xlink_ions_beta="4" prescore="0.11625" prescore_alpha="0.08108" prescore_beta="0.16667"
        // prot1="sp|O14126|PRS6A_SCHPO" prot2="decoy_reverse_sp|Q9UUB6|UBLH2_SCHPO" score="8.93" search_hit_rank="2" seq1="KSKTLQYFA" seq2="KQYSAKAK" series_score_mean="2.48843" structure="KSKTLQYFA-KQYSAKAK" topology="a1-b1" type="xlink" wTIC="0.01521" weighted_matchodds_mean="1.31713728336586" weighted_matchodds_sum="0.658568641682928" xcorrall="0.00000" xcorrb="0.05442" xcorrx="0.11647" xlinkermass="138.0680796" xlinkions_matched="" xlinkposition="1,1">

        PeptideIdentification peptide_identification;

        PeptideHit peptide_hit_alpha;
        PeptideHit peptide_hit_beta;
        vector<PeptideHit> peptide_hits;

        String seq1 = String(this->attributeAsString_(attributes, "seq1"));
        if (!this->is_openpepxl_)
        {
          seq1 = seq1.substitute("X", "M(Oxidation)");
        }

        // XL Type, determined by "type"
        String xlink_type_string = this->attributeAsString_(attributes, "type");

        AASequence alpha_seq = AASequence::fromString(seq1);
        std::pair<SignedSize, SignedSize> positions;
        this->getLinkPosition_(attributes, positions);
        int xl_pos = positions.first - 1;

        std::vector<String> mods;

        // xQuest uses the non-standard character "\u2212" for the minus in negative numbers. This can happen for zero-length cross-linkers.
        // Replace it with a proper "-" (minus), if there is one, to be able to convert it to a negative double.
        String xlinkermass_string = this->attributeAsString_(attributes, "xlinkermass");
        xlinkermass_string.substitute("\u2212", "-");
        double xl_mass = DataValue(xlinkermass_string.toDouble());

        ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, xl_mass, 0.01, alpha_seq[xl_pos].getOneLetterCode(), ResidueModification::ANYWHERE);

        if (xlink_type_string == "monolink")
        {
          if (!mods.empty())
          {
            bool mod_set = false;
            for (const String& mod : mods)
            {
              if (mod.hasSubstring(this->cross_linker_name_))
              {
                alpha_seq.setModification(xl_pos, mod);
                mod_set = true;
                // do not break to prioritize mods later in the list
                // the XLMODS results are further down the list and should be prioritized instead of unimod
                //break;
              }
            }
            if (!mod_set)
            {
              cout << "Set default mono-link: " << mods[0] << endl;
              alpha_seq.setModification(xl_pos, mods[0]);
            }
          }
        }
        else
        {
          peptide_hit_alpha.setMetaValue("xl_mod", this->cross_linker_name_);
        }

        peptide_hit_alpha.setSequence(alpha_seq);

        UInt charge = this->attributeAsInt_(attributes, "charge");
        peptide_hit_alpha.setCharge(charge);

        peptide_hit_alpha.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectrum_index_light_);
        peptide_hit_alpha.setMetaValue("spectrum_index", spectrum_index_light_);
        peptide_hit_alpha.setMetaValue("spectrum_input_file", spectrum_input_file_);

        String specIDs;
        if (spectrum_index_light_ != spectrum_index_heavy_)
        {
          specIDs = String(spectrum_index_light_) + "," + String(spectrum_index_heavy_);

          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT, this->rt_heavy_);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ, this->mz_heavy_);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF, spectrum_index_heavy_);
          peptide_hit_alpha.setMetaValue("spectrum_index_heavy", spectrum_index_heavy_);
        }
        else
        {
          specIDs = String(spectrum_index_light_);
        }
        peptide_identification.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, specIDs);

        // Set xl_chain meta value for alpha
        peptide_hit_alpha.setMetaValue("xl_chain", "MS:1002509");

        // Set Attributes of Peptide Identification
        peptide_identification.setMZ(this->mz_light_);
        peptide_identification.setRT(this->rt_light_);
        peptide_identification.setScoreType("OpenPepXL:score"); // Needed, since hard-coded in MzIdentMLHandler

        String prot1_string = this->attributeAsString_(attributes, "prot1");

        // Decide if decoy for alpha
        DataValue target_decoy = DataValue(prot1_string.hasSubstring(this->decoy_string_) ? "decoy" : "target");
        peptide_hit_alpha.setMetaValue("target_decoy", target_decoy);

        // Attributes of peptide_hit_alpha
        double score = this->attributeAsDouble_(attributes, "score");
        DataValue xlinkermass = DataValue(xl_mass);

        // Set minscore and maxscore encountered
        if (score < this->min_score_)
        {
          this->min_score_ = score;
        }
        if (score > this->max_score_)
        {
          this->max_score_ = score;
        }
        peptide_hit_alpha.setScore(score);

        peptide_hit_alpha.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM, DataValue(this->attributeAsDouble_(attributes, "error_rel")));

        // Get common attributes of Peptide Identification
        this->peptide_id_meta_values_["OpenPepXL:id"] = DataValue(this->attributeAsString_(attributes, "id"));
        this->peptide_id_meta_values_["OpenPepXL:xlinkermass"] = xlinkermass;
        this->peptide_id_meta_values_[Constants::UserParam::OPENPEPXL_XL_RANK] = DataValue(this->attributeAsInt_(attributes, "search_hit_rank"));
        this->peptide_id_meta_values_["OpenPepXL:score"] = DataValue(score);
        this->peptide_id_meta_values_["OpenPepXL:structure"] = DataValue(this->attributeAsString_(attributes, "structure"));

        // get scores (which might be optional)
        String wTIC, TIC, intsum, match_odds, fdr;
        if (this->optionalAttributeAsString_(wTIC, attributes, "wTIC") && !wTIC.empty())
        {
          this->peptide_id_meta_values_["OpenPepXL:wTIC"] = DataValue(wTIC.toDouble());
        }
        if (this->optionalAttributeAsString_(TIC, attributes, "TIC") && !TIC.empty())
        {
          this->peptide_id_meta_values_["OpenPepXL:percTIC"] = DataValue(TIC.toDouble());
        }

        if (this->optionalAttributeAsString_(intsum, attributes, "intsum") && !intsum.empty())
        {
          this->peptide_id_meta_values_["OpenPepXL:intsum"] = DataValue(intsum.toDouble());
        }
        if (this->optionalAttributeAsString_(match_odds, attributes, "match_odds") && !match_odds.empty())
        {
          this->peptide_id_meta_values_["OpenPepXL:match-odds"] = DataValue(match_odds.toDouble());
        }
        if (this->optionalAttributeAsString_(fdr, attributes, "fdr") && !fdr.empty())
        {
          this->peptide_id_meta_values_["XFDR:FDR"] = DataValue(fdr.toDouble());
        }

        String xprophet_f;
        if (this->optionalAttributeAsString_(xprophet_f, attributes, "xprophet_f") && !xprophet_f.empty())
        {
          peptide_hit_alpha.setMetaValue("XFDR:used_for_FDR", xprophet_f.toInt());
        }
        String fdr_type;
        if (this->optionalAttributeAsString_(xprophet_f, attributes, "fdr_type") && !fdr_type.empty())
        {
          peptide_hit_alpha.setMetaValue("XFDR:fdr_type", fdr_type);
        }

        assert(this->peptide_id_meta_values_["OpenPepXL:id"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenPepXL:xlinkermass"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_[Constants::UserParam::OPENPEPXL_XL_RANK] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenPepXL:score"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenPepXL:structure"] != DataValue::EMPTY);

        this->addMetaValues_(peptide_hit_alpha);

        // Store specific stuff for peptide hit alpha
        peptide_hit_alpha.setMetaValue("matched_linear_alpha",
                                        DataValue(this->attributeAsInt_(attributes, "num_of_matched_common_ions_alpha")));
        peptide_hit_alpha.setMetaValue("matched_xlink_alpha",
                                        DataValue(this->attributeAsInt_(attributes, "num_of_matched_xlink_ions_alpha")));
        peptide_hit_alpha.setMetaValue("matched_linear_beta",
                                        DataValue(this->attributeAsInt_(attributes, "num_of_matched_common_ions_beta")));
        peptide_hit_alpha.setMetaValue("matched_xlink_beta",
                                        DataValue(this->attributeAsInt_(attributes, "num_of_matched_xlink_ions_beta")));

        peptide_hit_alpha.setMetaValue("prot1", DataValue(prot1_string));
        peptide_hit_alpha.setMetaValue("prot2", DataValue("-"));
        peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS, DataValue("-"));
        peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS, xlinkermass);

        // Set peptide Evidences for Alpha (need one for each accession in the prot1_string)
        this->setPeptideEvidence_(prot1_string, peptide_hit_alpha);

        // Switch on Cross-link type
        if (xlink_type_string == "xlink")
        {
          // Set the cross Link Mass
          ProteinIdentification::SearchParameters search_params((*this->prot_ids_)[0].getSearchParameters());
          if (!search_params.metaValueExists("cross_link:mass"))
          {
            search_params.setMetaValue("cross_link:mass", DataValue(xlinkermass));
          }
          (*this->prot_ids_)[0].setSearchParameters(search_params);

          peptide_hit_beta.setScore(score);

          peptide_hit_beta.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM, DataValue(this->attributeAsDouble_(attributes, "error_rel")));

          String seq2 = String(this->attributeAsString_(attributes, "seq2"));
          if (!this->is_openpepxl_)
          {
            seq2 = seq2.substitute("X", "M(Oxidation)");
          }
          peptide_hit_beta.setSequence(AASequence::fromString(seq2));
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE, seq2);
          peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE, seq2);

          peptide_hit_beta.setCharge(charge);

          peptide_hit_beta.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectrum_index_light_);

          if (spectrum_index_light_ != spectrum_index_heavy_)
          {
            peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT, this->rt_heavy_);
            peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ, this->mz_heavy_);
            peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF, spectrum_index_heavy_);
            peptide_hit_beta.setMetaValue("spectrum_index_heavy", spectrum_index_heavy_);
          }

          this->addMetaValues_(peptide_hit_beta);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE, DataValue("cross-link"));
          peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE, DataValue("cross-link"));

          // Set xl positions, depends on xl_type
          this->getLinkPosition_(attributes, positions);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1, DataValue(positions.first - 1));
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, DataValue(positions.second - 1));
          peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1, DataValue(positions.first - 1));
          peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, DataValue(positions.second - 1));

          String term_spec_alpha("ANYWHERE"), term_spec_beta("ANYWHERE");
          StringList aarequired1 = search_params.getMetaValue("cross_link:residue1");
          StringList aarequired2 = search_params.getMetaValue("cross_link:residue2");

          // StringListUtils::searchSuffix(aarequired1, seq1[positions.first-1])
          if (positions.first == 1 && StringListUtils::searchSuffix(aarequired1, seq1[0]) == aarequired1.end())
          {
            term_spec_alpha = "N_TERM";
          }
          if (positions.second == 1 && StringListUtils::searchSuffix(aarequired2, seq2[0]) == aarequired2.end())
          {
            term_spec_beta = "N_TERM";
          }
          if (positions.first == SignedSize(seq1.size()) && StringListUtils::searchSuffix(aarequired1, seq1[positions.first-1]) == aarequired1.end())
          {
            term_spec_alpha = "C_TERM";
          }
          if (positions.second == SignedSize(seq2.size()) && StringListUtils::searchSuffix(aarequired2, seq2[positions.second-1]) == aarequired2.end())
          {
            term_spec_beta = "C_TERM";
          }
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA, term_spec_alpha);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA, term_spec_beta);

          // Protein
          String prot2_string = this->attributeAsString_(attributes, "prot2");

          // Decide if decoy for beta
          target_decoy = DataValue(prot2_string.hasSubstring(this->decoy_string_) ? "decoy" : "target");
          peptide_hit_beta.setMetaValue(Constants::UserParam::TARGET_DECOY, target_decoy);

          //  Set xl_chain meta value for beta
          peptide_hit_beta.setMetaValue("xl_chain", "MS:1002510");

          // Set peptide_hit specific stuff
          peptide_hit_beta.setMetaValue("matched_linear_alpha",
                                          DataValue(this->attributeAsInt_(attributes, "num_of_matched_common_ions_alpha")));
          peptide_hit_beta.setMetaValue("matched_xlink_alpha",
                                          DataValue(this->attributeAsInt_(attributes, "num_of_matched_xlink_ions_alpha")));
          peptide_hit_beta.setMetaValue("matched_linear_beta",
                                          DataValue(this->attributeAsInt_(attributes, "num_of_matched_common_ions_beta")));
          peptide_hit_beta.setMetaValue("matched_xlink_beta",
                                          DataValue(this->attributeAsInt_(attributes, "num_of_matched_xlink_ions_beta")));

          peptide_hit_alpha.setMetaValue("prot2", DataValue(prot2_string));
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS, DataValue(prot2_string));
          peptide_hit_beta.setMetaValue("prot1", DataValue(prot1_string));
          peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS, DataValue(prot2_string));
          peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS, xlinkermass);
          peptide_hit_beta.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD, this->cross_linker_name_);

          // Set Peptide Evidences for Beta
          this->setPeptideEvidence_(prot2_string, peptide_hit_beta);

          // TODO Determine if protein is intra/inter protein, check all protein ID combinations
          // StringList prot1_list;
          // prot1_string.split(",", prot1_list);
          // StringList prot2_list;
          // prot2_string.split( ",", prot2_list);

        }
        else if (xlink_type_string == "intralink")
        {
          // xl type
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE, DataValue("loop-link"));

          // Set xl positions, depends on xl_type
          this->getLinkPosition_(attributes, positions);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1, DataValue(positions.first - 1));
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, DataValue(positions.second - 1));
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS, xlinkermass);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD, this->cross_linker_name_);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE, "-");
        }
        else if (xlink_type_string == "monolink")
        {
          // TODO Set the xl_mass and xl_mod MetaValues instead
          // this->monolinks_masses_.insert(this->attributeAsDouble_(attributes, "xlinkermass"));

          // xl_type
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE, DataValue("mono-link"));

          std::pair< SignedSize, SignedSize > xlink_pos;
          this->getLinkPosition_(attributes, xlink_pos);
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1, DataValue(xlink_pos.first - 1));
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, DataValue("-"));
          peptide_hit_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE, "-");
        }
        else
        {
          OPENMS_LOG_ERROR << "ERROR: Unsupported Cross-Link type: " << xlink_type_string << endl;
          throw std::exception();
        }

        // Finalize this record
        peptide_hits.push_back(peptide_hit_alpha);

        if (peptide_hit_beta.metaValueExists(Constants::UserParam::OPENPEPXL_XL_POS1))
        {
          peptide_hits.push_back(peptide_hit_beta);
        }

        peptide_identification.setHits(peptide_hits);
        this->peptide_id_meta_values_.clear();
        this->pep_ids_->push_back(peptide_identification);
        this->n_hits_++;
      }
    }

    StringList XQuestResultXMLHandler::splitByNth(const String& input, const char separator, const Size n)
    {
      StringList list;
      Size current_index = 0;

      Size count = 0;
      while (current_index < input.size() && count < n)
      {
        current_index++;
        if (input.at(current_index) == separator)
        {
          count++;
        }
      }

      list.push_back(input.prefix(current_index));
      list.push_back(input.suffix(input.size()-current_index-1));

      return list;
    }

    StringList XQuestResultXMLHandler::splitByMiddle(const String& input, const char separator)
    {
      // count separators
      Size n = std::count(input.begin(), input.end(), separator);

      if ( (n != 0) && (n % 2 == 1) )
      {
        n = (n / 2) + 1;
        return splitByNth(input, separator, n);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The separator has to occur in the input string an uneven number of times (and at least once).");
      }
    }

    void XQuestResultXMLHandler::writeTo(std::ostream& os)
    {
      ProteinIdentification::SearchParameters search_params;
      search_params = (*this->cpro_id_)[0].getSearchParameters();

      String input_filename;
      if (search_params.metaValueExists("input_mzML"))
      {
        input_filename = search_params.getMetaValue("input_mzML");
      }
      String spec_xml_name = search_params.getMetaValue("out_xquest_specxml");

      os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
      os << "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>\n";

      DateTime time= DateTime::now();
      String timestring = time.getDate() + " " + time.getTime();

      String mono_masses = search_params.getMetaValue("cross_link:mass_monolink");
      mono_masses = mono_masses.substr(1).chop(1);

      String precursor_mass_tolerance_unit = search_params.precursor_mass_tolerance_ppm ? "ppm" : "Da";
      double precursor_mass_tolerance = search_params.precursor_mass_tolerance;
      String fragment_mass_tolerance_unit = search_params.fragment_mass_tolerance_ppm ? "ppm" : "Da";
      double fragment_mass_tolerance = search_params.fragment_mass_tolerance;
      double fragment_mass_tolerance_xlinks = search_params.getMetaValue("fragment:mass_tolerance_xlinks");

      String cross_link_name = search_params.getMetaValue("cross_link:name");
      double cross_link_mass_light = search_params.getMetaValue("cross_link:mass");
      double cross_link_mass_iso_shift = 0;
      if (search_params.metaValueExists("cross_link:mass_isoshift"))
      {
        cross_link_mass_iso_shift = search_params.getMetaValue("cross_link:mass_isoshift");
      }
      String aarequired1, aarequired2;
      aarequired1 = search_params.getMetaValue("cross_link:residue1");
      aarequired1 = aarequired1.substr(1).chop(1);
      aarequired2 = search_params.getMetaValue("cross_link:residue2");
      aarequired2 = aarequired2.substr(1).chop(1);
      bool ntermxlinkable = aarequired1.hasSubstring("N-term") || aarequired2.hasSubstring("N-term");

      String in_fasta = search_params.db;
      String in_decoy_fasta = search_params.getMetaValue("input_decoys");
      String enzyme_name = search_params.digestion_enzyme.getName();
      int missed_cleavages = search_params.missed_cleavages;

      StringList variable_mod_list = search_params.variable_modifications;
      String variable_mods;
      for (Size i = 0; i < variable_mod_list.size(); ++i)
      {
        variable_mods += variable_mod_list[i] + ",";
      }
      variable_mods = variable_mods.chop(1);

      StringList fixed_mod_list = search_params.fixed_modifications;
      String fixed_mods;
      for (Size i = 0; i < fixed_mod_list.size(); ++i)
      {
        fixed_mods += fixed_mod_list[i] + ",";
      }
      fixed_mods = fixed_mods.chop(1);

      String decoy_prefix = search_params.getMetaValue("decoy_prefix").toString();
      String decoy_string = search_params.getMetaValue("decoy_string").toString();

      String searched_charges = search_params.charges;
      StringList ms_runs;
      (*this->cpro_id_)[0].getPrimaryMSRunPath(ms_runs);
      String ms_runs_string = ListUtils::concatenate(ms_runs, ",");

      os << R"(<xquest_results xquest_version="OpenPepXL 1.0" date=")" << timestring <<
               R"(" author="Eugen Netz" tolerancemeasure_ms1=")" << precursor_mass_tolerance_unit  <<
               "\" tolerancemeasure_ms2=\"" << fragment_mass_tolerance_unit << "\" ms1tolerance=\"" << precursor_mass_tolerance <<
               "\" ms2tolerance=\"" << fragment_mass_tolerance << "\" xlink_ms2tolerance=\"" << fragment_mass_tolerance_xlinks <<
               "\" crosslinkername=\"" << cross_link_name << "\" xlinkermw=\"" << cross_link_mass_light <<
               "\" monolinkmw=\"" << mono_masses << "\" database=\"" << in_fasta << "\" database_dc=\"" << in_decoy_fasta <<
               R"(" xlinktypes="1111" AArequired1=")" << aarequired1 << "\" AArequired2=\"" << aarequired2 <<  "\" cp_isotopediff=\"" << cross_link_mass_iso_shift <<
               "\" enzyme_name=\"" << enzyme_name << "\" outputpath=\"" << spec_xml_name <<
               "\" missed_cleavages=\"" << missed_cleavages <<
               "\" ntermxlinkable=\"" << ntermxlinkable << "\" CID_match2ndisotope=\"1" <<
               "\" variable_mod=\"" << variable_mods << "\" fixed_mod=\"" << fixed_mods <<
               "\" decoy_prefix=\"" << decoy_prefix << "\" decoy_string=\"" << decoy_string <<
               "\" charges=\"" << searched_charges << "\" run_path=\"" << ms_runs_string <<
               R"(" nocutatxlink="1">)" << std::endl;

      String current_spectrum_light("");
      String current_spectrum_heavy("");

      for (const auto& current_pep_id : *cpep_id_)
      {
        std::vector< PeptideHit > pep_hits = current_pep_id.getHits();
        if (pep_hits.empty())
        {
          continue;
        }
        for (PeptideHit ph : pep_hits)
        {
          // TODO write the spectrum_search entry for this ph

          double precursor_mz = current_pep_id.getMZ();
          int precursor_charge = ph.getCharge();
          double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

          bool new_spectrum(false);
          new_spectrum = (ph.getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE) != current_spectrum_light) ||
              (ph.metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF) && ph.getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF) != current_spectrum_heavy);

          if (new_spectrum)
          {
            if (!current_spectrum_light.empty())
            {
              os << "</spectrum_search>" << std::endl;
            }
            current_spectrum_light = ph.getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE);
            current_spectrum_heavy = "";
            if (ph.metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF))
            {
              current_spectrum_heavy = ph.getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF);
            }

            vector<String> input_split_dir;
            vector<String> input_split;
            String base_name;
            if (!input_filename.empty())
            {
              input_filename.split(String("/"), input_split_dir);
              input_split_dir[input_split_dir.size()-1].split(String("."), input_split);
              base_name = input_split[0];
            }
            else if (ph.metaValueExists("spectrum_input_file"))
            {
              base_name = ph.getMetaValue("spectrum_input_file");
            }

            Size scan_index_light = ph.getMetaValue("spectrum_index");
            Size scan_index_heavy = scan_index_light;
            if (ph.metaValueExists("spectrum_index_heavy"))
            {
              scan_index_heavy = ph.getMetaValue("spectrum_index_heavy");
            }
            String spectrum_light_name = base_name + ".light." + scan_index_light;
            String spectrum_heavy_name = base_name + ".heavy." + scan_index_heavy;

            String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;

            String rt_scans = String(current_pep_id.getRT()) + ":";
            String mz_scans = String(precursor_mz) + ":";
            String scantype = "light_heavy";

            if (scan_index_light == scan_index_heavy)
            {
              scantype = "light";
              rt_scans += String(current_pep_id.getRT());
              mz_scans += String(precursor_mz);
            }
            else
            {
              rt_scans += ph.getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT).toString();
              mz_scans += ph.getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ).toString();
            }


            os << "<spectrum_search spectrum=\"" << spectrum_name << "\" mz_precursor=\"" << precursor_mz << "\" scantype=\"" << scantype << "\" charge_precursor=\"" << precursor_charge
                << "\" Mr_precursor=\"" << precursor_mass <<  "\" rtsecscans=\"" << rt_scans << "\" mzscans=\"" << mz_scans
                << "\" scan_index_light=\"" << scan_index_light << "\" scan_index_heavy=\"" << scan_index_heavy
                << "\" >" << std::endl;

            // TODO values missing, most of them probably unimportant:
            // mean_ionintensity = mean ion intensity of each MS2 spectrum
            // ionintensity_stdev = ion intensity spectrum_index_heavy
            // addedMass = ???
            // iontag_ncandidates = number of candidates extracted per ion tag
            // apriori_pmatch_common, apriori_pmatch_xlink = a priori probs from match-odds probability
            // ncommonions = number of common ions
            // nxlinkions = number of xlinked ions

          }
          // one of "cross-link", "mono-link" or "loop-link"
          String xltype_OPXL = ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE);
          String xltype = "monolink";


          String structure = ph.getSequence().toUnmodifiedString();
          String letter_first = structure.substr( Int(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1)), 1);

          double weight = ph.getSequence().getMonoWeight();
          int alpha_pos = Int(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1)) + 1;
          int beta_pos = 0;

          String topology = String("a") + alpha_pos;
          String id("");
          String seq_beta("");

          if (xltype_OPXL == "cross-link")
          {
            xltype = "xlink";
            beta_pos = Int(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2)) + 1;
            AASequence beta_aaseq = AASequence::fromString(ph.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE));
            structure += "-" + beta_aaseq.toUnmodifiedString();
            topology += String("-b") + beta_pos;
            weight += beta_aaseq.getMonoWeight() + double(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS));
            id = structure + "-" + topology;
            seq_beta = ph.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE);
          }
          else if (xltype_OPXL == "loop-link")
          {
            xltype = "intralink";
            beta_pos = Int(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2)) + 1;
            topology += String("-b") + beta_pos;
            String letter_second = structure.substr(beta_pos-1, 1);
            id = structure + String("-") + letter_first + alpha_pos + String("-") + letter_second + beta_pos;
            weight += cross_link_mass_light;
          }
          else // mono-link
          {
            if (ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD).toString().hasPrefix("unknown"))
            {
              weight += double(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS));
            }
            id = structure + String("-") + letter_first + alpha_pos + String("-") + Int(double(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS)));
          }

          // Precursor error calculation, rel_error is read from the metaValue for consistency, but an absolute error is also used in the xQuest format
          // use the formula, if the MetaValue is unavailable
          double theo_mz = (weight + (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / static_cast<double>(precursor_charge);
          double error = precursor_mz - theo_mz;
          double rel_error = 0;
          if (ph.metaValueExists(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM))
          {
            rel_error = double(ph.getMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM));
          }
          else
          {
            rel_error = (error / theo_mz) / 1e-6;
          }

          // Protein accessions
          String prot_alpha = ph.getPeptideEvidences()[0].getProteinAccession();
          if (ph.getPeptideEvidences().size() > 1)
          {
            for (Size i = 1; i < ph.getPeptideEvidences().size(); ++i)
            {
              prot_alpha = prot_alpha + "," + ph.getPeptideEvidences()[i].getProteinAccession();
            }
          }

          String prot_beta = "";
          if (ph.metaValueExists(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS) && ph.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS) != "-")
          {
            prot_beta = ph.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS);
          }

          String xlinkposition = String(alpha_pos);
          if (beta_pos > 0)
          {
            xlinkposition += "," + String(beta_pos);
          }

          os << "<search_hit search_hit_rank=\"" << ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK).toString() << "\" id=\"" << id << "\" type=\"" << xltype << "\" structure=\"" << structure << "\" seq1=\"" << ph.getSequence().toString() << "\" seq2=\"" << seq_beta
                << "\" prot1=\"" << prot_alpha << "\" prot2=\"" << prot_beta << "\" topology=\"" << topology << "\" xlinkposition=\"" << xlinkposition
                << "\" Mr=\"" << weight << "\" mz=\"" << theo_mz << "\" charge=\"" << precursor_charge << "\" xlinkermass=\"" << ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS).toString()
                << "\" measured_mass=\"" << precursor_mass << "\" error=\"" << error << "\" error_rel=\"" << rel_error
                << "\" xlinkions_matched=\"" << (Int(ph.getMetaValue("matched_xlink_alpha")) + Int(ph.getMetaValue("matched_xlink_beta"))) << "\" backboneions_matched=\"" << (Int(ph.getMetaValue("matched_linear_alpha")) + Int(ph.getMetaValue("matched_linear_beta")))
                << "\" xcorrx=\"" << ph.getMetaValue("OpenPepXL:xcorr xlink").toString() << "\" xcorrb=\"" << ph.getMetaValue("OpenPepXL:xcorr common").toString() << "\" match_odds=\"" << ph.getMetaValue("OpenPepXL:match-odds").toString() << "\" prescore=\"" << ph.getMetaValue("OpenPepXL:prescore").toString()
                << "\" num_of_matched_ions_alpha=\"" << (Int(ph.getMetaValue("matched_linear_alpha")) + Int(ph.getMetaValue("matched_xlink_alpha")))
                << "\" num_of_matched_ions_beta=\"" << (Int(ph.getMetaValue("matched_xlink_beta")) + Int(ph.getMetaValue("matched_linear_beta")))
                << "\" num_of_matched_common_ions_alpha=\"" << ph.getMetaValue("matched_linear_alpha").toString() << "\" num_of_matched_common_ions_beta=\"" << ph.getMetaValue("matched_linear_beta").toString()
                << "\" num_of_matched_xlink_ions_alpha=\"" << ph.getMetaValue("matched_xlink_alpha").toString() << "\" num_of_matched_xlink_ions_beta=\"" << ph.getMetaValue("matched_xlink_beta").toString()
                << "\" TIC=\"" << ph.getMetaValue("OpenPepXL:TIC").toString() << "\" wTIC=\"" << ph.getMetaValue("OpenPepXL:wTIC").toString() << "\" intsum=\"" << ph.getMetaValue("OpenPepXL:intsum").toString();

          if (ph.metaValueExists("XFDR:used_for_FDR"))
          {
            os << "\" xprophet_f=\"" << ph.getMetaValue("XFDR:used_for_FDR");
          }
          if (ph.metaValueExists("XFDR:fdr_type"))
          {
            os << "\" fdr_type=\"" << ph.getMetaValue("XFDR:fdr_type");
          }
          if (ph.metaValueExists(Constants::UserParam::XFDR_FDR))
          {
            os << "\" fdr=\"" << ph.getMetaValue(Constants::UserParam::XFDR_FDR);
          }

          // remove MetaValues, that were already used and written out with a different key.
          ph.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS);
          ph.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK);
          ph.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1);
          ph.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE);
          if (ph.metaValueExists(Constants::UserParam::OPENPEPXL_XL_POS2))
          {
            ph.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2);
          }
          ph.removeMetaValue("matched_xlink_alpha");
          ph.removeMetaValue("matched_linear_alpha");
          ph.removeMetaValue("matched_xlink_beta");
          ph.removeMetaValue("matched_linear_beta");
          ph.removeMetaValue("OpenPepXL:xcorr xlink");
          ph.removeMetaValue("OpenPepXL:xcorr common");
          ph.removeMetaValue("OpenPepXL:match-odds");
          ph.removeMetaValue("OpenPepXL:prescore");
          ph.removeMetaValue("OpenPepXL:TIC");
          ph.removeMetaValue("OpenPepXL:wTIC");
          ph.removeMetaValue("OpenPepXL:intsum");
          ph.removeMetaValue("spectrum_reference");
          ph.removeMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF);
          ph.removeMetaValue("spectrum_index");
          ph.removeMetaValue("spectrum_index_heavy");
          ph.removeMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT);
          ph.removeMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ);
          ph.removeMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM);
          ph.removeMetaValue(Constants::UserParam::XFDR_FDR);

          // also remove MetaValues, that we do not need in xquestXML
          ph.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD);
          // ph.removeMetaValue("xl_chain");

          // these metaValues can be present, e.g. if the data came from loading a xquest.xml file
          // since they are already generated by other methods, they should not be duplicated in the output
          ph.removeMetaValue("prot1");
          ph.removeMetaValue("prot2");
          ph.removeMetaValue("OpenPepXL:id");
          ph.removeMetaValue("OpenPepXL:percTIC");
          ph.removeMetaValue("OpenPepXL:score");
          ph.removeMetaValue("OpenPepXL:structure");
          ph.removeMetaValue("OpenPepXL:xlinkermass");

          // automate writing out any additional MetaValues
          std::vector<String> keys;
          ph.getKeys(keys);

          for (const String& key : keys)
          {
            os << "\" " << key << "=\"" << ph.getMetaValue(key).toString();
          }
          // score, end of the line and closing tag for this hit
          os << "\" annotated_spec=\"" << "" << "\" score=\"" << ph.getScore() << "\" >\n</search_hit>" << std::endl;

          // TODO values missing, most of them probably unimportant:
          // weighted_matchodds_mean = a weighted version of match-odds?
          // weighted_matchodds_sum
          // match_error_mean = is this per peak error?
          // match_error_stdev = is this per peak error?
          // prescore_alpha, prescore_beta
          // match_odds_alphacommon, match_odds_betacommon, match_odds_alphaxlink, match_odds_betaxlink
          // xcorrall = xcorr for the whole combined theo spectrum?
          // TIC_alpha, TIC_beta
          // apriori_match_probs
          // apriori_match_probs_log
        }
      }
      os << "</spectrum_search>" << std::endl;
      os << "</xquest_results>" << std::endl;
    }
} // namespace OpenMS   // namespace Internal
