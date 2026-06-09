// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzIdentMLHandler.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CHEMISTRY/CrossLinksDB.h>

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

using namespace std;

namespace OpenMS::Internal
  {

    void IdentificationHit::setId(const std::string& id) noexcept 
    {
      id_ = id;
    }

    const std::string& IdentificationHit::getId() const noexcept 
    {
      return id_;
    }

    void IdentificationHit::setCharge(int charge) noexcept 
    {
      charge_ = charge;
    }

    int IdentificationHit::getCharge() const noexcept 
    {
      return charge_;
    }

    void IdentificationHit::setCalculatedMassToCharge(double mz) noexcept 
    {
      calculated_mass_to_charge_ = mz;
    }

    double IdentificationHit::getCalculatedMassToCharge() const noexcept 
    {
      return calculated_mass_to_charge_;
    }

    void IdentificationHit::setExperimentalMassToCharge(double mz) noexcept 
    {
      experimental_mass_to_charge_ = mz;
    }

    double IdentificationHit::getExperimentalMassToCharge() const noexcept 
    {
      return experimental_mass_to_charge_;
    }

    void IdentificationHit::setName(const std::string& name) noexcept 
    {
      name_ = name;
    }

    const std::string& IdentificationHit::getName() const noexcept 
    {
      return name_;
    }

    void IdentificationHit::setPassThreshold(bool pass) noexcept 
    {
      pass_threshold_ = pass;
    }

    bool IdentificationHit::getPassThreshold() const noexcept 
    {
      return pass_threshold_;
    }

    void IdentificationHit::setRank(int rank) noexcept 
    {
      rank_ = rank;
    }

    int IdentificationHit::getRank() const noexcept 
    {
      return rank_;
    }

    bool IdentificationHit::operator==(const IdentificationHit& rhs) const noexcept 
    {
      return MetaInfoInterface::operator==(rhs)
          && id_ == rhs.id_
          && charge_ == rhs.charge_
          && calculated_mass_to_charge_ == rhs.calculated_mass_to_charge_
          && experimental_mass_to_charge_ == rhs.experimental_mass_to_charge_
          && name_ == rhs.name_
          && pass_threshold_ == rhs.pass_threshold_
          && rank_ == rhs.rank_;
    }

    bool IdentificationHit::operator!=(const IdentificationHit& rhs) const noexcept 
    {
      return !(*this == rhs);
    }

  SpectrumIdentification::~SpectrumIdentification() = default;

  // Equality operator
  bool SpectrumIdentification::operator==(const SpectrumIdentification & rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && hits_ == rhs.hits_;
  }

  // Inequality operator
  bool SpectrumIdentification::operator!=(const SpectrumIdentification & rhs) const
  {
    return !(*this == rhs);
  }

  void SpectrumIdentification::setHits(const vector<IdentificationHit> & hits)
  {
    hits_ = hits;
  }

  void SpectrumIdentification::addHit(const IdentificationHit & hit)
  {
    hits_.push_back(hit);
  }

  const vector<IdentificationHit> & SpectrumIdentification::getHits() const
  {
    return hits_;
  }

  Identification::~Identification() = default;

  // Equality operator
  bool Identification::operator==(const Identification & rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && creation_date_ == rhs.creation_date_
           && spectrum_identifications_ == rhs.spectrum_identifications_;
  }

  // Inequality operator
  bool Identification::operator!=(const Identification & rhs) const
  {
    return !(*this == rhs);
  }

  void Identification::setCreationDate(const DateTime & date)
  {
    creation_date_ = date;
  }

  const DateTime & Identification::getCreationDate() const
  {
    return creation_date_;
  }

  void Identification::setSpectrumIdentifications(const vector<SpectrumIdentification> & ids)
  {
    spectrum_identifications_ = ids;
  }

  void Identification::addSpectrumIdentification(const SpectrumIdentification & id)
  {
    spectrum_identifications_.push_back(id);
  }

  const vector<SpectrumIdentification> & Identification::getSpectrumIdentifications() const
  {
    return spectrum_identifications_;
  }

    MzIdentMLHandler::MzIdentMLHandler(const std::vector<ProteinIdentification>& pro_id, const PeptideIdentificationList& pep_id, const std::string& filename, const std::string& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      //~ ms_exp_(0),
      pro_id_(nullptr),
      pep_id_(nullptr),
      cpro_id_(&pro_id),
      cpep_id_(&pep_id)
    {
      initCvCaches_();
    }

    MzIdentMLHandler::MzIdentMLHandler(std::vector<ProteinIdentification>& pro_id, PeptideIdentificationList& pep_id, const std::string& filename, const std::string& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      //~ ms_exp_(0),
      pro_id_(&pro_id),
      pep_id_(&pep_id),
      cpro_id_(nullptr),
      cpep_id_(nullptr)
    {
      initCvCaches_();
    }

    void MzIdentMLHandler::initCvCaches_()
    {
      cv_.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      unimod_.loadFromOBO("PSI-MS", File::find("/CV/unimod.obo"));
      // descendants only: the parent term MS:1001143 ("PSM-level search engine specific statistic")
      // is non-numeric and must not be serialized as a scored cvParam (it is handled via the userParam fallback)
      cv_.getAllChildTerms(peptide_result_details_, "MS:1001143");
    }

    //~ TODO create MzIdentML instances from MSExperiment which contains much of the information yet needed
    //~ MzIdentMLHandler(const PeakMap& mx, const std::string& filename, const std::string& version, const ProgressLogger& logger)
    //~ : XMLHandler(filename, version),
    //~ logger_(logger),
    //~ ms_exp_(mx),
    //~ pro_id_(0),
    //~ pepid_(0),
    //~ cpepid_(0),
    //~ cpro_id_(0)
    //~ {
    //~ cv_.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
    //~ unimod_.loadFromOBO("PSI-MS",File::find("/CV/unimod.obo"));
    //~ }

    MzIdentMLHandler::~MzIdentMLHandler()
    = default;

    void MzIdentMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      tag_ = sm_.convert(qname);
      open_tags_.push_back(tag_);

      static set<std::string> to_ignore;
      if (to_ignore.empty())
      {
        to_ignore.insert("peptideSequence");
      }

      if (to_ignore.contains(tag_))
      {
        return;
      }

      //determine parent tag
      std::string parent_tag;
      if (open_tags_.size() > 1)
      {
        parent_tag = *(open_tags_.end() - 2);
      }
      std::string parent_parent_tag;
      if (open_tags_.size() > 2)
      {
        parent_parent_tag = *(open_tags_.end() - 3);
      }

      if (tag_ == "cvParam")
      {
        static const XMLCh* s_value = xercesc::XMLString::transcode("value");
        static const XMLCh* s_unit_accession = xercesc::XMLString::transcode("unitAccession");
        static const XMLCh* s_cv_ref = xercesc::XMLString::transcode("cvRef");
        //~ static const XMLCh* s_name = xercesc::XMLString::transcode("name");
        static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");

        std::string value, unit_accession, cv_ref;
        optionalAttributeAsString_(value, attributes, s_value);
        optionalAttributeAsString_(unit_accession, attributes, s_unit_accession);
        optionalAttributeAsString_(cv_ref, attributes, s_cv_ref);
        handleCVParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_accession), /* attributeAsString_(attributes, s_name), value, */ attributes, cv_ref /*,  unit_accession */);
        return;
      }

      if (tag_ == "MzIdentML")
      {
        // TODO handle version with mzid 1.2 release
        return;
      }

      if (tag_ == "Peptide")
      {
        // start new peptide
        actual_peptide_ = AASequence();

        // name attribute (opt)
        std::string name;
        if (optionalAttributeAsString_(name, attributes, "name"))
        {
          // TODO save name in AASequence
        }

        return;
      }

      if (tag_ == "Modification")
      {
        // average mass delta attribute (opt)
        // TODO

        // location attribute (opt)
        Int mod_location = -1;
        if (optionalAttributeAsInt_(mod_location, attributes, "location"))
        {
          current_mod_location_ = mod_location;
        }
        else
        {
          current_mod_location_ = -1;
        }

        // monoisotopic mass delta attribute (opt)
        // TODO

        // residues attribute (opt)
        // TODO
        return;
      }

      if (tag_ == "SpectrumIdentificationList")
      {

        return;
      }

      if (tag_ == "SpectrumIdentificationResult")
      {

        return;
      }

      if (tag_ == "SpectrumIdentificationItem")
      {
        //  <SpectrumIdentificationItem id="SII_1_1"  calculatedMassToCharge="670.86261" chargeState="2" experimentalMassToCharge="671.9" Peptide_ref="peptide_1_1" rank="1" passThreshold="true">
        // required attributes
        current_id_hit_.setId((attributeAsString_(attributes, "id")));
        current_id_hit_.setPassThreshold(asBool_(attributeAsString_(attributes, "passThreshold")));
        int rank = attributeAsInt_(attributes, "rank");
        current_id_hit_.setRank(rank - 1); // rank starts at 1 in mzid,OpenMS 0-based

        // optional attributes
        double double_value(0);
        if (optionalAttributeAsDouble_(double_value, attributes, "calculatedMassToCharge"))
        {
          current_id_hit_.setCalculatedMassToCharge(double_value);
        }

        Int int_value(0);
        if (optionalAttributeAsInt_(int_value, attributes, "chargeState"))
        {
          current_id_hit_.setCharge(int_value);
        }

        if (optionalAttributeAsDouble_(double_value, attributes, "experimentalMassToCharge"))
        {
          current_id_hit_.setExperimentalMassToCharge(double_value);
        }

        if (optionalAttributeAsDouble_(double_value, attributes, "calculatedMassToCharge"))
        {
          current_id_hit_.setCalculatedMassToCharge(double_value);
        }

        std::string string_value;
        if (optionalAttributeAsString_(string_value, attributes, "name"))
        {
          current_id_hit_.setName(string_value);
        }

        // TODO PeptideEvidence, pf:cvParam, pf:userParam, Fragmentation

        return;
      }
      error(LOAD, "MzIdentMLHandler::startElement: Unknown element found: '" + tag_ + "' in tag '" + parent_tag + "', ignoring.");
    }

    void MzIdentMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
    {
      if (tag_ == "Customizations")
      {
        std::string customizations = sm_.convert(chars);
        // TODO write customizations to Software
        return;
      }

      if (tag_ == "seq")
      {
        std::string seq = sm_.convert(chars);
        actual_protein_.setSequence(seq);
        return;
      }

      if (tag_ == "peptideSequence")
      {
        std::string pep = sm_.convert(chars);
        actual_peptide_ = AASequence::fromString(pep);
        return;
      }

      //error(LOAD, "MzIdentMLHandler::characters: Unknown character section found: '" + tag_ + "', ignoring.");
    }

    void MzIdentMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      static set<std::string> to_ignore;
      if (to_ignore.empty())
      {
        to_ignore.insert("mzIdentML");
        to_ignore.insert("cvParam");
      }

      tag_ = sm_.convert(qname);
      open_tags_.pop_back();

      if (to_ignore.contains(tag_))
      {
        return;
      }

      if (tag_ == "DataCollection")
      {
        return;
      }

      if (tag_ == "AnalysisData")
      {
        return;
      }

      if (tag_ == "ProteinDetectionList")
      {
        return;
      }

      if (tag_ == "SpectrumIdentificationList")
      {
        return;
      }

      if (tag_ == "SpectrumIdentificationResult")
      {
        return;
      }

      if (tag_ == "SpectrumIdentificationItem")
      {
        current_spectrum_id_.addHit(current_id_hit_);
        current_id_hit_ = IdentificationHit();
        return;
      }
      error(LOAD, "MzIdentMLHandler::endElement: Unknown element found: '" + tag_ + "', ignoring.");
    }

    void MzIdentMLHandler::handleCVParam_(const std::string& /* parent_parent_tag*/, const std::string& parent_tag, const std::string& accession, /* const std::string& name, */ /* const std::string& value, */ const xercesc::Attributes& attributes, const std::string& cv_ref /* , const std::string& unit_accession */)
    {
      if (parent_tag == "Modification")
      {
        if (cv_ref == "UNIMOD")
        {
          set<const ResidueModification*> mods;
          Int loc = numeric_limits<Int>::max();
          if (optionalAttributeAsInt_(loc, attributes, "location"))
          {
            std::string uni_mod_id = StringUtils::suffix(accession, ':');
            std::string residues;
            if (optionalAttributeAsString_(residues, attributes, "residues"))
            {
              // TODO handle ambiguous/multiple residues
            }
            if (loc == 0)
            {
              ModificationsDB::getInstance()->searchModifications(mods, uni_mod_id, "", ResidueModification::N_TERM);
            }
            else if (loc == (Int)actual_peptide_.size())
            {
              ModificationsDB::getInstance()->searchModifications(mods, uni_mod_id, "", ResidueModification::C_TERM);
            }
            else
            {
              ModificationsDB::getInstance()->searchModifications(mods, uni_mod_id, residues, ResidueModification::ANYWHERE);
            }
          }
          else
          {
            warning(LOAD, "location of modification not defined!");
          }
          if (mods.empty())
          {
            std::string message =std::string("Modification '") + accession + "' is unknown.";
            throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
          }
        }
      }
    }

    void MzIdentMLHandler::writeTo(std::ostream& os)
    {
      std::string cv_ns = cv_.name();
      std::string inputs_element;
      std::map<std::string, std::string> /* peps, pepevis, */ sil_map, sil_2_date;
      std::set<std::string> sen_set, sof_set, sip_set;
      std::map<std::string, std::string> sdb_ids, sen_ids, sof_ids, sdat_ids, pep_ids;
      //std::map<std::string, std::string> pep_pairs_ppxl;
      std::map<std::string, double> pp_identifier_2_thresh;
      //std::vector< std::pair<std::string, std::string> > pepid_pairs_ppxl;

      // file type-specific definitions needed for SpectraData element:
      std::map<FileTypes::Type, std::pair<std::string, std::string> > formats_map;
      formats_map[FileTypes::MZML] = make_pair("mzML format", "mzML unique identifier");
      formats_map[FileTypes::MZXML] = make_pair("ISB mzXML format", "scan number only nativeID format");
      formats_map[FileTypes::MZDATA] = make_pair("PSI mzData format", "spectrum identifier nativeID format");
      formats_map[FileTypes::MGF] = make_pair("Mascot MGF format", "multiple peak list nativeID format");


      //TODO if constructed with a msexperiment - not yet implemented
      //~ if(ms_exp_ == 0)
      //~ {
      //~ synthesize spectrum references
      //~ }
      //~ else
      //~ {
      //~ extract peptide and proteinid from msexperiment
      //~ genereate spectrum references from msexperiment foreach peptideidentification
      //~ }

      /*---------------------------------------------------------------------
      DataCollection:
      +Inputs
      -AnalysisData collected in sidlist --> unclosed element string
      ---------------------------------------------------------------------*/
      inputs_element +=std::string("\t<Inputs>\n");
      std::string spectra_data, search_database;

      /*
      1st: iterate over proteinidentification vector
      */
      //TODO read type of crosslink reagent from settings
      bool is_ppxl = false;
      for (std::vector<ProteinIdentification>::const_iterator it = cpro_id_->begin(); it != cpro_id_->end(); ++it)
      {
        //~ collect analysissoftware in this loop - does not go into inputelement
        std::string sof_id;
        std::string sof_name =std::string(it->getSearchEngine());
        std::map<std::string, std::string>::iterator soit = sof_ids.find(sof_name);
        std::string osecv;
        if (sof_name == "OMSSA")
        {
          osecv = "OMSSA";
        }
        else if (sof_name == "Mascot")
        {
          osecv = "Mascot";
        }
        else if (sof_name == "XTandem")
        {
          osecv = "X\\!Tandem";
        }
        else if (sof_name == "SEQUEST")
        {
          osecv = "Sequest";
        }
        else if (sof_name == "MS-GF+")
        {
          osecv = "MS-GF+";
        }
        else if (sof_name == "Percolator")
        {
          osecv = "Percolator";
        }
        else if (sof_name == "OpenPepXL")
        {
          osecv = "OpenPepXL";
        }
        else if (cv_.hasTermWithName(sof_name))
        {
          osecv = sof_name;
        }
        else
        {
          osecv = "analysis software";
        }

        if (soit == sof_ids.end())
        {
          sof_id = "SOF_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
          //~ TODO consider not only searchengine but also version!
          std::string sost =std::string("\t<AnalysisSoftware version=\"") + std::string(it->getSearchEngineVersion()) + "\" name=\"" + sof_name +  "\" id=\"" + sof_id + "\">\n" + std::string("\t\t<SoftwareName>\n");
          sost += "\t\t\t" + cv_.getTermByName(osecv).toXMLString(cv_ns);
          sost +=std::string("\n\t\t</SoftwareName>\n\t</AnalysisSoftware>\n");
          sof_set.insert(sost);
          sof_ids.insert(make_pair(sof_name, sof_id));
        }
        else
        {
          sof_id = soit->second;
        }

        if (it->metaValueExists("is_cross_linking_experiment") ||
                (it->metaValueExists("SpectrumIdentificationProtocol") &&
                it->getMetaValue("SpectrumIdentificationProtocol") == "MS:1002494"))
        {
          is_ppxl = true;  //needed as incoming data is structured differently and output deviates as well for ppxl
          // ppxl is like (1PeptideIdentification, 1-2 PeptideHits) but there might be more PeptideIdentifications for one spectrum
        }

        std::string thcv;
        pp_identifier_2_thresh.insert(make_pair(it->getIdentifier(),it->getSignificanceThreshold()));
        if (it->getSignificanceThreshold() != 0.0)
        {
          thcv = cv_.getTermByName("PSM-level statistical threshold").toXMLString(cv_ns,StringUtils::toStr(it->getSignificanceThreshold()));
        }
        else
        {
          thcv = cv_.getTermByName("no threshold").toXMLString(cv_ns);
        }
        // TODO add other software than searchengine for evidence trace

        // get a map from identifier to match OpenMS Protein/PeptideIdentification match string;
        std::string sil_id =  "SIL_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
        pp_identifier_2_sil_.insert(make_pair(it->getIdentifier(), sil_id));

        //~ collect SpectrumIdentificationProtocol for analysisprotocol in this loop - does not go into inputelement
        std::string sip_id = "SIP_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
        sil_2_sip_.insert(make_pair(sil_id, sip_id));

        std::string sip = "\t<SpectrumIdentificationProtocol id=\"" + std::string(sip_id) + "\" analysisSoftware_ref=\"" + std::string(sof_id) + "\">\n";
        sip += "\t\t<SearchType>\n\t\t\t" + cv_.getTermByName("ms-ms search").toXMLString(cv_ns) + "\n\t\t</SearchType>\n";
        sip += "\t\t<AdditionalSearchParams>\n";
        if (is_ppxl)
        {
          sip += "\n\t\t\t" + cv_.getTermByName("crosslinking search").toXMLString(cv_ns) + "\n";
        }
        //remove MS:1001029 written if present in <SearchDatabase> as of SearchDatabase_may rule
        ProteinIdentification::SearchParameters search_params = it->getSearchParameters();
        search_params.removeMetaValue("MS:1001029");
        writeMetaInfos_(sip, search_params, 3);
        sip +=std::string(3, '\t') + R"(<userParam name="charges" unitName="xsd:string" value=")" + search_params.charges + "\"/>\n";
//        sip +=std::string(3, '\t') + "<userParam name=\"" + "missed_cleavages" + "\" unitName=\"" + "xsd:integer" + "\" value=\"" + StringUtils::toStr(it->getSearchParameters().missed_cleavages) + "\"/>" + "\n";
        sip += "\t\t</AdditionalSearchParams>\n";
        // modifications:
        if (search_params.fixed_modifications.empty() &&
            search_params.variable_modifications.empty()
            && (!is_ppxl)) // TODO some OpenPepXL modifications are not covered by the unimod.obo and cause problems in the search_params
        {
          // no modifications used or are they just missing from the parameters?
          ModificationDefinitionsSet mod_defs;
          mod_defs.inferFromPeptides(*cpep_id_);
          mod_defs.getModificationNames(search_params.fixed_modifications,
                                        search_params.variable_modifications);
        }
        if (!search_params.fixed_modifications.empty() ||
            !search_params.variable_modifications.empty())
        {
          sip += "\t\t<ModificationParams>\n";
          writeModParam_(sip, search_params.fixed_modifications, true, 2);
          writeModParam_(sip, search_params.variable_modifications, false, 2);
          sip += "\t\t</ModificationParams>\n";
        }
        writeEnzyme_(sip, search_params.digestion_enzyme, search_params.missed_cleavages, 2);
        // TODO MassTable section
        sip +=std::string("\t\t<FragmentTolerance>\n");
        std::string unit_str = R"(unitCvRef="UO" unitName="dalton" unitAccession="UO:0000221")";
        if (search_params.fragment_mass_tolerance_ppm)
        {
          unit_str = R"(unitCvRef="UO" unitName="parts per million" unitAccession="UO:0000169")";
        }
        sip +=std::string(3, '\t') + R"(<cvParam accession="MS:1001412" name="search tolerance plus value" )" + unit_str + R"( cvRef="PSI-MS" value=")" + StringUtils::toStr(search_params.fragment_mass_tolerance) + "\"/>\n";
        sip +=std::string(3, '\t') + R"(<cvParam accession="MS:1001413" name="search tolerance minus value" )" + unit_str + R"( cvRef="PSI-MS" value=")" + StringUtils::toStr(search_params.fragment_mass_tolerance) + "\"/>\n";
        sip +=std::string("\t\t</FragmentTolerance>\n");
        sip +=std::string("\t\t<ParentTolerance>\n");
        unit_str = R"(unitCvRef="UO" unitName="dalton" unitAccession="UO:0000221")";
        if (search_params.precursor_mass_tolerance_ppm)
        {
          unit_str = R"(unitCvRef="UO" unitName="parts per million" unitAccession="UO:0000169")";
        }
        sip +=std::string(3, '\t') + R"(<cvParam accession="MS:1001412" name="search tolerance plus value" )" + unit_str + R"( cvRef="PSI-MS" value=")" + StringUtils::toStr(search_params.precursor_mass_tolerance) + "\"/>\n";
        sip +=std::string(3, '\t') + R"(<cvParam accession="MS:1001413" name="search tolerance minus value" )" + unit_str + R"( cvRef="PSI-MS" value=")" + StringUtils::toStr(search_params.precursor_mass_tolerance) + "\"/>\n";
        sip +=std::string("\t\t</ParentTolerance>\n");
        sip +=std::string("\t\t<Threshold>\n\t\t\t") + thcv + "\n";
        sip +=std::string("\t\t</Threshold>\n");
        sip +=std::string("\t</SpectrumIdentificationProtocol>\n");
        sip_set.insert(sip);
        
        // empty date would lead to XML schema validation error:
        DateTime date_time = it->getDateTime();
        if (!date_time.isValid()) 
        { 
          date_time = DateTime::now(); 
        }
        sil_2_date.insert(make_pair(sil_id,std::string(date_time.getDate() + "T" + date_time.getTime())));

        //~ collect SpectraData element for each ProteinIdentification
        std::string sdat_id;
        StringList sdat_files;
        std::string sdat_file(StringUtils::toStr(it->getMetaValue("spectra_data")));

        if (sdat_file.empty())
        {
          sdat_file =std::string("UNKNOWN");
        }
        else
        {
          sdat_file = trimOpenMSfileURI(sdat_file);
        }
        std::map<std::string, std::string>::iterator sdit = sdat_ids.find(sdat_file); //this part is strongly connected to AnalysisCollection write part
        if (sdit == sdat_ids.end())
        {
          sdat_id = "SDAT_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());

          FileTypes::Type type = FileHandler::getTypeByFileName(sdat_file);
          if (!formats_map.contains(type)) type = FileTypes::MZML; // default

          //xml
          spectra_data +=std::string("\t\t<SpectraData location=\"") + sdat_file + "\" id=\"" + sdat_id + "\">";
          spectra_data +=std::string("\n\t\t\t<FileFormat>\n");
          spectra_data +=std::string(4, '\t') + cv_.getTermByName(formats_map[type].first).toXMLString(cv_ns);
          spectra_data +=std::string("\n\t\t\t</FileFormat>\n\t\t\t<SpectrumIDFormat>\n");
          spectra_data +=std::string(4, '\t') + cv_.getTermByName(formats_map[type].second).toXMLString(cv_ns);
          spectra_data +=std::string("\n\t\t\t</SpectrumIDFormat>\n\t\t</SpectraData>\n");

          sdat_ids.insert(make_pair(sdat_file, sdat_id));
          ph_2_sdat_.insert(make_pair(it->getIdentifier(), sdat_id));
        }
        else
        {
          sdat_id = sdit->second;
        }
        sil_2_sdat_.insert(make_pair(sil_id,  sdat_id));

        //~ collect SearchDatabase element for each ProteinIdentification
        std::string sdb_id;
        std::string sdb_file(search_params.db); //TODO @mths for several IdentificationRuns this must be something else, otherwise for two of the same db just one will be created
        std::map<std::string, std::string>::iterator dbit = sdb_ids.find(sdb_file);
        if (dbit == sdb_ids.end())
        {
          sdb_id = "SDB_"+ StringUtils::toStr(UniqueIdGenerator::getUniqueId());

          search_database +=std::string("\t\t<SearchDatabase ");
          search_database +=std::string("location=\"") + sdb_file + "\" ";
          if (!std::string(search_params.db_version).empty())
          {
            search_database +=std::string("version=\"") + std::string(search_params.db_version) + "\" ";
          }
          search_database +=std::string("id=\"") + sdb_id + "\">\n\t\t\t<FileFormat>\n";
          //TODO Searchdb file format type cvParam handling
          search_database +=std::string(4, '\t') + cv_.getTermByName("FASTA format").toXMLString(cv_ns);
          search_database +=std::string("\n\t\t\t</FileFormat>\n\t\t\t<DatabaseName>\n\t\t\t\t<userParam name=\"") + sdb_file + "\"/>\n\t\t\t</DatabaseName>\n";
          // "MS:1001029" was removed from the "search_params" copy!
          if (it->getSearchParameters().metaValueExists("MS:1001029"))
          {
            search_database +=std::string(3, '\t') + cv_.getTerm("MS:1001029").toXMLString(cv_ns, it->getSearchParameters().getMetaValue("MS:1001029")) + std::string("\n");
          }
          search_database += "\t\t</SearchDatabase>\n";

          sdb_ids.insert(make_pair(sdb_file, sdb_id));
        }
        else
        {
          sdb_id = dbit->second;
        }
        sil_2_sdb_.insert(make_pair(sil_id, sdb_id));


        for (std::vector<ProteinHit>::const_iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
        {
          std::string enid;
          std::map<std::string, std::string>::iterator enit = sen_ids.find(std::string(jt->getAccession()));
          if (enit == sen_ids.end())
          {
            std::string entry;
            enid = "PROT_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId()); //TODO IDs from metadata or where its stored at read in;
            std::string enst(jt->getAccession());

            entry += "\t<DBSequence accession=\"" + enst + "\" ";
            entry += "searchDatabase_ref=\"" + sdb_id + "\" ";
            std::string s =std::string(jt->getSequence());
            if (!s.empty())
            {
              entry += "length=\"" + StringUtils::toStr(jt->getSequence().length()) + "\" ";
            }
            entry +=std::string("id=\"") + std::string(enid) + "\">\n";
            if (!s.empty())
            {
              entry += "\t<Seq>" + s + "</Seq>\n";
            }
            entry += "\t\t" + cv_.getTermByName("protein description").toXMLString(cv_ns, enst);
            entry += "\n\t</DBSequence>\n";

            sen_ids.insert(std::pair<std::string, std::string>(enst, enid));
            sen_set.insert(entry);

          }
          else
          {
            enid = enit->second;
          }
        }
      }
      inputs_element += search_database;
      inputs_element += spectra_data;
      inputs_element += "\t</Inputs>\n";

      /*
      2nd: iterate over peptideidentification vector
      */
      //TODO ppxl - write here "MS:1002511" Cross-linked spectrum identification item linking the other spectrum
      //          PeptideIdentification represents xl pair.
      //          PeptideHit score_type is then the final score of xQuest_cpp.
      //          top5 ids -> 5 PeptideIdentification for one (pair) spectra. SIR with 5 entries and ranks
      std::map<std::string, std::string> ppxl_specref_2_element; //where the SII will get added for one spectrum reference
      std::map<std::string, std::vector<std::string> > pep_evis; //maps the sequence to the corresponding evidence elements for the next scope
      for (PeptideIdentificationList::const_iterator it = cpep_id_->begin(); it != cpep_id_->end(); ++it)
      {
        std::string emz = StringUtils::toStr(it->getMZ());
        const double rt = it->getRT();
        std::string ert = rt == rt ? StringUtils::toStr(rt) : "nan";

        std::string sid = StringUtils::toStr(it->getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE));
        if (sid.empty())
        {
          sid =StringUtils::toStr(it->getMetaValue("spectrum_id"));
          if (sid.empty())
          {
              if (it->getMZ() != it->getMZ())
            {
              emz = "nan";
              OPENMS_LOG_WARN << "Found no spectrum reference and no m/z position of identified spectrum! You are probably converting from an old format with insufficient data provision. Setting 'nan' - downstream applications might fail unless you set the references right.\n";
            }
            if (it->getRT() != it->getRT())
            {
              ert = "nan";
              OPENMS_LOG_WARN << "Found no spectrum reference and no RT position of identified spectrum! You are probably converting from an old format with insufficient data provision. Setting 'nan' - downstream applications might fail unless you set the references right.\n";
            }
            sid =std::string("MZ:") + emz + std::string("@RT:") + ert;
          }
        }
        if (is_ppxl && it->metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF))
        {
          sid.append("," + StringUtils::toStr(it->getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF)));
        }
        std::string sidres;
        std::string sir = "SIR_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
        std::string sdr = sdat_ids.begin()->second;
        std::map<std::string, std::string>::iterator pfo = ph_2_sdat_.find(it->getIdentifier());
        if (pfo != ph_2_sdat_.end())
        {
          sdr = pfo->second;
        }
        else
        {
          OPENMS_LOG_WARN << "Falling back to referencing first spectrum file given because file or identifier could not be mapped.\n";
        }

        sidres +=std::string("\t\t\t<SpectrumIdentificationResult spectraData_ref=\"")
        //multi identification runs lookup from file_origin here
                + sdr + "\" spectrumID=\""
                + sid + "\" id=\"" + sir + "\">\n";

        if (is_ppxl)
        {
          if (!ppxl_specref_2_element.contains(sid))
          {
            ppxl_specref_2_element[sid] = sidres;
          }
        }

        //map.begin access ok here because make sure at least one "UNKOWN" element is in the sdats_ids map

        double ppxl_crosslink_mass(0);
        if (is_ppxl)
        {
          ProteinIdentification::SearchParameters search_params = cpro_id_->front().getSearchParameters();
          // use a default so a missing/empty cross_link:mass meta value does not throw a ConversionError
          // (e.g. on a store after a load where the search parameter was not preserved)
          ppxl_crosslink_mass = StringUtils::toDouble(StringUtils::toStr(search_params.getMetaValue("cross_link:mass", "0")));
        }

        for (std::vector<PeptideHit>::const_iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
        {
          if (!is_ppxl)
          {
            MzIdentMLHandler::writePeptideHit(*jt, it, pep_ids, cv_ns, sen_set, sen_ids, pep_evis, pp_identifier_2_thresh, sidres);
          }
          else
          {
            std::string ppxl_linkid = StringUtils::toStr(UniqueIdGenerator::getUniqueId());
            MzIdentMLHandler::writeXLMSPeptideHit(*jt, it, ppxl_linkid, pep_ids, cv_ns, sen_set, sen_ids, pep_evis, pp_identifier_2_thresh, ppxl_crosslink_mass, ppxl_specref_2_element, sid, true);
            // XL-MS IDs from OpenPepXL can have two Peptides and SpectrumIdentifications, but with practically the same data except for the sequence and its modifications
            if (jt->metaValueExists(Constants::UserParam::OPENPEPXL_XL_TYPE) && jt->getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "cross-link")
            {
              MzIdentMLHandler::writeXLMSPeptideHit(*jt, it, ppxl_linkid, pep_ids, cv_ns, sen_set, sen_ids, pep_evis, pp_identifier_2_thresh, ppxl_crosslink_mass, ppxl_specref_2_element, sid, false);
            }
          }
        }
        if (!ert.empty() && ert != "nan" && ert != "NaN" && !is_ppxl)
        {
          DataValue rtcv(ert);
          rtcv.setUnit(10); // id: UO:0000010 name: second
          rtcv.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
          sidres +=  "\t\t\t\t" + cv_.getTermByName("retention time").toXMLString(cv_ns, rtcv) + "\n";
        }
        if (!is_ppxl)
        {
          sidres += "\t\t\t</SpectrumIdentificationResult>\n";
          std::map<std::string, std::string>::const_iterator ps_it = pp_identifier_2_sil_.find(it->getIdentifier());
          if (ps_it != pp_identifier_2_sil_.end())
          {
            std::map<std::string, std::string>::iterator sil_it = sil_map.find(ps_it->second);
            if (sil_it != sil_map.end())
            {
              sil_it->second.append(sidres);
            }
            else
            {
              sil_map.insert(make_pair(ps_it->second,sidres));
            }
          }
          else
          {
            //encountered a PeptideIdentification which is not linked to any ProteinIdentification
            OPENMS_LOG_ERROR << "encountered a PeptideIdentification which is not linked to any ProteinIdentification\n";
          }
        }

      }
      // ppxl - write spectrumidentificationresult closing tags!
      if (is_ppxl)
      {
        for (std::map<std::string, std::string>::iterator it = ppxl_specref_2_element.begin(); it != ppxl_specref_2_element.end(); ++it)
        {
          it->second += "\t\t\t</SpectrumIdentificationResult>\n";
          std::map<std::string, std::string>::const_iterator ps_it = pp_identifier_2_sil_.begin();

          if (ps_it != pp_identifier_2_sil_.end())
          {
            std::map<std::string, std::string>::iterator sil_it = sil_map.find(ps_it->second);
            if (sil_it != sil_map.end())
            {
              sil_it->second.append(it->second);
            }
            else
            {
              sil_map.insert(make_pair(ps_it->second,it->second));
            }
          }
          else
          {
            //encountered a PeptideIdentification which is not linked to any ProteinIdentification
            OPENMS_LOG_ERROR << "encountered a PeptideIdentification crosslink information which is not linked to any ProteinIdentification\n";
          }
        }
      }

      //--------------------------------------------------------------------------------------------
      // XML header
      //--------------------------------------------------------------------------------------------
      std::string v_s = "1.3.0";
      // namespace uses only major.minor (e.g. "1.3" from "1.3.0"); derive it by dropping the patch component
      std::string v_short = StringUtils::prefix(v_s, v_s.rfind('.'));
      os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
         << "<MzIdentML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
         << "\txsi:schemaLocation=\"http://psidev.info/psi/pi/mzIdentML/" << v_short << " "
         << "https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML" << v_s << ".xsd\"\n"
         << "\txmlns=\"http://psidev.info/psi/pi/mzIdentML/" << v_short << "\"\n"
         << "\tversion=\"" << v_s << "\"\n";
      os << "\tid=\"OpenMS_" << StringUtils::toStr(UniqueIdGenerator::getUniqueId()) << "\"\n"
         << "\tcreationDate=\"" << DateTime::now().getDate() << "T" << DateTime::now().getTime() << "\">\n";

      //--------------------------------------------------------------------------------------------
      // CV list
      //--------------------------------------------------------------------------------------------
      os << "<cvList>\n"
         << "\t<cv id=\"PSI-MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Vocabularies\" "
         << "uri=\"http://purl.obolibrary.org/obo/ms/psi-ms.obo\" "
         << "version=\"3.15.0\"></cv>\n "
         << "\t<cv id=\"UNIMOD\" fullName=\"UNIMOD\" uri=\"http://www.unimod.org/obo/unimod.obo\"></cv>\n"
         << "\t<cv id=\"UO\"     fullName=\"UNIT-ONTOLOGY\" "
         << "uri=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"></cv>\n";
      if (is_ppxl)
      {
          os << "\t<cv id=\"XLMOD\" fullName=\"PSI cross-link modifications\" "
             << "uri=\"https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD-1.0.0.obo\"></cv>\n";
      }
      os << "</cvList>\n";

      //--------------------------------------------------------------------------------------------
      // AnalysisSoftwareList
      //--------------------------------------------------------------------------------------------
      os << "<AnalysisSoftwareList>\n";
      for (std::set<std::string>::const_iterator sof = sof_set.begin(); sof != sof_set.end(); ++sof)
      {
        os << *sof;
      }

      std::map<std::string, std::string>::iterator soit = sof_ids.find("TOPP software");
      if (soit == sof_ids.end())
      {
        os << "\t<AnalysisSoftware version=\"OpenMS TOPP v"<< VersionInfo::getVersion() <<R"(" name="TOPP software" id=")" << std::string("SOF_") << StringUtils::toStr(UniqueIdGenerator::getUniqueId()) << "\">\n"
           << "\t\t<SoftwareName>\n\t\t\t" << cv_.getTermByName("TOPP software").toXMLString(cv_ns) << "\n\t\t</SoftwareName>\n\t</AnalysisSoftware>\n";
      }
      os << "</AnalysisSoftwareList>\n";

      //--------------------------------------------------------------------------------------------
      // SequenceCollection
      //--------------------------------------------------------------------------------------------
      os << "<SequenceCollection>\n";
      for (std::set<std::string>::const_iterator sen = sen_set.begin(); sen != sen_set.end(); ++sen)
      {
        os << *sen;
      }
      os << "</SequenceCollection>\n";

      //--------------------------------------------------------------------------------------------
      // AnalysisCollection:
      // + SpectrumIdentification
      // TODO ProteinDetection
      //--------------------------------------------------------------------------------------------
      os << "<AnalysisCollection>\n";
      for (std::map<std::string, std::string>::const_iterator pp2sil_it = pp_identifier_2_sil_.begin(); pp2sil_it != pp_identifier_2_sil_.end(); ++pp2sil_it)
      {
          std::string entry =std::string("\t<SpectrumIdentification id=\"SI_") + pp2sil_it->first + "\" spectrumIdentificationProtocol_ref=\""
                           + sil_2_sip_[pp2sil_it->second] + "\" spectrumIdentificationList_ref=\"" + pp2sil_it->second
                           + "\" activityDate=\"" + sil_2_date[pp2sil_it->second]
                           + "\">\n"
                            //if crosslink +cvparam crosslink search performed
                           + "\t\t<InputSpectra spectraData_ref=\"" + sil_2_sdat_[pp2sil_it->second] + "\"/>\n" // spd_ids.insert(std::pair<std::string, UInt64>(sdst, sdid));
                           + "\t\t<SearchDatabaseRef searchDatabase_ref=\"" + sil_2_sdb_[pp2sil_it->second] + "\"/>\n"
                           + "\t</SpectrumIdentification>\n";
          os <<   entry;
      }
      os << "</AnalysisCollection>\n";

      //--------------------------------------------------------------------------------------------
      // AnalysisProtocolCollection
      //+ SpectrumIdentificationProtocol + SearchType + Threshold
      //--------------------------------------------------------------------------------------------
      os << "<AnalysisProtocolCollection>\n";
      for (std::set<std::string>::const_iterator sip = sip_set.begin(); sip != sip_set.end(); ++sip)
      {
        os << *sip;
      }
      os << "</AnalysisProtocolCollection>\n";

      //--------------------------------------------------------------------------------------------
      // DataCollection
      //+Inputs
      //+AnalysisData
      //--------------------------------------------------------------------------------------------
      os << "<DataCollection>\n"
         << inputs_element;
      os << "\t<AnalysisData>\n";
      for (std::map<std::string, std::string>::const_iterator sil_it = sil_map.begin(); sil_it != sil_map.end(); ++sil_it)
      {
        os << "\t\t<SpectrumIdentificationList id=\"" << sil_it->first << "\">\n";
        os << "\t\t\t<FragmentationTable>\n"
           << "\t\t\t\t<Measure id=\"Measure_mz\">\n"
           << "\t\t\t\t\t<cvParam accession=\"MS:1001225\" cvRef=\"PSI-MS\" unitCvRef=\"PSI-MS\" unitName=\"m/z\" unitAccession=\"MS:1000040\" name=\"product ion m/z\"/>\n"
           << "\t\t\t\t</Measure>\n"
           << "\t\t\t\t<Measure id=\"Measure_int\">\n"
           << "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001226\" name=\"product ion intensity\" unitAccession=\"MS:1000131\" unitCvRef=\"UO\" unitName=\"number of detector counts\"/>\n"
           << "\t\t\t\t</Measure>\n"
           << "\t\t\t\t<Measure id=\"Measure_error\">\n"
           << "\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001227\" name=\"product ion m/z error\" unitAccession=\"MS:1000040\" unitCvRef=\"PSI-MS\" unitName=\"m/z\"/>\n"
           << "\t\t\t\t</Measure>\n";
        if (is_ppxl)
        {
            os << "<!-- userParam cross-link_chain will contain a list of chain type corresponding to the indexed ion [alpha|beta] -->\n";
            os << "<!-- userParam cross-link_ioncategory will contain a list of ion category corresponding to the indexed ion [xi|ci] -->\n";
        }
        os << "\t\t\t</FragmentationTable>\n";
        os << sil_it->second;
        os << "\t\t</SpectrumIdentificationList>\n";
      }
      os << "\t</AnalysisData>\n</DataCollection>\n";

      //--------------------------------------------------------------------------------------------
      // close XML header
      //--------------------------------------------------------------------------------------------
      os << "</MzIdentML>\n";
    }

    void MzIdentMLHandler::writeMetaInfos_(std::string& s, const MetaInfoInterface& meta, UInt indent) const
    {
      //TODO @mths: write those metas with their name in the cvs loaded as CVs!
      if (meta.isMetaEmpty())
      {
        return;
      }
      std::vector<std::string> keys;
      meta.getKeys(keys);

      for (Size i = 0; i != keys.size(); ++i)
      {
        if (cv_.exists(keys[i]))
        {
          ControlledVocabulary::CVTerm a = cv_.getTerm(keys[i]);
          s +=std::string(indent, '\t') + a.toXMLString("PSI-MS", StringUtils::toStr((meta.getMetaValue(keys[i])))) + "\n";
        }
        else
        {
          s +=std::string(indent, '\t') + "<userParam name=\"" + keys[i] + "\" unitName=\"";

          const DataValue& d = meta.getMetaValue(keys[i]);
          //determine type
          if (d.valueType() == DataValue::INT_VALUE)
          {
            s += "xsd:integer";
          }
          else if (d.valueType() == DataValue::DOUBLE_VALUE)
          {
            s += "xsd:double";
          }
          else //string or lists are converted to string
          {
            s += "xsd:string";
          }
          s += "\" value=\"" + StringUtils::toStr(d) + "\"/>\n";
        }
      }
    }

    void MzIdentMLHandler::writeEnzyme_(std::string& s, const DigestionEnzymeProtein& enzy, UInt miss, UInt indent) const
    {
      std::string cv_ns = cv_.name();
      s +=std::string(indent, '\t') + "<Enzymes independent=\"false\">\n";
      s +=std::string(indent + 1, '\t') + "<Enzyme missedCleavages=\"" + StringUtils::toStr(miss) + "\" id=\"" + std::string("ENZ_") + StringUtils::toStr(UniqueIdGenerator::getUniqueId()) + "\">\n";
      s +=std::string(indent + 2, '\t') + "<EnzymeName>\n";
      const std::string& enzymename = enzy.getName();
      if (cv_.hasTermWithName(enzymename))
      {
        s +=std::string(indent + 3, '\t') + cv_.getTermByName(enzymename).toXMLString(cv_ns) + "\n";
      }
      else if (enzymename == "no cleavage")
      {
        s +=std::string(indent + 3, '\t') + cv_.getTermByName("NoEnzyme").toXMLString(cv_ns) + "\n";
      }
      else
      {
        s +=std::string(indent + 3, '\t') + cv_.getTermByName("cleavage agent details").toXMLString(cv_ns) + "\n";
      }
      s +=std::string(indent + 2, '\t') + "</EnzymeName>\n";
      s +=std::string(indent + 1, '\t') + "</Enzyme>\n";
      s +=std::string(indent, '\t') + "</Enzymes>\n";
    }

    void MzIdentMLHandler::writeModParam_(std::string& s, const std::vector<std::string>& mod_names, bool fixed, UInt indent) const
    {
      std::string cv_ns = unimod_.name();
      for (std::vector<std::string>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
      {
        std::set<const ResidueModification*> mods;
        ModificationsDB::getInstance()->searchModifications(mods, *it);
        if (!mods.empty())
        {
          // @TODO: if multiple mods match, we write all of them?
          for (std::set<const ResidueModification*>::const_iterator mt = mods.begin(); mt != mods.end(); ++mt)
          {
            char origin = (*mt)->getOrigin();
            if (origin == 'X') origin = '.'; // terminal without res. spec.

            s +=std::string(indent + 1, '\t') + "<SearchModification fixedMod=\"" + (fixed ? "true" : "false") + "\" massDelta=\"" + StringUtils::toStr((*mt)->getDiffMonoMass()) + "\" residues=\"" + origin + "\">\n";

            // @TODO: handle protein C-term/N-term
            ResidueModification::TermSpecificity spec = (*mt)->getTermSpecificity();
            if ((spec == ResidueModification::C_TERM) || (spec == ResidueModification::N_TERM))
            {
              const std::string& cv_name = "modification specificity peptide " + (*mt)->getTermSpecificityName();
              s +=std::string(indent + 2, '\t') + "<SpecificityRules>\n";
              s +=std::string(indent + 3, '\t') + cv_.getTermByName(cv_name).toXMLString(cv_ns) + "\n";
              s +=std::string(indent + 2, '\t') + "</SpecificityRules>\n";
            }

            std::string ac = (*mt)->getUniModAccession();
            if (StringUtils::hasPrefix(ac, "UniMod:")) ac = "UNIMOD:" + StringUtils::suffix(ac, ':');
            if (!ac.empty())
            {
              s +=std::string(indent + 2, '\t') + unimod_.getTerm(ac).toXMLString(cv_ns) + "\n";
            }
            else
            {
              s +=std::string(indent + 2, '\t') + "<cvParam cvRef=\"MS\" accession=\"MS:1001460\" name=\"unknown modification\"/>\n";
            }
            s +=std::string(indent + 1, '\t') + "</SearchModification>\n";
          }
        }
        else
        {
          std::string message =std::string("Registered ") + (fixed ? "fixed" : "variable") + " modification '" + *it + "' is unknown and will be ignored.";
          throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
        }
      }
    }

    void MzIdentMLHandler::writeFragmentAnnotations_(std::string& s, const std::vector<PeptideHit::PeakAnnotation>& annotations, UInt indent, bool is_ppxl) const
    {
      std::map<UInt,std::map<std::string,std::vector<StringList> > > annotation_map;
      for (const PeptideHit::PeakAnnotation& pep : annotations)
      {// string coding example: [alpha|ci$y3-H2O-NH3]5+
        // static const boost::regex frag_regex("\\[(?:([\\|\\w]+)\\$)*([abcxyz])(\\d+)((?:[\\+\\-\\w])*)\\](\\d+)\\+"); // this will fetch the complete loss/gain part as one
        static const boost::regex frag_regex_tweak(R"(\[(?:([\|\w]+)\$)*([abcxyz])(\d+)(?:-(H2O|NH3))*\][(\d+)\+]*)"); // this will only fetch the last loss - and is preferred for now, as only these extra cv params are present
        std::string ionseries_index;
        std::string iontype;
        //std::string loss_gain;
        std::string loss;
        StringList extra;
        boost::smatch str_matches;
        if (boost::regex_match(pep.annotation, str_matches, frag_regex_tweak))
        {
          StringUtils::split(std::string(str_matches[1]), "|", extra);
          iontype = std::string(str_matches[2]);
          ionseries_index = std::string(str_matches[3]);
          loss = std::string(str_matches[4]);
        }
        else
        {
          // since PeakAnnotations are very flexible and not all of them fit into the limited mzid fragment structure,
          // this would happen quite often and flood the output, but we still need them for other output formats
          // TODO find ways to represent additional fragment types or filter out known incompatible types

          // OPENMS_LOG_WARN << "Well, fudge you very much, there is no matching annotation. ";
          // OPENMS_LOG_WARN << pep.annotation << '\n';
          continue;
        }
        std::string lt = "frag: " + iontype + " ion";
        if (!loss.empty())
        {
          lt += " - "+loss;
        }
        auto& charge_map = annotation_map[pep.charge];
        auto& lt_vec = charge_map[lt];
        if (lt_vec.empty())
        {
          lt_vec.resize(3);
          if (is_ppxl)
          {
            lt_vec.push_back(StringList());
            lt_vec.push_back(StringList());
          }
        }
        lt_vec[0].push_back(ionseries_index);
        lt_vec[1].emplace_back(StringUtils::toStr(pep.mz));
        lt_vec[2].emplace_back(StringUtils::toStr(pep.intensity));
        if (is_ppxl)
        {
          std::string ab = ListUtils::contains<std::string>(extra ,std::string("alpha")) ? std::string("alpha"):std::string("beta");
          std::string cx = ListUtils::contains<std::string>(extra ,std::string("ci")) ? std::string("ci"):std::string("xi");
          lt_vec[3].push_back(ab);
          lt_vec[4].push_back(cx);
        }
      }

      // stop and return, if no mzid compatible fragments were found
      if (annotation_map.empty())
      {
        return;
      }
      //double map: charge + ion type; collect in StringList: index + annotations; write:
      s +=std::string(indent, '\t') + "<Fragmentation>\n";
      for (std::map<UInt,std::map<std::string,std::vector<StringList> > >::iterator i=annotation_map.begin();
           i!=annotation_map.end(); ++i)
      {
        for (std::map<std::string,std::vector<StringList> >::iterator j=i->second.begin();
             j!= i->second.end(); ++j)
        {
          s +=std::string(indent+1, '\t') + "<IonType charge=\"" + StringUtils::toStr(i->first) +"\""
                    + " index=\"" + ListUtils::concatenate(j->second[0], " ") + "\">\n";
          s +=std::string(indent+2, '\t') + "<FragmentArray measure_ref=\"Measure_mz\""
                    + " values=\"" + ListUtils::concatenate(j->second[1], " ") + "\"/>\n";
          s +=std::string(indent+2, '\t') + "<FragmentArray measure_ref=\"Measure_int\""
                    + " values=\"" + ListUtils::concatenate(j->second[2], " ") + "\"/>\n";
          if (is_ppxl)
          {
              s +=std::string(indent+2, '\t') + "<userParam name=\"cross-link_chain\"" + " unitName=\"xsd:string\""
                        + " value=\"" + ListUtils::concatenate(j->second[3], " ") + "\"/>\n";
              s +=std::string(indent+2, '\t') + "<userParam name=\"cross-link_ioncategory\"" + " unitName=\"xsd:string\""
                        + " value=\"" + ListUtils::concatenate(j->second[4], " ") + "\"/>\n";
          }
          s +=std::string(indent+2, '\t') + cv_.getTermByName(j->first).toXMLString("PSI-MS") + "\n";
          s +=std::string(indent+1, '\t') + "</IonType>\n";
        }
      }
      s +=std::string(indent, '\t') + "</Fragmentation>\n";
//<Fragmentation>
//    <IonType charge="1" index="2 3 5 6 5 6 10">
//        <FragmentArray measure_ref="Measure_MZ" values="363.908 511.557 754.418 853.489 377.941 427.477 633.674"/>
//        <FragmentArray measure_ref="Measure_Int" values="208.52 2034.9 1098.44 239.26 3325.34 3028.33 335.63"/>
//        <FragmentArray measure_ref="Measure_Error" values="-0.255 0.326 0.101 0.104 0.285 0.287 0.369"/>
//        <cvParam accession="MS:1001118" cvRef="PSI-MS" name="param: b ion"/>
//    </IonType>
//    <IonType charge="1" index="0 1 3 4 5 6 4 10">
//        <FragmentArray measure_ref="Measure_MZ" values="175.202 246.812 474.82 587.465 686.52 814.542 294.235 652.206"/>
//        <FragmentArray measure_ref="Measure_Int" values="84.44999 90.26999 143.95 3096.84 815.34 1.15999 612.52 18.79999"/>
//        <FragmentArray measure_ref="Measure_Error" values="0.083 0.656 0.553 0.114 0.101 0.064 0.061 0.375"/>
//        <cvParam accession="MS:1001262" cvRef="PSI-MS" name="param: y ion"/>
//    </IonType>
//</Fragmentation>
    }

    std::string MzIdentMLHandler::trimOpenMSfileURI(const std::string& file) const
    {
      std::string r = file;
      if (StringUtils::hasPrefix(r, "["))
        r = StringUtils::substr(r, 1);
      if (StringUtils::hasSuffix(r, "]"))
        r = StringUtils::substr(r, 0,r.size()-1);
      StringUtils::substitute(r, "\\","/");
      return r;
    }

    void MzIdentMLHandler::writePeptideHit(const PeptideHit& hit,
                                                PeptideIdentificationList::const_iterator& it,
                                                std::map<std::string, std::string>& pep_ids,
                                                const std::string& cv_ns, std::set<std::string>& sen_set,
                                                std::map<std::string, std::string>& sen_ids,
                                                std::map<std::string, std::vector<std::string> >& pep_evis,
                                                std::map<std::string, double>& pp_identifier_2_thresh,
                                                std::string& sidres)
    {
        std::string pepid =  "PEP_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
        std::string pepi = hit.getSequence().toString();

        bool duplicate = false;
        std::map<std::string, std::string>::iterator pit;

        // avoid duplicates
        pit = pep_ids.find(pepi);
        if (pit != pep_ids.end())
        {
          duplicate = true;
        }

        if (!duplicate)
        {
          std::string p;
          //~ TODO simplify mod cv param write
          // write peptide id with conversion to universal, "human-readable" bracket string notation
          p +=std::string("\t<Peptide id=\"") + pepid + "\" name=\"" +
                hit.getSequence().toBracketString(false) + "\">\n\t\t<PeptideSequence>" + hit.getSequence().toUnmodifiedString() + std::string("</PeptideSequence>\n");
          if (hit.getSequence().isModified())
          {
            const ResidueModification* n_term_mod = hit.getSequence().getNTerminalModification();
            if (n_term_mod != nullptr)
            {
              p += "\t\t<Modification location=\"0\">\n";
              std::string acc = n_term_mod->getUniModAccession();
              p += "\t\t\t<cvParam accession=\"UNIMOD:" + StringUtils::suffix(acc, ':');
              p += "\" name=\"" + n_term_mod->getId();
              p += R"(" cvRef="UNIMOD"/>)";
              p += "\n\t\t</Modification>\n";
            }
            const ResidueModification* c_term_mod = hit.getSequence().getCTerminalModification();
            if (c_term_mod != nullptr)
            {
              p += "\t\t<Modification location=\"" + StringUtils::toStr(hit.getSequence().size()) + "\">\n";
              std::string acc = c_term_mod->getUniModAccession();
              p += "\t\t\t<cvParam accession=\"UNIMOD:" + StringUtils::suffix(acc, ':');
              p += "\" name=\"" + c_term_mod->getId();
              p += R"(" cvRef="UNIMOD"/>)";
              p += "\n\t\t</Modification>\n";
            }
            for (Size i = 0; i < hit.getSequence().size(); ++i)
            {
              const ResidueModification* mod = hit.getSequence()[i].getModification(); // "UNIMOD:" prefix??
              if (mod != nullptr)
              {
                //~ p += hit.getSequence()[i].getModification() + "\t" +  hit.getSequence()[i].getOneLetterCode()  + "\t" +  x +   "\n" ;
                p += "\t\t<Modification location=\"" + StringUtils::toStr(i + 1);
                p += "\" residues=\"" + hit.getSequence()[i].getOneLetterCode();
                std::string acc = mod->getUniModAccession();
                if (!acc.empty())
                {
                  p += "\">\n\t\t\t<cvParam accession=\"UNIMOD:" + StringUtils::suffix(acc, ':'); //TODO @all: do not exclusively use unimod ...
                  p += "\" name=\"" + mod->getId();
                  p += R"(" cvRef="UNIMOD"/>)";
                  p += "\n\t\t</Modification>\n";
                }
                else
                {
                  // We have an unknown modification, so lets write unknown
                  // and at least try to write down the delta mass.
                  if (mod->getDiffMonoMass() != 0.0)
                  {
                    double diffmass = mod->getDiffMonoMass();
                    p += "\" monoisotopicMassDelta=\"" + StringUtils::toStr(diffmass);
                  }
                  else if (mod->getMonoMass() > 0.0)
                  {
                    double diffmass = mod->getMonoMass() - hit.getSequence()[i].getMonoWeight();
                    p += "\" monoisotopicMassDelta=\"" + StringUtils::toStr(diffmass);
                  }
                  p += "\">\n\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001460\" name=\"unknown modification\"/>";
                  p += "\n\t\t</Modification>\n";
                }
              }
            }
          }
          p += "\t</Peptide>\n ";
          sen_set.insert(p);
          pep_ids.insert(std::make_pair(pepi, pepid));
        }
        else
        {
          pepid = pit->second;
        }

        std::vector<std::string> pevid_ids;
        if (!duplicate)
        {
          std::vector<PeptideEvidence> peptide_evidences = hit.getPeptideEvidences();
          // TODO idXML allows peptide hits without protein references! Fails in that case - run PeptideIndexer first
          for (std::vector<PeptideEvidence>::const_iterator pe = peptide_evidences.begin(); pe != peptide_evidences.end(); ++pe)
          {
            std::string pevid =  "PEV_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
            std::string dBSequence_ref;
            map<std::string, std::string>::const_iterator pos = sen_ids.find(pe->getProteinAccession());
            if (pos != sen_ids.end())
            {
              dBSequence_ref = pos->second;
            }
            else
            {
              OPENMS_LOG_ERROR << "Error: Missing or invalid protein reference for peptide '" << pepi << "': '" << pe->getProteinAccession() << "' - skipping." << endl;
              continue;
            }
            std::string idec;
            if (hit.metaValueExists(Constants::UserParam::TARGET_DECOY))
            {
              idec =std::string(boost::lexical_cast<std::string>(StringUtils::hasSubstring(hit.getMetaValue(Constants::UserParam::TARGET_DECOY).toString(), "decoy")));
            }

            std::string e;
            std::string nc_termini = "-";    // character for N- and C-termini as specified in mzIdentML
            e += "\t<PeptideEvidence id=\"" + pevid + "\" peptide_ref=\"" + pepid + "\" dBSequence_ref=\"" + dBSequence_ref + "\"";

            if (pe->getAAAfter() != PeptideEvidence::UNKNOWN_AA)
            {
              e += " post=\"" + (pe->getAAAfter() == PeptideEvidence::C_TERMINAL_AA ? nc_termini : StringUtils::toStr(pe->getAAAfter())) + "\"";
            }
            if (pe->getAABefore() != PeptideEvidence::UNKNOWN_AA)
            {
              e += " pre=\"" + (pe->getAABefore() == PeptideEvidence::N_TERMINAL_AA ? nc_termini : StringUtils::toStr(pe->getAABefore())) + "\"";
            }
            if (pe->getStart() != PeptideEvidence::UNKNOWN_POSITION)
            {
              e += " start=\"" + StringUtils::toStr(pe->getStart() + 1) + "\"";
            }
            else if (hit.metaValueExists("start"))
            {
              e += " start=\"" + StringUtils::toStr( int(hit.getMetaValue("start")) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no start position of peptide hit in protein sequence.\n";
            }
            if (pe->getEnd() != PeptideEvidence::UNKNOWN_POSITION)
            {
              e += " end=\"" + StringUtils::toStr(pe->getEnd() + 1) + "\"";
            }
            else if (hit.metaValueExists("end"))
            {
              e += " end=\"" + StringUtils::toStr( int(hit.getMetaValue("end")) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no end position of peptide hit in protein sequence.\n";
            }
            if (!idec.empty())
            {
              e += " isDecoy=\"" + std::string(idec)+ "\"";
            }
            e += "/>\n";
            sen_set.insert(e);
            pevid_ids.push_back(pevid);
          }
          pep_evis.insert(make_pair(pepi, pevid_ids));
        }
        else
        {
          pevid_ids =  pep_evis[pepi];
        }

        std::string cmz = StringUtils::toStr(hit.getSequence().getMZ(hit.getCharge())); //calculatedMassToCharge
        std::string r = StringUtils::toStr(hit.getRank() + 1); //rank
        std::string sc = StringUtils::toStr(hit.getScore());

        if (sc.empty())
        {
          sc = "NA";
          OPENMS_LOG_WARN << "No score assigned to this PSM: " /*<< hit.getSequence().toString()*/ << '\n';
        }
        std::string c = StringUtils::toStr(hit.getCharge()); //charge

        std::string pte;
        if (hit.metaValueExists("pass_threshold"))
        {
          pte = boost::lexical_cast<std::string>(hit.getMetaValue("pass_threshold"));
        }
        else
        {
          auto thresh_it = pp_identifier_2_thresh.find(it->getIdentifier());
          if (thresh_it != pp_identifier_2_thresh.end() && thresh_it->second != 0.0)
          {
            double th = thresh_it->second;
            //threshold was 'set' in proteinIdentification (!= default value of member, now check pass
            pte = boost::lexical_cast<std::string>(it->isHigherScoreBetter() ? hit.getScore() > th : hit.getScore() < th); //passThreshold-eval
          }
          else
          {
            pte = "1";
          }
        }

        //write SpectrumIdentificationItem elements
        std::string emz = StringUtils::toStr(it->getMZ());
        std::string sii = "SII_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
        std::string sii_tmp;
        sii_tmp +=std::string("\t\t\t\t<SpectrumIdentificationItem passThreshold=\"")
                + pte + "\" rank=\"" + r + "\" peptide_ref=\""
                + pepid + "\" calculatedMassToCharge=\"" + cmz
                + "\" experimentalMassToCharge=\"" + emz
                + "\" chargeState=\"" + c +  "\" id=\""
                + sii + "\">\n";

        if (pevid_ids.empty())
        {
          OPENMS_LOG_WARN << "PSM without peptide evidence registered in the given search database found. This will cause an invalid mzIdentML file (which OpenMS can still consume).\n";
        }
        for (std::vector<std::string>::const_iterator pevref = pevid_ids.begin(); pevref != pevid_ids.end(); ++pevref)
        {
          sii_tmp += "\t\t\t\t\t<PeptideEvidenceRef peptideEvidence_ref=\"" +  std::string(*pevref) + "\"/>\n";
        }

        if (! hit.getPeakAnnotations().empty())
        {
          writeFragmentAnnotations_(sii_tmp, hit.getPeakAnnotations(), 5, false);
        }

        MetaInfoInterface copy_hit = hit;
        std::string st(it->getScoreType()); //scoretype

        // Only consume the score type here if it maps to a PSM-level search-engine statistic (MS:1001143 subtree).
        // Otherwise fall through to the dedicated alias branches below (e.g. Mascot/OMSSA), which would
        // otherwise be unreachable for score types that happen to be valid (non-score) PSI-MS term names.
        if (const auto* term = cv_.checkAndGetTermByName(st);
            term != nullptr && peptide_result_details_.contains(term->id))
        {
          (sii_tmp += "\t\t\t\t\t") += term->toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term->id);
        }
        else if (cv_.exists(st) && peptide_result_details_.contains(st))
        {
          const auto& term = cv_.getTerm(st);
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else if (st == "q-value" || st == "FDR")
        {
          const auto& term = cv_.getTermByName("PSM-level q-value");
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else if (st == "Posterior Error Probability")
        {
          const auto& term = cv_.getTermByName("percolator:PEP"); // 'percolaror' was not a typo in the code but in the cv.
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else if (st == "OMSSA")
        {
          const auto& term = cv_.getTermByName("OMSSA:evalue");
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else if (st == "Mascot")
        {
          const auto& term = cv_.getTermByName("Mascot:score");
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else if (st == "XTandem")
        {
          const auto& term = cv_.getTermByName("X\\!Tandem:hyperscore");
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else if (st == "SEQUEST")
        {
          const auto& term = cv_.getTermByName("Sequest:xcorr");
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else if (st == "MS-GF+")
        {
          const auto& term = cv_.getTermByName("MS-GF:RawScore");
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else if (st == Constants::UserParam::OPENPEPXL_SCORE)
        {
          const auto& term = cv_.getTermByName(st);
          (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(term.id);
        }
        else
        {
          std::string score_name_placeholder = st.empty()?"PSM-level search engine specific statistic":st;
          sii_tmp +=std::string(5, '\t') + cv_.getTermByName("PSM-level search engine specific statistic").toXMLString(cv_ns);
          sii_tmp += "\n" + std::string(5, '\t') + "<userParam name=\"" + score_name_placeholder
                       + "\" unitName=\"" + "xsd:double" + "\" value=\"" + sc + "\"/>";
          OPENMS_LOG_WARN << "Converting unknown score type to PSM-level search engine specific statistic from PSI controlled vocabulary.\n";
        }
        sii_tmp += "\n";

        copy_hit.removeMetaValue("calcMZ");
        copy_hit.removeMetaValue(Constants::UserParam::TARGET_DECOY);
        writeMetaInfos_(sii_tmp, copy_hit, 5);

        //~ sidres += "<cvParam accession=\"MS:1000796\" cvRef=\"PSI-MS\" value=\"55.835.842.3.dta\" name=\"spectrum title\"/>";
        sii_tmp += "\t\t\t\t</SpectrumIdentificationItem>\n";
        sidres += sii_tmp;
    }

    void MzIdentMLHandler::writeXLMSPeptideHit(const PeptideHit& hit,
                                                PeptideIdentificationList::const_iterator& it,
                                                const std::string& ppxl_linkid, std::map<std::string, std::string>& pep_ids,
                                                const std::string& cv_ns, std::set<std::string>& sen_set,
                                                std::map<std::string, std::string>& sen_ids,
                                                std::map<std::string, std::vector<std::string> >& pep_evis,
                                                std::map<std::string, double>& pp_identifier_2_thresh,
                                                double ppxl_crosslink_mass,
                                                std::map<std::string, std::string>& ppxl_specref_2_element,
                                                std::string& sid, bool alpha_peptide)
    {

      std::string pepid =  "PEP_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());

      AASequence peptide_sequence;
      if (alpha_peptide)
      {
        peptide_sequence = hit.getSequence();
      }
      else
      {
        peptide_sequence = AASequence::fromString(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE));
      }
      std::string pepi = peptide_sequence.toString();

      // The same peptide sequence (including mods and link position) can be used several times in different pairs
      // make pepi unique enough, so that PeptideEvidences are written for each case

      if (alpha_peptide)
      {
        pepi += "_MS:1002509";
      }
      else
      {
        pepi += "_MS:1002510";
      }
      if (hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) != "mono-link")  //sequence may contain more than one linker anchors; also code position linked
      {
        if (alpha_peptide)
        {
          pepi += "_" + hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1).toString();
        }
        else
        {
          pepi += "_" + hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2).toString();
        }
      }
      pepi += ppxl_linkid;

      bool duplicate = false;
      std::map<std::string, std::string>::iterator pit;
      // avoid duplicates in case with only one peptide
      if (hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) != "cross-link")
      {
        pit = pep_ids.find(pepi);
        if (pit != pep_ids.end())
        {
          duplicate = true;
        }
      }

      // TODO another criterion for ppxl: the same "donor" pep_id can only be reused in combination with the same "acceptor" pep_id,
      // for now I will just make an exemption for ppxl data, otherwise we get an invalid file as we have missing peptide pairings.
      // the redundancy should not increase too much, since pairings between exactly the same peptides for different spectra
      // should be a minority

      // avoid duplicate pairs in ppxl data
      // TODO access to pep_ids for both peptides from each pair necessary for PSMs
      // below code does not work at all yet

      // if (hit.metaValueExists("xl_chain") && it->getHits().size() == 2)
      //          {
      //            std::vector<PeptideHit> peps = it->getHits();
      //            std::string pepi2 = peps[1].getSequence().toString();
      //            std::map<std::string, std::string>::iterator pit = pep_pairs_ppxl.find(pepi);

      //            // last entry in vector
      //            std::pair<std::string, std::string> pepid_pairs_ppxl[pepid_pairs_ppxl.size()-1];

      //            if (pit != pep_pairs_ppxl.end())
      //            {
      //              if (pit->second == pep2)
      //              {
      //                duplicate = true;
      //              }
      //            }
      //          }

      if (!duplicate)
      {
        std::string p;
        //~ TODO simplify mod cv param write
        // write peptide id with conversion to universal, "human-readable" bracket string notation
        p +=std::string("\t<Peptide id=\"") + pepid + "\" name=\"" +
              peptide_sequence.toBracketString(false) + "\">\n\t\t<PeptideSequence>" + peptide_sequence.toUnmodifiedString() + std::string("</PeptideSequence>\n");

        const ResidueModification* n_term_mod = peptide_sequence.getNTerminalModification();
        if (n_term_mod != nullptr)
        {
          p += "\t\t<Modification location=\"0\">\n";
          std::string acc = n_term_mod->getUniModAccession();
          bool unimod = true;
          if (!acc.empty())
          {
            p += "\t\t\t<cvParam accession=\"UNIMOD:" + StringUtils::suffix(acc, ':');
          }
          else
          {
            acc = n_term_mod->getPSIMODAccession();
            p += "\t\t\t<cvParam accession=\"XLMOD:" + StringUtils::suffix(acc, ':');
            unimod = false;
          }
          p += "\" name=\"" + n_term_mod->getId();
          if (unimod)
          {
            p += R"(" cvRef="UNIMOD"/>)";
          }
          else
          {
            p += R"(" cvRef="XLMOD"/>)";
          }
          p += "\n\t\t</Modification>\n";
        }
        const ResidueModification* c_term_mod = peptide_sequence.getCTerminalModification();
        if (c_term_mod != nullptr)
        {
          p += "\t\t<Modification location=\"" + StringUtils::toStr(peptide_sequence.size()) + "\">\n";
          std::string acc = c_term_mod->getUniModAccession();
          bool unimod = true;
          if (!acc.empty())
          {
            p += "\t\t\t<cvParam accession=\"UNIMOD:" + StringUtils::suffix(acc, ':');
          }
          else
          {
            acc = c_term_mod->getPSIMODAccession();
            p += "\t\t\t<cvParam accession=\"XLMOD:" + StringUtils::suffix(acc, ':');
            unimod = false;
          }
          p += "\" name=\"" + c_term_mod->getId();
          if (unimod)
          {
            p += R"(" cvRef="UNIMOD"/>)";
          }
          else
          {
            p += R"(" cvRef="XLMOD"/>)";
          }
          p += "\n\t\t</Modification>\n";
        }
        for (Size i = 0; i < peptide_sequence.size(); ++i)
        {
          const ResidueModification* mod = peptide_sequence[i].getModification(); // "UNIMOD:" prefix??
          if (mod != nullptr)
          {
            p += "\t\t<Modification location=\"" + StringUtils::toStr(i + 1);
            p += "\" residues=\"" + peptide_sequence[i].getOneLetterCode();
            std::string acc = mod->getUniModAccession();
            if (!acc.empty())
            {
              p += "\">\n\t\t\t<cvParam accession=\"UNIMOD:" + StringUtils::suffix(acc, ':'); //TODO @all: do not exclusively use unimod ...
              p += "\" name=\"" + mod->getId();
              p += R"(" cvRef="UNIMOD"/>)";
              p += "\n\t\t</Modification>\n";
            }
            else
            {
              acc = mod->getPSIMODAccession();
              if (!acc.empty())
              {
                p += "\">\n\t\t\t<cvParam accession=\"XLMOD:" + StringUtils::suffix(acc, ':');
                p += "\" name=\"" +  mod->getId();
                p += R"(" cvRef="XLMOD"/>)";
                p += "\n\t\t</Modification>\n";
              }
              else
              {
                // We have an unknown modification, so lets write unknown
                // and at least try to write down the delta mass.
                if (mod->getDiffMonoMass() != 0.0)
                {
                  double diffmass = mod->getDiffMonoMass();
                  p += "\" monoisotopicMassDelta=\"" + StringUtils::toStr(diffmass);
                }
                else if (mod->getMonoMass() > 0.0)
                {
                  double diffmass = mod->getMonoMass() - peptide_sequence[i].getMonoWeight();
                  p += "\" monoisotopicMassDelta=\"" + StringUtils::toStr(diffmass);
                }
                p += "\">\n\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001460\" name=\"unknown modification\"/>";
                p += "\n\t\t</Modification>\n";
              }
            }
          }
          // Failsafe, if someone uses a new cross-linker (given these conditions, there MUST be a linker at this position, but it does not have a Unimod or XLMOD entry)
          else if (alpha_peptide && hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_MOD) && (static_cast<Size>(StringUtils::toInt32(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1).toString())) == i) && (hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE).toString() == "mono-link") )
          {
            p += "\t\t<Modification location=\"" + StringUtils::toStr(i + 1);
            p += "\" residues=\"" + std::string(peptide_sequence[i].getOneLetterCode());
            p += "\">\n\t\t\t<cvParam accession=\"XLMOD:XXXXX";
            p += "\" name=\"" +  hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD).toString();
            p += R"(" cvRef="XLMOD"/>)";
            p += "\n\t\t</Modification>\n";
          }
        }

        if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TYPE) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) != "mono-link")
        {
          int i = StringUtils::toInt32(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1).toString());
          if (alpha_peptide)
          {
            const CrossLinksDB* xl_db = CrossLinksDB::getInstance();
            std::vector<std::string> mods;
            if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA) == "N_TERM")
            {
              xl_db->searchModificationsByDiffMonoMass(mods, double(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS)), 0.0001, "", ResidueModification::N_TERM);
              if (!mods.empty())
              {
                p += "\t\t<Modification location=\"0";
              }
            }
            else if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA) == "C_TERM")
            {
              xl_db->searchModificationsByDiffMonoMass(mods, double(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS)), 0.0001, "", ResidueModification::C_TERM);
              if (!mods.empty())
              {
                p += "\t\t<Modification location=\"" + StringUtils::toStr(i + 2);
              }
            }
            else
            {
              xl_db->searchModificationsByDiffMonoMass(mods, double(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS)), 0.0001,std::string(hit.getSequence()[i].getOneLetterCode()), ResidueModification::ANYWHERE);
              if (!mods.empty())
              {
                p += "\t\t<Modification location=\"" + StringUtils::toStr(i + 1);
              }
            }

            std::string acc;
            std::string name;
            for (Size s = 0; s < mods.size(); ++s)
            {
              if (StringUtils::hasSubstring(mods[s], hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD).toString()))
              {
                const ResidueModification* mod = nullptr;
                try
                {
                  mod = xl_db->getModification(mods[s], peptide_sequence[i].getOneLetterCode(), ResidueModification::ANYWHERE);
                }
                catch (...)
                {
                  if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA) == "N_TERM")
                  {
                    mod = xl_db->getModification(mods[s], "", ResidueModification::N_TERM);
                  }
                  else if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA) == "C_TERM")
                  {
                    mod = xl_db->getModification(mods[s], "", ResidueModification::C_TERM);
                  }
                }
                // mod should never be null, but gcc complains (-Werror=maybe-uninitialized)
                if (mod != nullptr)
                {
                  acc = mod->getPSIMODAccession();
                }
                if (mod != nullptr)
                {
                  name = mod->getId();
                }
              }
              if (!acc.empty())
              {
                break;
              }
            }
            if ( acc.empty() && (!mods.empty()) ) // If ambiguity can not be resolved by xl_mod, just take one with the same mass diff from the database
            {
              const ResidueModification* mod = xl_db->getModification(std::string(peptide_sequence[i].getOneLetterCode()), mods[0], ResidueModification::ANYWHERE);
              acc = mod->getPSIMODAccession();
              name = mod->getId();
            }
            if (!acc.empty())
            {
              p += "\" residues=\"" + std::string(peptide_sequence[i].getOneLetterCode());
              p += "\" monoisotopicMassDelta=\"" + hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS).toString() + "\">\n";
              p += "\t\t\t<cvParam accession=\"" + acc + R"(" cvRef="XLMOD" name=")" + name + "\"/>\n";
            }
            else // if there is no matching modification in the database, write out a placeholder
            {
              p += "\t\t<Modification location=\"" + StringUtils::toStr(i + 1);
              p += "\" residues=\"" + std::string(peptide_sequence[i].getOneLetterCode());
              p += "\" monoisotopicMassDelta=\"" + hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS).toString() + "\">\n";
              p += "\t\t\t<cvParam accession=\"XLMOD:XXXXX\" cvRef=\"XLMOD\" name=\"" + hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD).toString() + "\"/>\n";
            }
          }
          else // xl_chain = "MS:1002510", acceptor, beta peptide
          {
            i = StringUtils::toInt32(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2).toString());
            if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA) == "N_TERM")
            {
              p += "\t\t<Modification location=\"0";
            }
            else if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA) == "C_TERM")
            {
              p += "\t\t<Modification location=\"" + StringUtils::toStr(peptide_sequence.size() + 2);
            }
            else
            {
              p += "\t\t<Modification location=\"" + StringUtils::toStr(i + 1);
            }
            p += "\" residues=\"" + std::string(peptide_sequence[i].getOneLetterCode());
            p += "\" monoisotopicMassDelta=\"0\">\n";
          }
          if (alpha_peptide)
          {
            p += "\t\t\t" + cv_.getTerm(std::string("MS:1002509")).toXMLString(cv_ns, DataValue(ppxl_linkid));
          }
          else
          {
            p += "\t\t\t" + cv_.getTerm(std::string("MS:1002510")).toXMLString(cv_ns, DataValue(ppxl_linkid));
          }
          p += "\n\t\t</Modification>\n";
        }
        if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TYPE) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "loop-link")
        {
          int i = StringUtils::toInt32(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2).toString());
          p += "\t\t<Modification location=\"" + StringUtils::toStr(i + 1);
          p += "\" residues=\"" + std::string(peptide_sequence[i].getOneLetterCode());
          p += "\" monoisotopicMassDelta=\"0";
          // ppxl crosslink loop xl_pos2 is always the acceptor ("MS:1002510")
          p += "\">\n\t\t\t" + cv_.getTerm("MS:1002510").toXMLString(cv_ns, DataValue(ppxl_linkid));
          p += "\n\t\t</Modification>\n";
        }
        p += "\t</Peptide>\n ";
        sen_set.insert(p);
        pep_ids.insert(std::make_pair(pepi, pepid));
      }
      else
      {
        pepid = pit->second;
      }

      std::vector<std::string> pevid_ids;
      if (!duplicate)
      {
        if (alpha_peptide)
        {
          std::vector<PeptideEvidence> peptide_evidences = hit.getPeptideEvidences();
          // TODO idXML allows peptide hits without protein references! Fails in that case - run PeptideIndexer first

          // TODO BETA PEPTIDE Protein Info
          for (std::vector<PeptideEvidence>::const_iterator pe = peptide_evidences.begin(); pe != peptide_evidences.end(); ++pe)
          {
            std::string pevid =  "PEV_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
            std::string dBSequence_ref;
            map<std::string, std::string>::const_iterator pos = sen_ids.find(pe->getProteinAccession());
            if (pos != sen_ids.end())
            {
              dBSequence_ref = pos->second;
            }
            else
            {
              OPENMS_LOG_ERROR << "Error: Missing or invalid protein reference for peptide '" << pepi << "': '" << pe->getProteinAccession() << "' - skipping." << endl;
              continue;
            }
            std::string idec;
            if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA))
            {
              idec =std::string(boost::lexical_cast<std::string>(StringUtils::hasSubstring(hit.getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA).toString(), "decoy")));
            }

            std::string e;
            std::string nc_termini = "-";    // character for N- and C-termini as specified in mzIdentML
            e += "\t<PeptideEvidence id=\"" + pevid + "\" peptide_ref=\"" + pepid + "\" dBSequence_ref=\"" + dBSequence_ref + "\"";

            if (pe->getAAAfter() != PeptideEvidence::UNKNOWN_AA)
            {
              e += " post=\"" + (pe->getAAAfter() == PeptideEvidence::C_TERMINAL_AA ? nc_termini : StringUtils::toStr(pe->getAAAfter())) + "\"";
            }
            if (pe->getAABefore() != PeptideEvidence::UNKNOWN_AA)
            {
              e += " pre=\"" + (pe->getAABefore() == PeptideEvidence::N_TERMINAL_AA ? nc_termini : StringUtils::toStr(pe->getAABefore())) + "\"";
            }
            if (pe->getStart() != PeptideEvidence::UNKNOWN_POSITION)
            {
              e += " start=\"" + StringUtils::toStr(pe->getStart() + 1) + "\"";
            }
            else if (hit.metaValueExists("start"))
            {
              e += " start=\"" + StringUtils::toStr( int(hit.getMetaValue("start")) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no start position of peptide hit in protein sequence.\n";
            }
            if (pe->getEnd() != PeptideEvidence::UNKNOWN_POSITION)
            {
              e += " end=\"" + StringUtils::toStr(pe->getEnd() + 1) + "\"";
            }
            else if (hit.metaValueExists("end"))
            {
              e += " end=\"" + StringUtils::toStr( int(hit.getMetaValue("end")) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no end position of peptide hit in protein sequence.\n";
            }
            if (!idec.empty())
            {
              e += " isDecoy=\"" + std::string(idec)+ "\"";
            }
            e += "/>\n";
            sen_set.insert(e);
            pevid_ids.push_back(pevid);
          }
          pep_evis.insert(make_pair(pepi, pevid_ids));
        }
        else // acceptor, beta peptide, does not have its own PeptideHit and PeptideEvidences
        {
          StringList prot = ListUtils::create<std::string>(StringUtils::toStr(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS)));
          StringList pre = ListUtils::create<std::string>(StringUtils::toStr(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_PRE)));
          StringList post = ListUtils::create<std::string>(StringUtils::toStr(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_POST)));
          StringList start = ListUtils::create<std::string>(StringUtils::toStr(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_START)));
          StringList end = ListUtils::create<std::string>(StringUtils::toStr(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_END)));
          for (Size ev = 0; ev < pre.size(); ++ev)
          {
            std::string pevid =  "PEV_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
            std::string dBSequence_ref;
            map<std::string, std::string>::const_iterator pos = sen_ids.find(prot[ev]);
            if (pos != sen_ids.end())
            {
              dBSequence_ref = pos->second;
            }
            else
            {
              OPENMS_LOG_ERROR << "Error: Missing or invalid protein reference for peptide '" << pepi << "': '" << prot[ev] << "' - skipping." << endl;
              continue;
            }
            std::string idec;
            if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA))
            {
              idec =std::string(boost::lexical_cast<std::string>(StringUtils::hasSubstring(hit.getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA).toString(), "decoy")));
            }

            std::string e;
            std::string nc_termini = "-";    // character for N- and C-termini as specified in mzIdentML
            e += "\t<PeptideEvidence id=\"" + pevid + "\" peptide_ref=\"" + pepid + "\" dBSequence_ref=\"" + dBSequence_ref + "\"";

            if (post[ev] !=StringUtils::toStr(PeptideEvidence::UNKNOWN_AA))
            {
              e += " post=\"" + (post[ev] ==StringUtils::toStr(PeptideEvidence::C_TERMINAL_AA) ? nc_termini : post[ev]) + "\"";
            }
            if (pre[ev] !=StringUtils::toStr(PeptideEvidence::UNKNOWN_AA))
            {
              e += " pre=\"" + (pre[ev] ==StringUtils::toStr(PeptideEvidence::N_TERMINAL_AA) ? nc_termini : pre[ev]) + "\"";
            }
            if (start[ev] !=StringUtils::toStr(PeptideEvidence::UNKNOWN_POSITION))
            {
              e += " start=\"" + StringUtils::toStr(StringUtils::toInt32(start[ev]) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no start position of peptide hit in protein sequence.\n";
            }
            if (end[ev] !=StringUtils::toStr(PeptideEvidence::UNKNOWN_POSITION))
            {
              e += " end=\"" + StringUtils::toStr(StringUtils::toInt32(end[ev]) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no end position of peptide hit in protein sequence.\n";
            }
            if (!idec.empty())
            {
              e += " isDecoy=\"" + std::string(idec)+ "\"";
            }
            e += "/>\n";
            sen_set.insert(e);
            pevid_ids.push_back(pevid);
          }
          pep_evis.insert(make_pair(pepi, pevid_ids));
        }
      }
      else
      {
        pevid_ids =  pep_evis[pepi];
      }

      std::string r = StringUtils::toStr(hit.getRank() + 1); //rank
      std::string sc = StringUtils::toStr(hit.getScore());
      if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_RANK))
      {
        r = hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK).toString();  // ppxl remove xl_rank later (in copy_hit)
      }

      //Calculated mass to charge for cross-linked is both peptides + linker
      double calc_ppxl_mass = hit.getSequence().getMonoWeight();
      if (hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "cross-link")
      {
        calc_ppxl_mass += ppxl_crosslink_mass + AASequence::fromString(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE)).getMonoWeight();
      }
      // if xl_mod and xl_mass MetaValues exist, then the mass of the mono-link could not be set as a AASequence modification and will not be considered by .getMonoWeight
      else if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_MOD) && hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_MASS))
      {
        calc_ppxl_mass += double(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS));
      }
      std::string cmz =StringUtils::toStr( ((calc_ppxl_mass) +  (hit.getCharge() * Constants::PROTON_MASS_U)) / hit.getCharge()); //calculatedMassToCharge

      if (sc.empty())
      {
        sc = "NA";
        OPENMS_LOG_WARN << "No score assigned to this PSM: " /*<< hit.getSequence().toString()*/ << '\n';
      }
      std::string c = StringUtils::toStr(hit.getCharge()); //charge

      std::string pte;
      if (hit.metaValueExists("pass_threshold"))
      {
        pte = boost::lexical_cast<std::string>(hit.getMetaValue("pass_threshold"));
      }
      else
      {
        auto thresh_it = pp_identifier_2_thresh.find(it->getIdentifier());
        if (thresh_it != pp_identifier_2_thresh.end() && thresh_it->second != 0.0)
        {
          double th = thresh_it->second;
          //threshold was 'set' in proteinIdentification (!= default value of member, now check pass
          pte = boost::lexical_cast<std::string>(it->isHigherScoreBetter() ? hit.getScore() > th : hit.getScore() < th); //passThreshold-eval
        }
        else
        {
          pte = "1";
        }
      }

      //write SpectrumIdentificationItem elements
      std::string emz = StringUtils::toStr(it->getMZ());
      std::string sii = "SII_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
      std::string sii_tmp;
      sii_tmp +=std::string("\t\t\t\t<SpectrumIdentificationItem passThreshold=\"")
              + pte + "\" rank=\"" + r + "\" peptide_ref=\""
              + pepid + "\" calculatedMassToCharge=\"" + cmz
              + "\" experimentalMassToCharge=\"" + emz
              + "\" chargeState=\"" + c +  "\" id=\""
              + sii + "\">\n";

      if (pevid_ids.empty())
      {
        OPENMS_LOG_WARN << "PSM without peptide evidence registered in the given search database found. This will cause an invalid mzIdentML file (which OpenMS can still consume).\n";
      }
      for (std::vector<std::string>::const_iterator pevref = pevid_ids.begin(); pevref != pevid_ids.end(); ++pevref)
      {
        sii_tmp += "\t\t\t\t\t<PeptideEvidenceRef peptideEvidence_ref=\"" +  std::string(*pevref) + "\"/>\n";
      }

      if (! hit.getPeakAnnotations().empty() && alpha_peptide)
      {
        // is_ppxl = true
        writeFragmentAnnotations_(sii_tmp, hit.getPeakAnnotations(), 5, true);
      }

      MetaInfoInterface copy_hit = hit;
      std::string st(it->getScoreType()); //scoretype

      // Only consume the score type here if it maps to a PSM-level search-engine statistic (MS:1001143 subtree);
      // otherwise fall through to the dedicated alias branches below.
      if (const auto* term = cv_.checkAndGetTermByName(st);
          term != nullptr && peptide_result_details_.contains(term->id))
      {
        (sii_tmp += "\t\t\t\t\t") += term->toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(term->id);
      }
      else if (cv_.exists(st) && peptide_result_details_.contains(st))
      {
        const auto& term = cv_.getTerm(st);
        (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(term.id);
      }
      else if (st == "q-value" || st == "FDR")
      {
        const auto& term = cv_.getTermByName("PSM-level q-value");
        (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(term.id);
      }
      else if (st == "Posterior Error Probability")
      {
        const auto& term = cv_.getTermByName("percolator:PEP"); // 'percolaror' was not a typo in the code but in the cv.
        (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(term.id);
      }
      else if (st == Constants::UserParam::OPENPEPXL_SCORE)
      {
        const auto& term = cv_.getTermByName(st);
        (sii_tmp += "\t\t\t\t\t") += term.toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(term.id);
      }
      else
      {
        std::string score_name_placeholder = st.empty()?"PSM-level search engine specific statistic":st;
        sii_tmp +=std::string(5, '\t') + cv_.getTermByName("PSM-level search engine specific statistic").toXMLString(cv_ns);
        sii_tmp += "\n" + std::string(5, '\t') + "<userParam name=\"" + score_name_placeholder
                     + "\" unitName=\"" + "xsd:double" + "\" value=\"" + sc + "\"/>";
        OPENMS_LOG_WARN << "Converting unknown score type to PSM-level search engine specific statistic from PSI controlled vocabulary.\n";
      }
      sii_tmp += "\n";

      copy_hit.removeMetaValue("calcMZ");
      copy_hit.removeMetaValue(Constants::UserParam::TARGET_DECOY);

      // TODO this would be the correct way, but need to adjust parsing as well
      if (copy_hit.metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ) || copy_hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "cross-link")
      {
        sii_tmp +=  "\t\t\t\t\t" + cv_.getTerm("MS:1002511").toXMLString(cv_ns, ppxl_linkid) + "\n"; // cross-linked spectrum identification item
      }

      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD);
      copy_hit.removeMetaValue("xl_chain");
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS);
      copy_hit.removeMetaValue("protein_references");
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF);
      copy_hit.removeMetaValue(Constants::UserParam::SPECTRUM_REFERENCE);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_PRE);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_POST);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_START);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_END);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE);
      copy_hit.removeMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS);

      writeMetaInfos_(sii_tmp, copy_hit, 5);

      //~ sidres += "<cvParam accession=\"MS:1000796\" cvRef=\"PSI-MS\" value=\"55.835.842.3.dta\" name=\"spectrum title\"/>";
      sii_tmp += "\t\t\t\t</SpectrumIdentificationItem>\n";

      const double rt = it->getRT();
      std::string ert = rt == rt ? StringUtils::toStr(rt) : "nan";
      DataValue rtcv(ert);
      rtcv.setUnit(10); // id: UO:0000010 name: second
      rtcv.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
      StringUtils::substitute(sii_tmp, "</SpectrumIdentificationItem>",
                              "\t" + cv_.getTermByName("retention time").toXMLString(cv_ns, rtcv) + "\n\t\t\t\t</SpectrumIdentificationItem>\n");
      ppxl_specref_2_element[sid] += sii_tmp;
      if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT) && hit.metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ))
      {
        // ppxl : TODO remove if existing <Fragmentation/>block
        { std::string from = std::string("experimentalMassToCharge=\"") + std::string(emz);
          std::string to = "experimentalMassToCharge=\"" + StringUtils::toStr(hit.getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ));
          StringUtils::substitute(sii_tmp, from, to); } // mz
        sii_tmp = StringUtils::substitute(sii_tmp, sii,std::string("SII_") + StringUtils::toStr(UniqueIdGenerator::getUniqueId())); // uid
        sii_tmp = StringUtils::substitute(sii_tmp, "value=\"" + ert, "value=\"" + StringUtils::toStr(hit.getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT)));

        ProteinIdentification::SearchParameters search_params = cpro_id_->front().getSearchParameters();
        // default avoids a ConversionError when the isoshift meta value is absent/empty
        double iso_shift = StringUtils::toDouble(StringUtils::toStr(search_params.getMetaValue("cross_link:mass_isoshift", "0")));
        double cmz_heavy = StringUtils::toDouble(cmz) + (iso_shift / hit.getCharge());

        { std::string from2 = std::string("calculatedMassToCharge=\"") + std::string(cmz);
          std::string to2 = "calculatedMassToCharge=\"" + StringUtils::toStr(cmz_heavy);
          StringUtils::substitute(sii_tmp, from2, to2); }

        ppxl_specref_2_element[sid] += sii_tmp;
      }
    }
  } // namespace OpenMS
