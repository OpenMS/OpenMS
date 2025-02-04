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

    MzIdentMLHandler::MzIdentMLHandler(const std::vector<ProteinIdentification>& pro_id, const std::vector<PeptideIdentification>& pep_id, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      //~ ms_exp_(0),
      pro_id_(nullptr),
      pep_id_(nullptr),
      cpro_id_(&pro_id),
      cpep_id_(&pep_id)
    {
      cv_.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      unimod_.loadFromOBO("PSI-MS", File::find("/CV/unimod.obo"));
    }

    MzIdentMLHandler::MzIdentMLHandler(std::vector<ProteinIdentification>& pro_id, std::vector<PeptideIdentification>& pep_id, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      //~ ms_exp_(0),
      pro_id_(&pro_id),
      pep_id_(&pep_id),
      cpro_id_(nullptr),
      cpep_id_(nullptr)
    {
      cv_.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      unimod_.loadFromOBO("PSI-MS", File::find("/CV/unimod.obo"));
    }

    //~ TODO create MzIdentML instances from MSExperiment which contains much of the information yet needed
    //~ MzIdentMLHandler(const PeakMap& mx, const String& filename, const String& version, const ProgressLogger& logger)
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

      static set<String> to_ignore;
      if (to_ignore.empty())
      {
        to_ignore.insert("peptideSequence");
      }

      if (to_ignore.find(tag_) != to_ignore.end())
      {
        return;
      }

      //determine parent tag
      String parent_tag;
      if (open_tags_.size() > 1)
      {
        parent_tag = *(open_tags_.end() - 2);
      }
      String parent_parent_tag;
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

        String value, unit_accession, cv_ref;
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
        String name;
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
        current_id_hit_.setRank(attributeAsInt_(attributes, "rank"));

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

        String string_value("");
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
        String customizations = sm_.convert(chars);
        // TODO write customizations to Software
        return;
      }

      if (tag_ == "seq")
      {
        String seq = sm_.convert(chars);
        actual_protein_.setSequence(seq);
        return;
      }

      if (tag_ == "peptideSequence")
      {
        String pep = sm_.convert(chars);
        actual_peptide_ = AASequence::fromString(pep);
        return;
      }

      //error(LOAD, "MzIdentMLHandler::characters: Unknown character section found: '" + tag_ + "', ignoring.");
    }

    void MzIdentMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      static set<String> to_ignore;
      if (to_ignore.empty())
      {
        to_ignore.insert("mzIdentML");
        to_ignore.insert("cvParam");
      }

      tag_ = sm_.convert(qname);
      open_tags_.pop_back();

      if (to_ignore.find(tag_) != to_ignore.end())
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

    void MzIdentMLHandler::handleCVParam_(const String& /* parent_parent_tag*/, const String& parent_tag, const String& accession, /* const String& name, */ /* const String& value, */ const xercesc::Attributes& attributes, const String& cv_ref /* , const String& unit_accession */)
    {
      if (parent_tag == "Modification")
      {
        if (cv_ref == "UNIMOD")
        {
          set<const ResidueModification*> mods;
          Int loc = numeric_limits<Int>::max();
          if (optionalAttributeAsInt_(loc, attributes, "location"))
          {
            String uni_mod_id = accession.suffix(':');
            String residues;
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
            String message = String("Modification '") + accession + "' is unknown.";
            throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
          }
        }
      }
    }

    void MzIdentMLHandler::writeTo(std::ostream& os)
    {
      String cv_ns = cv_.name();
      String inputs_element;
      std::map<String,String> /* peps, pepevis, */ sil_map, sil_2_date;
      std::set<String> sen_set, sof_set, sip_set;
      std::map<String, String> sdb_ids, sen_ids, sof_ids, sdat_ids, pep_ids;
      //std::map<String, String> pep_pairs_ppxl;
      std::map<String, double> pp_identifier_2_thresh;
      //std::vector< std::pair<String, String> > pepid_pairs_ppxl;

      // file type-specific definitions needed for SpectraData element:
      std::map<FileTypes::Type, std::pair<String, String> > formats_map;
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
      inputs_element += String("\t<Inputs>\n");
      String spectra_data, search_database;

      /*
      1st: iterate over proteinidentification vector
      */
      //TODO read type of crosslink reagent from settings
      bool is_ppxl = false;
      for (std::vector<ProteinIdentification>::const_iterator it = cpro_id_->begin(); it != cpro_id_->end(); ++it)
      {
        //~ collect analysissoftware in this loop - does not go into inputelement
        String sof_id;
        String sof_name = String(it->getSearchEngine());
        std::map<String, String>::iterator soit = sof_ids.find(sof_name);
        String osecv;
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
          sof_id = "SOF_" + String(UniqueIdGenerator::getUniqueId());
          //~ TODO consider not only searchengine but also version!
          String sost = String("\t<AnalysisSoftware version=\"") + String(it->getSearchEngineVersion()) + String("\" name=\"") + sof_name +  String("\" id=\"") + sof_id + String("\">\n") + String("\t\t<SoftwareName>\n");
          sost += "\t\t\t" + cv_.getTermByName(osecv).toXMLString(cv_ns);
          sost += String("\n\t\t</SoftwareName>\n\t</AnalysisSoftware>\n");
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

        String thcv;
        pp_identifier_2_thresh.insert(make_pair(it->getIdentifier(),it->getSignificanceThreshold()));
        if (it->getSignificanceThreshold() != 0.0)
        {
          thcv = cv_.getTermByName("PSM-level statistical threshold").toXMLString(cv_ns, String(it->getSignificanceThreshold()));
        }
        else
        {
          thcv = cv_.getTermByName("no threshold").toXMLString(cv_ns);
        }
        // TODO add other software than searchengine for evidence trace

        // get a map from identifier to match OpenMS Protein/PeptideIdentification match string;
        String sil_id =  "SIL_" + String(UniqueIdGenerator::getUniqueId());
        pp_identifier_2_sil_.insert(make_pair(it->getIdentifier(), sil_id));

        //~ collect SpectrumIdentificationProtocol for analysisprotocol in this loop - does not go into inputelement
        String sip_id = "SIP_" + String(UniqueIdGenerator::getUniqueId());
        sil_2_sip_.insert(make_pair(sil_id, sip_id));

        String sip = "\t<SpectrumIdentificationProtocol id=\"" + String(sip_id) + "\" analysisSoftware_ref=\"" + String(sof_id) + "\">\n";
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
        sip += String(3, '\t') + R"(<userParam name="charges" unitName="xsd:string" value=")" + search_params.charges + "\"/>\n";
//        sip += String(3, '\t') + "<userParam name=\"" + "missed_cleavages" + "\" unitName=\"" + "xsd:integer" + "\" value=\"" + String(it->getSearchParameters().missed_cleavages) + "\"/>" + "\n";
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
        sip += String("\t\t<FragmentTolerance>\n");
        String unit_str = R"(unitCvRef="UO" unitName="dalton" unitAccession="UO:0000221")";
        if (search_params.fragment_mass_tolerance_ppm)
        {
          unit_str = R"(unitCvRef="UO" unitName="parts per million" unitAccession="UO:0000169")";
        }
        sip += String(3, '\t') + R"(<cvParam accession="MS:1001412" name="search tolerance plus value" )" + unit_str + R"( cvRef="PSI-MS" value=")" + String(search_params.fragment_mass_tolerance) + "\"/>\n";
        sip += String(3, '\t') + R"(<cvParam accession="MS:1001413" name="search tolerance minus value" )" + unit_str + R"( cvRef="PSI-MS" value=")" + String(search_params.fragment_mass_tolerance) + "\"/>\n";
        sip += String("\t\t</FragmentTolerance>\n");
        sip += String("\t\t<ParentTolerance>\n");
        unit_str = R"(unitCvRef="UO" unitName="dalton" unitAccession="UO:0000221")";
        if (search_params.precursor_mass_tolerance_ppm)
        {
          unit_str = R"(unitCvRef="UO" unitName="parts per million" unitAccession="UO:0000169")";
        }
        sip += String(3, '\t') + R"(<cvParam accession="MS:1001412" name="search tolerance plus value" )" + unit_str + R"( cvRef="PSI-MS" value=")" + String(search_params.precursor_mass_tolerance) + "\"/>\n";
        sip += String(3, '\t') + R"(<cvParam accession="MS:1001413" name="search tolerance minus value" )" + unit_str + R"( cvRef="PSI-MS" value=")" + String(search_params.precursor_mass_tolerance) + "\"/>\n";
        sip += String("\t\t</ParentTolerance>\n");
        sip += String("\t\t<Threshold>\n\t\t\t") + thcv + "\n";
        sip += String("\t\t</Threshold>\n");
        sip += String("\t</SpectrumIdentificationProtocol>\n");
        sip_set.insert(sip);
        
        // empty date would lead to XML schema validation error:
        DateTime date_time = it->getDateTime();
        if (!date_time.isValid()) 
        { 
          date_time = DateTime::now(); 
        }
        sil_2_date.insert(make_pair(sil_id, String(date_time.getDate() + "T" + date_time.getTime())));

        //~ collect SpectraData element for each ProteinIdentification
        String sdat_id;
        StringList sdat_files;
        String sdat_file(it->getMetaValue("spectra_data"));

        if (sdat_file.empty())
        {
          sdat_file = String("UNKNOWN");
        }
        else
        {
          sdat_file = trimOpenMSfileURI(sdat_file);
        }
        std::map<String, String>::iterator sdit = sdat_ids.find(sdat_file); //this part is strongly connected to AnalysisCollection write part
        if (sdit == sdat_ids.end())
        {
          sdat_id = "SDAT_" + String(UniqueIdGenerator::getUniqueId());

          FileTypes::Type type = FileHandler::getTypeByFileName(sdat_file);
          if (formats_map.find(type) == formats_map.end()) type = FileTypes::MZML; // default

          //xml
          spectra_data += String("\t\t<SpectraData location=\"") + sdat_file + String("\" id=\"") + sdat_id + String("\">");
          spectra_data += String("\n\t\t\t<FileFormat>\n");
          spectra_data += String(4, '\t') + cv_.getTermByName(formats_map[type].first).toXMLString(cv_ns);
          spectra_data += String("\n\t\t\t</FileFormat>\n\t\t\t<SpectrumIDFormat>\n");
          spectra_data += String(4, '\t') + cv_.getTermByName(formats_map[type].second).toXMLString(cv_ns);
          spectra_data += String("\n\t\t\t</SpectrumIDFormat>\n\t\t</SpectraData>\n");

          sdat_ids.insert(make_pair(sdat_file, sdat_id));
          ph_2_sdat_.insert(make_pair(it->getIdentifier(), sdat_id));
        }
        else
        {
          sdat_id = sdit->second;
        }
        sil_2_sdat_.insert(make_pair(sil_id,  sdat_id));

        //~ collect SearchDatabase element for each ProteinIdentification
        String sdb_id;
        String sdb_file(search_params.db); //TODO @mths for several IdentificationRuns this must be something else, otherwise for two of the same db just one will be created
        std::map<String, String>::iterator dbit = sdb_ids.find(sdb_file);
        if (dbit == sdb_ids.end())
        {
          sdb_id = "SDB_"+ String(UniqueIdGenerator::getUniqueId());

          search_database += String("\t\t<SearchDatabase ");
          search_database += String("location=\"") + sdb_file + "\" ";
          if (!String(search_params.db_version).empty())
          {
            search_database += String("version=\"") + String(search_params.db_version) + "\" ";
          }
          search_database += String("id=\"") + sdb_id + String("\">\n\t\t\t<FileFormat>\n");
          //TODO Searchdb file format type cvParam handling
          search_database += String(4, '\t') + cv_.getTermByName("FASTA format").toXMLString(cv_ns);
          search_database += String("\n\t\t\t</FileFormat>\n\t\t\t<DatabaseName>\n\t\t\t\t<userParam name=\"") + sdb_file + String("\"/>\n\t\t\t</DatabaseName>\n");
          // "MS:1001029" was removed from the "search_params" copy!
          if (it->getSearchParameters().metaValueExists("MS:1001029"))
          {
            search_database += String(3, '\t') + cv_.getTerm("MS:1001029").toXMLString(cv_ns, it->getSearchParameters().getMetaValue("MS:1001029")) + String("\n");
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
          String enid;
          std::map<String, String>::iterator enit = sen_ids.find(String(jt->getAccession()));
          if (enit == sen_ids.end())
          {
            String entry;
            enid = "PROT_" + String(UniqueIdGenerator::getUniqueId()); //TODO IDs from metadata or where its stored at read in;
            String enst(jt->getAccession());

            entry += "\t<DBSequence accession=\"" + enst + "\" ";
            entry += "searchDatabase_ref=\"" + sdb_id + "\" ";
            String s = String(jt->getSequence());
            if (!s.empty())
            {
              entry += "length=\"" + String(jt->getSequence().length()) + "\" ";
            }
            entry += String("id=\"") + String(enid) + String("\">\n");
            if (!s.empty())
            {
              entry += "\t<Seq>" + s + "</Seq>\n";
            }
            entry += "\t\t" + cv_.getTermByName("protein description").toXMLString(cv_ns, enst);
            entry += "\n\t</DBSequence>\n";

            sen_ids.insert(std::pair<String, String>(enst, enid));
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
      std::map<String, String> ppxl_specref_2_element; //where the SII will get added for one spectrum reference
      std::map<String, std::vector<String> > pep_evis; //maps the sequence to the corresponding evidence elements for the next scope
      for (std::vector<PeptideIdentification>::const_iterator it = cpep_id_->begin(); it != cpep_id_->end(); ++it)
      {
        String emz(it->getMZ());
        const double rt = it->getRT();
        String ert = rt == rt ? String(rt) : "nan";

        String sid = it->getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE);
        if (sid.empty())
        {
          sid = String(it->getMetaValue("spectrum_id"));
          if (sid.empty())
          {
              if (it->getMZ() != it->getMZ())
            {
              emz = "nan";
              OPENMS_LOG_WARN << "Found no spectrum reference and no m/z position of identified spectrum! You are probably converting from an old format with insufficient data provision. Setting 'nan' - downstream applications might fail unless you set the references right." << std::endl;
            }
            if (it->getRT() != it->getRT())
            {
              ert = "nan";
              OPENMS_LOG_WARN << "Found no spectrum reference and no RT position of identified spectrum! You are probably converting from an old format with insufficient data provision. Setting 'nan' - downstream applications might fail unless you set the references right." << std::endl;
            }
            sid = String("MZ:") + emz + String("@RT:") + ert;
          }
        }
        if (is_ppxl && it->metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF))
        {
          sid.append("," + String(it->getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF)));
        }
        String sidres;
        String sir = "SIR_" + String(UniqueIdGenerator::getUniqueId());
        String sdr = sdat_ids.begin()->second;
        std::map<String, String>::iterator pfo = ph_2_sdat_.find(it->getIdentifier());
        if (pfo != ph_2_sdat_.end())
        {
          sdr = pfo->second;
        }
        else
        {
          OPENMS_LOG_WARN << "Falling back to referencing first spectrum file given because file or identifier could not be mapped." << std::endl;
        }

        sidres += String("\t\t\t<SpectrumIdentificationResult spectraData_ref=\"")
        //multi identification runs lookup from file_origin here
                + sdr + String("\" spectrumID=\"")
                + sid + String("\" id=\"") + sir + String("\">\n");

        if (is_ppxl)
        {
          if (ppxl_specref_2_element.find(sid)==ppxl_specref_2_element.end())
          {
            ppxl_specref_2_element[sid] = sidres;
          }
        }

        //map.begin access ok here because make sure at least one "UNKOWN" element is in the sdats_ids map

        double ppxl_crosslink_mass(0);
        if (is_ppxl)
        {
          ProteinIdentification::SearchParameters search_params = cpro_id_->front().getSearchParameters();
          ppxl_crosslink_mass = String(search_params.getMetaValue("cross_link:mass")).toDouble();
        }

        for (std::vector<PeptideHit>::const_iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
        {
          if (!is_ppxl)
          {
            MzIdentMLHandler::writePeptideHit(*jt, it, pep_ids, cv_ns, sen_set, sen_ids, pep_evis, pp_identifier_2_thresh, sidres);
          }
          else
          {
            String ppxl_linkid = UniqueIdGenerator::getUniqueId();
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
          std::map<String, String>::const_iterator ps_it = pp_identifier_2_sil_.find(it->getIdentifier());
          if (ps_it != pp_identifier_2_sil_.end())
          {
            std::map<String, String>::iterator sil_it = sil_map.find(ps_it->second);
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
            OPENMS_LOG_ERROR << "encountered a PeptideIdentification which is not linked to any ProteinIdentification" << std::endl;
          }
        }

      }
      // ppxl - write spectrumidentificationresult closing tags!
      if (is_ppxl)
      {
        for (std::map<String, String>::iterator it = ppxl_specref_2_element.begin(); it != ppxl_specref_2_element.end(); ++it)
        {
          it->second += "\t\t\t</SpectrumIdentificationResult>\n";
          std::map<String, String>::const_iterator ps_it = pp_identifier_2_sil_.begin();

          if (ps_it != pp_identifier_2_sil_.end())
          {
            std::map<String, String>::iterator sil_it = sil_map.find(ps_it->second);
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
            OPENMS_LOG_ERROR << "encountered a PeptideIdentification crosslink information which is not linked to any ProteinIdentification" << std::endl;
          }
        }
      }

      //--------------------------------------------------------------------------------------------
      // XML header
      //--------------------------------------------------------------------------------------------
      String v_s = "1.1.0";
      if (is_ppxl)
      {
         v_s = "1.2.0";
      }
      os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
         << "<MzIdentML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
         << "\txsi:schemaLocation=\"http://psidev.info/psi/pi/mzIdentML/"<< v_s.substr(0,v_s.size()-2) <<" "
         << "https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML"<< v_s <<".xsd\"\n"
         << "\txmlns=\"http://psidev.info/psi/pi/mzIdentML/"<< v_s.substr(0,v_s.size()-2) <<"\"\n"
         << "\tversion=\"" << v_s << "\"\n";
      os << "\tid=\"OpenMS_" << String(UniqueIdGenerator::getUniqueId()) << "\"\n"
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
      for (std::set<String>::const_iterator sof = sof_set.begin(); sof != sof_set.end(); ++sof)
      {
        os << *sof;
      }

      std::map<String, String>::iterator soit = sof_ids.find("TOPP software");
      if (soit == sof_ids.end())
      {
        os << "\t<AnalysisSoftware version=\"OpenMS TOPP v"<< VersionInfo::getVersion() <<R"(" name="TOPP software" id=")" << String("SOF_") << String(UniqueIdGenerator::getUniqueId()) << "\">\n"
           << "\t\t<SoftwareName>\n\t\t\t" << cv_.getTermByName("TOPP software").toXMLString(cv_ns) << "\n\t\t</SoftwareName>\n\t</AnalysisSoftware>\n";
      }
      os << "</AnalysisSoftwareList>\n";

      //--------------------------------------------------------------------------------------------
      // SequenceCollection
      //--------------------------------------------------------------------------------------------
      os << "<SequenceCollection>\n";
      for (std::set<String>::const_iterator sen = sen_set.begin(); sen != sen_set.end(); ++sen)
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
      for (std::map<String,String>::const_iterator pp2sil_it = pp_identifier_2_sil_.begin(); pp2sil_it != pp_identifier_2_sil_.end(); ++pp2sil_it)
      {
          String entry = String("\t<SpectrumIdentification id=\"SI_") + pp2sil_it->first + String("\" spectrumIdentificationProtocol_ref=\"")
                           + sil_2_sip_[pp2sil_it->second] + String("\" spectrumIdentificationList_ref=\"") + pp2sil_it->second
                           + String("\" activityDate=\"") + sil_2_date[pp2sil_it->second]
                           + String("\">\n")
                            //if crosslink +cvparam crosslink search performed
                           + "\t\t<InputSpectra spectraData_ref=\"" + sil_2_sdat_[pp2sil_it->second] + "\"/>\n" // spd_ids.insert(std::pair<String, UInt64>(sdst, sdid));
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
      for (std::set<String>::const_iterator sip = sip_set.begin(); sip != sip_set.end(); ++sip)
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
      for (std::map<String,String>::const_iterator sil_it = sil_map.begin(); sil_it != sil_map.end(); ++sil_it)
      {
        os << "\t\t<SpectrumIdentificationList id=\"" << sil_it->first << String("\">\n");
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

    void MzIdentMLHandler::writeMetaInfos_(String& s, const MetaInfoInterface& meta, UInt indent) const
    {
      //TODO @mths: write those metas with their name in the cvs loaded as CVs!
      if (meta.isMetaEmpty())
      {
        return;
      }
      std::vector<String> keys;
      meta.getKeys(keys);

      for (Size i = 0; i != keys.size(); ++i)
      {
        if (cv_.exists(keys[i]))
        {
          ControlledVocabulary::CVTerm a = cv_.getTerm(keys[i]);
          s += String(indent, '\t') + a.toXMLString("PSI-MS", (String)(meta.getMetaValue(keys[i]))) + "\n";
        }
        else
        {
          s += String(indent, '\t') + "<userParam name=\"" + keys[i] + "\" unitName=\"";

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
          s += "\" value=\"" + (String)(d) + "\"/>\n";
        }
      }
    }

    void MzIdentMLHandler::writeEnzyme_(String& s, const DigestionEnzymeProtein& enzy, UInt miss, UInt indent) const
    {
      String cv_ns = cv_.name();
      s += String(indent, '\t') + "<Enzymes independent=\"false\">\n";
      s += String(indent + 1, '\t') + "<Enzyme missedCleavages=\"" + String(miss) + "\" id=\"" + String("ENZ_") + String(UniqueIdGenerator::getUniqueId()) + "\">\n";
      s += String(indent + 2, '\t') + "<EnzymeName>\n";
      const String& enzymename = enzy.getName();
      if (cv_.hasTermWithName(enzymename))
      {
        s += String(indent + 3, '\t') + cv_.getTermByName(enzymename).toXMLString(cv_ns) + "\n";
      }
      else if (enzymename == "no cleavage")
      {
        s += String(indent + 3, '\t') + cv_.getTermByName("NoEnzyme").toXMLString(cv_ns) + "\n";
      }
      else
      {
        s += String(indent + 3, '\t') + cv_.getTermByName("cleavage agent details").toXMLString(cv_ns) + "\n";
      }
      s += String(indent + 2, '\t') + "</EnzymeName>\n";
      s += String(indent + 1, '\t') + "</Enzyme>\n";
      s += String(indent, '\t') + "</Enzymes>\n";
    }

    void MzIdentMLHandler::writeModParam_(String& s, const std::vector<String>& mod_names, bool fixed, UInt indent) const
    {
      String cv_ns = unimod_.name();
      for (std::vector<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
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

            s += String(indent + 1, '\t') + "<SearchModification fixedMod=\"" + (fixed ? "true" : "false") + "\" massDelta=\"" + String((*mt)->getDiffMonoMass()) + "\" residues=\"" + origin + "\">\n";

            // @TODO: handle protein C-term/N-term
            ResidueModification::TermSpecificity spec = (*mt)->getTermSpecificity();
            if ((spec == ResidueModification::C_TERM) || (spec == ResidueModification::N_TERM))
            {
              const String& cv_name = "modification specificity peptide " + (*mt)->getTermSpecificityName();
              s += String(indent + 2, '\t') + "<SpecificityRules>\n";
              s += String(indent + 3, '\t') + cv_.getTermByName(cv_name).toXMLString(cv_ns) + "\n";
              s += String(indent + 2, '\t') + "</SpecificityRules>\n";
            }

            String ac = (*mt)->getUniModAccession();
            if (ac.hasPrefix("UniMod:")) ac = "UNIMOD:" + ac.suffix(':');
            if (!ac.empty())
            {
              s += String(indent + 2, '\t') + unimod_.getTerm(ac).toXMLString(cv_ns) + "\n";
            }
            else
            {
              s += String(indent + 2, '\t') + "<cvParam cvRef=\"MS\" accession=\"MS:1001460\" name=\"unknown modification\"/>\n";
            }
            s += String(indent + 1, '\t') + "</SearchModification>\n";
          }
        }
        else
        {
          String message = String("Registered ") + (fixed ? "fixed" : "variable") + " modification '" + *it + "' is unknown and will be ignored.";
          throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
        }
      }
    }

    void MzIdentMLHandler::writeFragmentAnnotations_(String& s, const std::vector<PeptideHit::PeakAnnotation>& annotations, UInt indent, bool is_ppxl) const
    {
      std::map<UInt,std::map<String,std::vector<StringList> > > annotation_map;
      for (const PeptideHit::PeakAnnotation& pep : annotations)
      {// string coding example: [alpha|ci$y3-H2O-NH3]5+
        // static const boost::regex frag_regex("\\[(?:([\\|\\w]+)\\$)*([abcxyz])(\\d+)((?:[\\+\\-\\w])*)\\](\\d+)\\+"); // this will fetch the complete loss/gain part as one
        static const boost::regex frag_regex_tweak(R"(\[(?:([\|\w]+)\$)*([abcxyz])(\d+)(?:-(H2O|NH3))*\][(\d+)\+]*)"); // this will only fetch the last loss - and is preferred for now, as only these extra cv params are present
        String ionseries_index;
        String iontype;
        //String loss_gain;
        String loss;
        StringList extra;
        boost::smatch str_matches;
        if (boost::regex_match(pep.annotation, str_matches, frag_regex_tweak))
        {
          String(str_matches[1]).split("|",extra);
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
          // OPENMS_LOG_WARN << pep.annotation << std::endl;
          continue;
        }
        String lt = "frag: " + iontype + " ion";
        if (!loss.empty())
        {
          lt += " - "+loss;
        }
        if (annotation_map.find(pep.charge) == annotation_map.end())
        {
          annotation_map[pep.charge] = std::map<String, std::vector<StringList> >();
        }
        if (annotation_map[pep.charge].find(lt) == annotation_map[pep.charge].end())
        {
          annotation_map[pep.charge][lt] = std::vector<StringList> (3);
          if (is_ppxl)
          {
            annotation_map[pep.charge][lt].push_back(StringList());  // alpha|beta
            annotation_map[pep.charge][lt].push_back(StringList());  // ci|xi
          }
        }
        annotation_map[pep.charge][lt][0].push_back(ionseries_index);
        annotation_map[pep.charge][lt][1].push_back(String(pep.mz));
        annotation_map[pep.charge][lt][2].push_back(String(pep.intensity));
        if (is_ppxl)
        {
          String ab = ListUtils::contains<String>(extra ,String("alpha")) ? String("alpha"):String("beta");
          String cx = ListUtils::contains<String>(extra ,String("ci")) ? String("ci"):String("xi");
          annotation_map[pep.charge][lt][3].push_back(ab);
          annotation_map[pep.charge][lt][4].push_back(cx);
        }
      }

      // stop and return, if no mzid compatible fragments were found
      if (annotation_map.empty())
      {
        return;
      }
      //double map: charge + ion type; collect in StringList: index + annotations; write:
      s += String(indent, '\t') + "<Fragmentation>\n";
      for (std::map<UInt,std::map<String,std::vector<StringList> > >::iterator i=annotation_map.begin();
           i!=annotation_map.end(); ++i)
      {
        for (std::map<String,std::vector<StringList> >::iterator j=i->second.begin();
             j!= i->second.end(); ++j)
        {
          s += String(indent+1, '\t') + "<IonType charge=\"" + String(i->first) +"\""
                    + " index=\"" + ListUtils::concatenate(j->second[0], " ") + "\">\n";
          s += String(indent+2, '\t') + "<FragmentArray measure_ref=\"Measure_mz\""
                    + " values=\"" + ListUtils::concatenate(j->second[1], " ") + "\"/>\n";
          s += String(indent+2, '\t') + "<FragmentArray measure_ref=\"Measure_int\""
                    + " values=\"" + ListUtils::concatenate(j->second[2], " ") + "\"/>\n";
          if (is_ppxl)
          {
              s += String(indent+2, '\t') + "<userParam name=\"cross-link_chain\"" + " unitName=\"xsd:string\""
                        + " value=\"" + ListUtils::concatenate(j->second[3], " ") + "\"/>\n";
              s += String(indent+2, '\t') + "<userParam name=\"cross-link_ioncategory\"" + " unitName=\"xsd:string\""
                        + " value=\"" + ListUtils::concatenate(j->second[4], " ") + "\"/>\n";
          }
          s += String(indent+2, '\t') + cv_.getTermByName(j->first).toXMLString("PSI-MS") + "\n";
          s += String(indent+1, '\t') + "</IonType>\n";
        }
      }
      s += String(indent, '\t') + "</Fragmentation>\n";
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

    String MzIdentMLHandler::trimOpenMSfileURI(const String& file) const
    {
      String r = file;
      if (r.hasPrefix("["))
        r = r.substr(1);
      if (r.hasSuffix("]"))
        r = r.substr(0,r.size()-1);
      r.substitute("\\","/");
      return r;
    }

    void MzIdentMLHandler::writePeptideHit(const PeptideHit& hit,
                                                std::vector<PeptideIdentification>::const_iterator& it,
                                                std::map<String, String>& pep_ids,
                                                const String& cv_ns, std::set<String>& sen_set,
                                                std::map<String, String>& sen_ids,
                                                std::map<String, std::vector<String> >& pep_evis,
                                                std::map<String, double>& pp_identifier_2_thresh,
                                                String& sidres)
    {
        String pepid =  "PEP_" + String(UniqueIdGenerator::getUniqueId());
        String pepi = hit.getSequence().toString();

        bool duplicate = false;
        std::map<String, String>::iterator pit;

        // avoid duplicates
        pit = pep_ids.find(pepi);
        if (pit != pep_ids.end())
        {
          duplicate = true;
        }

        if (!duplicate)
        {
          String p;
          //~ TODO simplify mod cv param write
          // write peptide id with conversion to universal, "human-readable" bracket string notation
          p += String("\t<Peptide id=\"") + pepid + String("\" name=\"") +
                hit.getSequence().toBracketString(false) + String("\">\n\t\t<PeptideSequence>") + hit.getSequence().toUnmodifiedString() + String("</PeptideSequence>\n");
          if (hit.getSequence().isModified())
          {
            const ResidueModification* n_term_mod = hit.getSequence().getNTerminalModification();
            if (n_term_mod != nullptr)
            {
              p += "\t\t<Modification location=\"0\">\n";
              String acc = n_term_mod->getUniModAccession();
              p += "\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':');
              p += "\" name=\"" + n_term_mod->getId();
              p += R"(" cvRef="UNIMOD"/>)";
              p += "\n\t\t</Modification>\n";
            }
            const ResidueModification* c_term_mod = hit.getSequence().getCTerminalModification();
            if (c_term_mod != nullptr)
            {
              p += "\t\t<Modification location=\"" + String(hit.getSequence().size()) + "\">\n";
              String acc = c_term_mod->getUniModAccession();
              p += "\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':');
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
                p += "\t\t<Modification location=\"" + String(i + 1);
                p += "\" residues=\"" + hit.getSequence()[i].getOneLetterCode();
                String acc = mod->getUniModAccession();
                if (!acc.empty())
                {
                  p += "\">\n\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':'); //TODO @all: do not exclusively use unimod ...
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
                    p += "\" monoisotopicMassDelta=\"" + String(diffmass);
                  }
                  else if (mod->getMonoMass() > 0.0)
                  {
                    double diffmass = mod->getMonoMass() - hit.getSequence()[i].getMonoWeight();
                    p += "\" monoisotopicMassDelta=\"" + String(diffmass);
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

        std::vector<String> pevid_ids;
        if (!duplicate)
        {
          std::vector<PeptideEvidence> peptide_evidences = hit.getPeptideEvidences();
          // TODO idXML allows peptide hits without protein references! Fails in that case - run PeptideIndexer first
          for (std::vector<PeptideEvidence>::const_iterator pe = peptide_evidences.begin(); pe != peptide_evidences.end(); ++pe)
          {
            String pevid =  "PEV_" + String(UniqueIdGenerator::getUniqueId());
            String dBSequence_ref;
            map<String, String>::const_iterator pos = sen_ids.find(pe->getProteinAccession());
            if (pos != sen_ids.end())
            {
              dBSequence_ref = pos->second;
            }
            else
            {
              OPENMS_LOG_ERROR << "Error: Missing or invalid protein reference for peptide '" << pepi << "': '" << pe->getProteinAccession() << "' - skipping." << endl;
              continue;
            }
            String idec;
            if (hit.metaValueExists(Constants::UserParam::TARGET_DECOY))
            {
              idec = String(boost::lexical_cast<std::string>((String(hit.getMetaValue(Constants::UserParam::TARGET_DECOY))).hasSubstring("decoy")));
            }

            String e;
            String nc_termini = "-";    // character for N- and C-termini as specified in mzIdentML
            e += "\t<PeptideEvidence id=\"" + pevid + "\" peptide_ref=\"" + pepid + "\" dBSequence_ref=\"" + dBSequence_ref + "\"";

            if (pe->getAAAfter() != PeptideEvidence::UNKNOWN_AA)
            {
              e += " post=\"" + (pe->getAAAfter() == PeptideEvidence::C_TERMINAL_AA ? nc_termini : String(pe->getAAAfter())) + "\"";
            }
            if (pe->getAABefore() != PeptideEvidence::UNKNOWN_AA)
            {
              e += " pre=\"" + (pe->getAABefore() == PeptideEvidence::N_TERMINAL_AA ? nc_termini : String(pe->getAABefore())) + "\"";
            }
            if (pe->getStart() != PeptideEvidence::UNKNOWN_POSITION)
            {
              e += " start=\"" + String(pe->getStart() + 1) + "\"";
            }
            else if (hit.metaValueExists("start"))
            {
              e += " start=\"" + String( int(hit.getMetaValue("start")) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no start position of peptide hit in protein sequence." << std::endl;
            }
            if (pe->getEnd() != PeptideEvidence::UNKNOWN_POSITION)
            {
              e += " end=\"" + String(pe->getEnd() + 1) + "\"";
            }
            else if (hit.metaValueExists("end"))
            {
              e += " end=\"" + String( int(hit.getMetaValue("end")) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no end position of peptide hit in protein sequence." << std::endl;
            }
            if (!idec.empty())
            {
              e += " isDecoy=\"" + String(idec)+ "\"";
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

        String cmz(hit.getSequence().getMZ(hit.getCharge())); //calculatedMassToCharge
        String r(hit.getRank()); //rank
        String sc(hit.getScore());

        if (sc.empty())
        {
          sc = "NA";
          OPENMS_LOG_WARN << "No score assigned to this PSM: " /*<< hit.getSequence().toString()*/ << std::endl;
        }
        String c(hit.getCharge()); //charge

        String pte;
        if (hit.metaValueExists("pass_threshold"))
        {
          pte = boost::lexical_cast<std::string>(hit.getMetaValue("pass_threshold"));
        }
        else if (pp_identifier_2_thresh.find(it->getIdentifier())!= pp_identifier_2_thresh.end() && pp_identifier_2_thresh.find(it->getIdentifier())->second != 0.0)
        {
          double th = pp_identifier_2_thresh.find(it->getIdentifier())->second;
          //threshold was 'set' in proteinIdentification (!= default value of member, now check pass
          pte = boost::lexical_cast<std::string>(it->isHigherScoreBetter() ? hit.getScore() > th : hit.getScore() < th); //passThreshold-eval
        }
        else
        {
          pte = true;
        }

        //write SpectrumIdentificationItem elements
        String emz(it->getMZ());
        String sii = "SII_" + String(UniqueIdGenerator::getUniqueId());
        String sii_tmp;
        sii_tmp += String("\t\t\t\t<SpectrumIdentificationItem passThreshold=\"")
                + pte + String("\" rank=\"") + r + String("\" peptide_ref=\"")
                + pepid + String("\" calculatedMassToCharge=\"") + cmz
                + String("\" experimentalMassToCharge=\"") + emz
                + String("\" chargeState=\"") + c +  String("\" id=\"")
                + sii + String("\">\n");

        if (pevid_ids.empty())
        {
          OPENMS_LOG_WARN << "PSM without peptide evidence registered in the given search database found. This will cause an invalid mzIdentML file (which OpenMS can still consume)." << std::endl;
        }
        for (std::vector<String>::const_iterator pevref = pevid_ids.begin(); pevref != pevid_ids.end(); ++pevref)
        {
          sii_tmp += "\t\t\t\t\t<PeptideEvidenceRef peptideEvidence_ref=\"" +  String(*pevref) + "\"/>\n";
        }

        if (! hit.getPeakAnnotations().empty())
        {
          writeFragmentAnnotations_(sii_tmp, hit.getPeakAnnotations(), 5, false);
        }

        std::set<String> peptide_result_details;
        cv_.getAllChildTerms(peptide_result_details, "MS:1001143"); // search engine specific score for PSMs
        MetaInfoInterface copy_hit = hit;
        String st(it->getScoreType()); //scoretype

        if (cv_.hasTermWithName(st) && peptide_result_details.find(cv_.getTermByName(st).id) != peptide_result_details.end())
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName(st).toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTermByName(st).id);
        }
        else if (cv_.exists(st) && peptide_result_details.find(st) != peptide_result_details.end())
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTerm(st).toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTerm(st).id);
        }
        else if (st == "q-value" || st == "FDR")
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("PSM-level q-value").toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTermByName("PSM-level q-value").id);
        }
        else if (st == "Posterior Error Probability")
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("percolator:PEP").toXMLString(cv_ns, sc); // 'percolaror' was not a typo in the code but in the cv.
          copy_hit.removeMetaValue(cv_.getTermByName("percolator:PEP").id);
        }
        else if (st == "OMSSA")
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("OMSSA:evalue").toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTermByName("OMSSA:evalue").id);
        }
        else if (st == "Mascot")
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("Mascot:score").toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTermByName("Mascot:score").id);
        }
        else if (st == "XTandem")
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("X\\!Tandem:hyperscore").toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTermByName("X\\!Tandem:hyperscore").id);
        }
        else if (st == "SEQUEST")
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("Sequest:xcorr").toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTermByName("Sequest:xcorr").id);
        }
        else if (st == "MS-GF+")
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("MS-GF:RawScore").toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTermByName("MS-GF:RawScore").id);
        }
        else if (st == Constants::UserParam::OPENPEPXL_SCORE)
        {
          sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName(st).toXMLString(cv_ns, sc);
          copy_hit.removeMetaValue(cv_.getTermByName(st).id);
        }
        else
        {
          String score_name_placeholder = st.empty()?"PSM-level search engine specific statistic":st;
          sii_tmp += String(5, '\t') + cv_.getTermByName("PSM-level search engine specific statistic").toXMLString(cv_ns);
          sii_tmp += "\n" + String(5, '\t') + "<userParam name=\"" + score_name_placeholder
                       + "\" unitName=\"" + "xsd:double" + "\" value=\"" + sc + "\"/>";
          OPENMS_LOG_WARN << "Converting unknown score type to PSM-level search engine specific statistic from PSI controlled vocabulary." << std::endl;
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
                                                std::vector<PeptideIdentification>::const_iterator& it,
                                                const String& ppxl_linkid, std::map<String, String>& pep_ids,
                                                const String& cv_ns, std::set<String>& sen_set,
                                                std::map<String, String>& sen_ids,
                                                std::map<String, std::vector<String> >& pep_evis,
                                                std::map<String, double>& pp_identifier_2_thresh,
                                                double ppxl_crosslink_mass,
                                                std::map<String, String>& ppxl_specref_2_element,
                                                String& sid, bool alpha_peptide)
    {

      String pepid =  "PEP_" + String(UniqueIdGenerator::getUniqueId());

      AASequence peptide_sequence;
      if (alpha_peptide)
      {
        peptide_sequence = hit.getSequence();
      }
      else
      {
        peptide_sequence = AASequence::fromString(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE));
      }
      String pepi = peptide_sequence.toString();

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
      std::map<String, String>::iterator pit;
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
      //            String pepi2 = peps[1].getSequence().toString();
      //            std::map<String, String>::iterator pit = pep_pairs_ppxl.find(pepi);

      //            // last entry in vector
      //            std::pair<String, String> pepid_pairs_ppxl[pepid_pairs_ppxl.size()-1];

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
        String p;
        //~ TODO simplify mod cv param write
        // write peptide id with conversion to universal, "human-readable" bracket string notation
        p += String("\t<Peptide id=\"") + pepid + String("\" name=\"") +
              peptide_sequence.toBracketString(false) + String("\">\n\t\t<PeptideSequence>") + peptide_sequence.toUnmodifiedString() + String("</PeptideSequence>\n");

        const ResidueModification* n_term_mod = peptide_sequence.getNTerminalModification();
        if (n_term_mod != nullptr)
        {
          p += "\t\t<Modification location=\"0\">\n";
          String acc = n_term_mod->getUniModAccession();
          bool unimod = true;
          if (!acc.empty())
          {
            p += "\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':');
          }
          else
          {
            acc = n_term_mod->getPSIMODAccession();
            p += "\t\t\t<cvParam accession=\"XLMOD:" + acc.suffix(':');
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
          p += "\t\t<Modification location=\"" + String(peptide_sequence.size()) + "\">\n";
          String acc = c_term_mod->getUniModAccession();
          bool unimod = true;
          if (!acc.empty())
          {
            p += "\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':');
          }
          else
          {
            acc = c_term_mod->getPSIMODAccession();
            p += "\t\t\t<cvParam accession=\"XLMOD:" + acc.suffix(':');
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
            p += "\t\t<Modification location=\"" + String(i + 1);
            p += "\" residues=\"" + peptide_sequence[i].getOneLetterCode();
            String acc = mod->getUniModAccession();
            if (!acc.empty())
            {
              p += "\">\n\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':'); //TODO @all: do not exclusively use unimod ...
              p += "\" name=\"" + mod->getId();
              p += R"(" cvRef="UNIMOD"/>)";
              p += "\n\t\t</Modification>\n";
            }
            else
            {
              acc = mod->getPSIMODAccession();
              if (!acc.empty())
              {
                p += "\">\n\t\t\t<cvParam accession=\"XLMOD:" + acc.suffix(':');
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
                  p += "\" monoisotopicMassDelta=\"" + String(diffmass);
                }
                else if (mod->getMonoMass() > 0.0)
                {
                  double diffmass = mod->getMonoMass() - peptide_sequence[i].getMonoWeight();
                  p += "\" monoisotopicMassDelta=\"" + String(diffmass);
                }
                p += "\">\n\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1001460\" name=\"unknown modification\"/>";
                p += "\n\t\t</Modification>\n";
              }
            }
          }
          // Failsafe, if someone uses a new cross-linker (given these conditions, there MUST be a linker at this position, but it does not have a Unimod or XLMOD entry)
          else if (alpha_peptide && hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_MOD) && (static_cast<Size>(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1).toString().toInt()) == i) && (hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE).toString() == "mono-link") )
          {
            p += "\t\t<Modification location=\"" + String(i + 1);
            p += "\" residues=\"" + String(peptide_sequence[i].getOneLetterCode());
            p += "\">\n\t\t\t<cvParam accession=\"XLMOD:XXXXX";
            p += "\" name=\"" +  hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD).toString();
            p += R"(" cvRef="XLMOD"/>)";
            p += "\n\t\t</Modification>\n";
          }
        }

        if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TYPE) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) != "mono-link")
        {
          int i = hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1).toString().toInt();
          if (alpha_peptide)
          {
            CrossLinksDB* xl_db = CrossLinksDB::getInstance();
            std::vector<String> mods;
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
                p += "\t\t<Modification location=\"" + String(i + 2);
              }
            }
            else
            {
              xl_db->searchModificationsByDiffMonoMass(mods, double(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS)), 0.0001, String(hit.getSequence()[i].getOneLetterCode()), ResidueModification::ANYWHERE);
              if (!mods.empty())
              {
                p += "\t\t<Modification location=\"" + String(i + 1);
              }
            }

            String acc = "";
            String name = "";
            for (Size s = 0; s < mods.size(); ++s)
            {
              if (mods[s].hasSubstring(hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD)))
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
              const ResidueModification* mod = xl_db->getModification( String(peptide_sequence[i].getOneLetterCode()), mods[0], ResidueModification::ANYWHERE);
              acc = mod->getPSIMODAccession();
              name = mod->getId();
            }
            if (!acc.empty())
            {
              p += "\" residues=\"" + String(peptide_sequence[i].getOneLetterCode());
              p += "\" monoisotopicMassDelta=\"" + hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS).toString() + "\">\n";
              p += "\t\t\t<cvParam accession=\"" + acc + R"(" cvRef="XLMOD" name=")" + name + "\"/>\n";
            }
            else // if there is no matching modification in the database, write out a placeholder
            {
              p += "\t\t<Modification location=\"" + String(i + 1);
              p += "\" residues=\"" + String(peptide_sequence[i].getOneLetterCode());
              p += "\" monoisotopicMassDelta=\"" + hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS).toString() + "\">\n";
              p += "\t\t\t<cvParam accession=\"XLMOD:XXXXX\" cvRef=\"XLMOD\" name=\"" + hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD).toString() + "\"/>\n";
            }
          }
          else // xl_chain = "MS:1002510", acceptor, beta peptide
          {
            i = hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2).toString().toInt();
            if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA) == "N_TERM")
            {
              p += "\t\t<Modification location=\"0";
            }
            else if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA) == "C_TERM")
            {
              p += "\t\t<Modification location=\"" + String(peptide_sequence.size() + 2);
            }
            else
            {
              p += "\t\t<Modification location=\"" + String(i + 1);
            }
            p += "\" residues=\"" + String(peptide_sequence[i].getOneLetterCode());
            p += "\" monoisotopicMassDelta=\"0\">\n";
          }
          if (alpha_peptide)
          {
            p += "\t\t\t" + cv_.getTerm(String("MS:1002509")).toXMLString(cv_ns, DataValue(ppxl_linkid));
          }
          else
          {
            p += "\t\t\t" + cv_.getTerm(String("MS:1002510")).toXMLString(cv_ns, DataValue(ppxl_linkid));
          }
          p += "\n\t\t</Modification>\n";
        }
        if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TYPE) && hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "loop-link")
        {
          int i = hit.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2).toString().toInt();
          p += "\t\t<Modification location=\"" + String(i + 1);
          p += "\" residues=\"" + String(peptide_sequence[i].getOneLetterCode());
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

      std::vector<String> pevid_ids;
      if (!duplicate)
      {
        if (alpha_peptide)
        {
          std::vector<PeptideEvidence> peptide_evidences = hit.getPeptideEvidences();
          // TODO idXML allows peptide hits without protein references! Fails in that case - run PeptideIndexer first

          // TODO BETA PEPTIDE Protein Info
          for (std::vector<PeptideEvidence>::const_iterator pe = peptide_evidences.begin(); pe != peptide_evidences.end(); ++pe)
          {
            String pevid =  "PEV_" + String(UniqueIdGenerator::getUniqueId());
            String dBSequence_ref;
            map<String, String>::const_iterator pos = sen_ids.find(pe->getProteinAccession());
            if (pos != sen_ids.end())
            {
              dBSequence_ref = pos->second;
            }
            else
            {
              OPENMS_LOG_ERROR << "Error: Missing or invalid protein reference for peptide '" << pepi << "': '" << pe->getProteinAccession() << "' - skipping." << endl;
              continue;
            }
            String idec;
            if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA))
            {
              idec = String(boost::lexical_cast<std::string>((String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA))).hasSubstring("decoy")));
            }

            String e;
            String nc_termini = "-";    // character for N- and C-termini as specified in mzIdentML
            e += "\t<PeptideEvidence id=\"" + pevid + "\" peptide_ref=\"" + pepid + "\" dBSequence_ref=\"" + dBSequence_ref + "\"";

            if (pe->getAAAfter() != PeptideEvidence::UNKNOWN_AA)
            {
              e += " post=\"" + (pe->getAAAfter() == PeptideEvidence::C_TERMINAL_AA ? nc_termini : String(pe->getAAAfter())) + "\"";
            }
            if (pe->getAABefore() != PeptideEvidence::UNKNOWN_AA)
            {
              e += " pre=\"" + (pe->getAABefore() == PeptideEvidence::N_TERMINAL_AA ? nc_termini : String(pe->getAABefore())) + "\"";
            }
            if (pe->getStart() != PeptideEvidence::UNKNOWN_POSITION)
            {
              e += " start=\"" + String(pe->getStart() + 1) + "\"";
            }
            else if (hit.metaValueExists("start"))
            {
              e += " start=\"" + String( int(hit.getMetaValue("start")) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no start position of peptide hit in protein sequence." << std::endl;
            }
            if (pe->getEnd() != PeptideEvidence::UNKNOWN_POSITION)
            {
              e += " end=\"" + String(pe->getEnd() + 1) + "\"";
            }
            else if (hit.metaValueExists("end"))
            {
              e += " end=\"" + String( int(hit.getMetaValue("end")) + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no end position of peptide hit in protein sequence." << std::endl;
            }
            if (!idec.empty())
            {
              e += " isDecoy=\"" + String(idec)+ "\"";
            }
            e += "/>\n";
            sen_set.insert(e);
            pevid_ids.push_back(pevid);
          }
          pep_evis.insert(make_pair(pepi, pevid_ids));
        }
        else // acceptor, beta peptide, does not have its own PeptideHit and PeptideEvidences
        {
          StringList prot = ListUtils::create<String>(String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS)));
          StringList pre = ListUtils::create<String>(String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_PRE)));
          StringList post = ListUtils::create<String>(String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_POST)));
          StringList start = ListUtils::create<String>(String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_START)));
          StringList end = ListUtils::create<String>(String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_END)));
          for (Size ev = 0; ev < pre.size(); ++ev)
          {
            String pevid =  "PEV_" + String(UniqueIdGenerator::getUniqueId());
            String dBSequence_ref;
            map<String, String>::const_iterator pos = sen_ids.find(prot[ev]);
            if (pos != sen_ids.end())
            {
              dBSequence_ref = pos->second;
            }
            else
            {
              OPENMS_LOG_ERROR << "Error: Missing or invalid protein reference for peptide '" << pepi << "': '" << prot[ev] << "' - skipping." << endl;
              continue;
            }
            String idec;
            if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA))
            {
              idec = String(boost::lexical_cast<std::string>((String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA))).hasSubstring("decoy")));
            }

            String e;
            String nc_termini = "-";    // character for N- and C-termini as specified in mzIdentML
            e += "\t<PeptideEvidence id=\"" + pevid + "\" peptide_ref=\"" + pepid + "\" dBSequence_ref=\"" + dBSequence_ref + "\"";

            if (post[ev] != String(PeptideEvidence::UNKNOWN_AA))
            {
              e += " post=\"" + (post[ev] == String(PeptideEvidence::C_TERMINAL_AA) ? nc_termini : post[ev]) + "\"";
            }
            if (pre[ev] != String(PeptideEvidence::UNKNOWN_AA))
            {
              e += " pre=\"" + (pre[ev] == String(PeptideEvidence::N_TERMINAL_AA) ? nc_termini : pre[ev]) + "\"";
            }
            if (start[ev] != String(PeptideEvidence::UNKNOWN_POSITION))
            {
              e += " start=\"" + String(start[ev].toInt() + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no start position of peptide hit in protein sequence." << std::endl;
            }
            if (end[ev] != String(PeptideEvidence::UNKNOWN_POSITION))
            {
              e += " end=\"" + String(end[ev].toInt() + 1) + "\"";
            }
            else
            {
              OPENMS_LOG_WARN << "Found no end position of peptide hit in protein sequence." << std::endl;
            }
            if (!idec.empty())
            {
              e += " isDecoy=\"" + String(idec)+ "\"";
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

      String r(hit.getRank()); //rank
      String sc(hit.getScore());
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
      String cmz = String( ((calc_ppxl_mass) +  (hit.getCharge() * Constants::PROTON_MASS_U)) / hit.getCharge()); //calculatedMassToCharge

      if (sc.empty())
      {
        sc = "NA";
        OPENMS_LOG_WARN << "No score assigned to this PSM: " /*<< hit.getSequence().toString()*/ << std::endl;
      }
      String c(hit.getCharge()); //charge

      String pte;
      if (hit.metaValueExists("pass_threshold"))
      {
        pte = boost::lexical_cast<std::string>(hit.getMetaValue("pass_threshold"));
      }
      else if (pp_identifier_2_thresh.find(it->getIdentifier())!= pp_identifier_2_thresh.end() && pp_identifier_2_thresh.find(it->getIdentifier())->second != 0.0)
      {
        double th = pp_identifier_2_thresh.find(it->getIdentifier())->second;
        //threshold was 'set' in proteinIdentification (!= default value of member, now check pass
        pte = boost::lexical_cast<std::string>(it->isHigherScoreBetter() ? hit.getScore() > th : hit.getScore() < th); //passThreshold-eval
      }
      else
      {
        pte = true;
      }

      //write SpectrumIdentificationItem elements
      String emz(it->getMZ());
      String sii = "SII_" + String(UniqueIdGenerator::getUniqueId());
      String sii_tmp;
      sii_tmp += String("\t\t\t\t<SpectrumIdentificationItem passThreshold=\"")
              + pte + String("\" rank=\"") + r + String("\" peptide_ref=\"")
              + pepid + String("\" calculatedMassToCharge=\"") + cmz
              + String("\" experimentalMassToCharge=\"") + emz
              + String("\" chargeState=\"") + c +  String("\" id=\"")
              + sii + String("\">\n");

      if (pevid_ids.empty())
      {
        OPENMS_LOG_WARN << "PSM without peptide evidence registered in the given search database found. This will cause an invalid mzIdentML file (which OpenMS can still consume)." << std::endl;
      }
      for (std::vector<String>::const_iterator pevref = pevid_ids.begin(); pevref != pevid_ids.end(); ++pevref)
      {
        sii_tmp += "\t\t\t\t\t<PeptideEvidenceRef peptideEvidence_ref=\"" +  String(*pevref) + "\"/>\n";
      }

      if (! hit.getPeakAnnotations().empty() && alpha_peptide)
      {
        // is_ppxl = true
        writeFragmentAnnotations_(sii_tmp, hit.getPeakAnnotations(), 5, true);
      }

      std::set<String> peptide_result_details;
      cv_.getAllChildTerms(peptide_result_details, "MS:1001143"); // search engine specific score for PSMs
      MetaInfoInterface copy_hit = hit;
      String st(it->getScoreType()); //scoretype

      if (cv_.hasTermWithName(st) && peptide_result_details.find(cv_.getTermByName(st).id) != peptide_result_details.end())
      {
        sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName(st).toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(cv_.getTermByName(st).id);
      }
      else if (cv_.exists(st) && peptide_result_details.find(st) != peptide_result_details.end())
      {
        sii_tmp +=  "\t\t\t\t\t" + cv_.getTerm(st).toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(cv_.getTerm(st).id);
      }
      else if (st == "q-value" || st == "FDR")
      {
        sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("PSM-level q-value").toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(cv_.getTermByName("PSM-level q-value").id);
      }
      else if (st == "Posterior Error Probability")
      {
        sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName("percolator:PEP").toXMLString(cv_ns, sc); // 'percolaror' was not a typo in the code but in the cv.
        copy_hit.removeMetaValue(cv_.getTermByName("percolator:PEP").id);
      }
      else if (st == Constants::UserParam::OPENPEPXL_SCORE)
      {
        sii_tmp +=  "\t\t\t\t\t" + cv_.getTermByName(st).toXMLString(cv_ns, sc);
        copy_hit.removeMetaValue(cv_.getTermByName(st).id);
      }
      else
      {
        String score_name_placeholder = st.empty()?"PSM-level search engine specific statistic":st;
        sii_tmp += String(5, '\t') + cv_.getTermByName("PSM-level search engine specific statistic").toXMLString(cv_ns);
        sii_tmp += "\n" + String(5, '\t') + "<userParam name=\"" + score_name_placeholder
                     + "\" unitName=\"" + "xsd:double" + "\" value=\"" + sc + "\"/>";
        OPENMS_LOG_WARN << "Converting unknown score type to PSM-level search engine specific statistic from PSI controlled vocabulary." << std::endl;
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
      String ert = rt == rt ? String(rt) : "nan";
      DataValue rtcv(ert);
      rtcv.setUnit(10); // id: UO:0000010 name: second
      rtcv.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
      sii_tmp = sii_tmp.substitute("</SpectrumIdentificationItem>",
                                   "\t" + cv_.getTermByName("retention time").toXMLString(cv_ns, rtcv) + "\n\t\t\t\t</SpectrumIdentificationItem>\n");
      ppxl_specref_2_element[sid] += sii_tmp;
      if (hit.metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT) && hit.metaValueExists(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ))
      {
        // ppxl : TODO remove if existing <Fragmentation/>block
        sii_tmp = sii_tmp.substitute(String("experimentalMassToCharge=\"") + String(emz),
                                     String("experimentalMassToCharge=\"") + String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ))); // mz
        sii_tmp = sii_tmp.substitute(sii, String("SII_") + String(UniqueIdGenerator::getUniqueId())); // uid
        sii_tmp = sii_tmp.substitute("value=\"" + ert, "value=\"" + String(hit.getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT)));

        ProteinIdentification::SearchParameters search_params = cpro_id_->front().getSearchParameters();
        double iso_shift = String(search_params.getMetaValue("cross_link:mass_isoshift")).toDouble();
        double cmz_heavy = cmz.toDouble() + (iso_shift / hit.getCharge());

        sii_tmp = sii_tmp.substitute(String("calculatedMassToCharge=\"") + String(cmz),
                                      String("calculatedMassToCharge=\"") + String(cmz_heavy));

        ppxl_specref_2_element[sid] += sii_tmp;
      }
    }
  } // namespace OpenMS
