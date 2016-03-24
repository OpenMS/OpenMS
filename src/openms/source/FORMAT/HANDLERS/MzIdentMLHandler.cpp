// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzIdentMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Enzyme.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <set>

#include <boost/lexical_cast.hpp>


using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    MzIdentMLHandler::MzIdentMLHandler(const Identification& id, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      //~ ms_exp_(0),
      id_(0),
      cid_(&id)
    {
      cv_.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      unimod_.loadFromOBO("PSI-MS", File::find("/CV/unimod.obo"));
    }

    MzIdentMLHandler::MzIdentMLHandler(Identification& id, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      //~ ms_exp_(0),
      id_(&id),
      cid_(0)
    {
      cv_.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      unimod_.loadFromOBO("PSI-MS", File::find("/CV/unimod.obo"));
    }

    MzIdentMLHandler::MzIdentMLHandler(const std::vector<ProteinIdentification>& pro_id, const std::vector<PeptideIdentification>& pep_id, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      //~ ms_exp_(0),
      pro_id_(0),
      pep_id_(0),
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
      cpro_id_(0),
      cpep_id_(0)
    {
      cv_.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
      unimod_.loadFromOBO("PSI-MS", File::find("/CV/unimod.obo"));
    }

    //~ TODO create MzIdentML instances from MSExperiment which contains much of the information yet needed
    //~ MzIdentMLHandler(const MSExperiment<>& mx, const String& filename, const String& version, const ProgressLogger& logger)
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
    {
    }

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
      error(LOAD, "MzIdentMLHandler::startElement: Unkown element found: '" + tag_ + "' in tag '" + parent_tag + "', ignoring.");
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

      //error(LOAD, "MzIdentMLHandler::characters: Unkown character section found: '" + tag_ + "', ignoring.");
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
      error(LOAD, "MzIdentMLHandler::endElement: Unkown element found: '" + tag_ + "', ignoring.");
    }

    void MzIdentMLHandler::handleCVParam_(const String& /* parent_parent_tag*/, const String& parent_tag, const String& accession, /* const String& name, */ /* const String& value, */ const xercesc::Attributes& attributes, const String& cv_ref /* , const String& unit_accession */)
    {
      if (parent_tag == "Modification")
      {
        if (cv_ref == "UNIMOD")
        {
          //void ModificationsDB::searchModifications(set<const ResidueModification*>& mods, const String& origin, const String& name, ResidueModification::Term_Specificity term_spec) const
          set<const ResidueModification*> mods;
          Int loc = numeric_limits<Size>::max();
          if (optionalAttributeAsInt_(loc, attributes, "location"))
          {
            String uni_mod_id = accession.suffix(':');
            // TODO handle ambiguous residues
            String residues;
            if (optionalAttributeAsString_(residues, attributes, "residues"))
            {

            }
            if (loc == 0)
            {
              ModificationsDB::getInstance()->searchTerminalModifications(mods, uni_mod_id, ResidueModification::N_TERM);
            }
            else if (loc == (Int)actual_peptide_.size())
            {
              ModificationsDB::getInstance()->searchTerminalModifications(mods, uni_mod_id, ResidueModification::C_TERM);
            }
            else
            {
              ModificationsDB::getInstance()->searchModifications(mods, residues, uni_mod_id, ResidueModification::ANYWHERE);
            }
          }
          else
          {
            warning(LOAD, "location of modification not defined!");
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
      std::map<String, double> pp_identifier_2_thresh;

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
        else
        {
          osecv = "analysis software";
        }

        if (soit == sof_ids.end())
        {
          sof_id = "SOF_" + String(UniqueIdGenerator::getUniqueId());
          //~ TODO consider not only searchengine but also version!
          String sost = String("\t<AnalysisSoftware version=\"") + String(it->getSearchEngineVersion()) + String("\" name=\"") + sof_name +  String("\" id=\"") + sof_id + String("\"> \n") + String("\t\t<SoftwareName> \n ");
          sost += "\t\t\t" + cv_.getTermByName(osecv).toXMLString(cv_ns);
          sost += String("\n\t\t</SoftwareName> \n\t</AnalysisSoftware> \n");
          sof_set.insert(sost);
          sof_ids.insert(make_pair(sof_name, sof_id));
        }
        else
        {
          sof_id = soit->second;
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

        String sip = String("\t<SpectrumIdentificationProtocol id=\"") + String(sip_id) + String("\" analysisSoftware_ref=\"") + String(sof_id) + String("\">");
        sip += String(" \n\t\t<SearchType>\n\t\t\t") + cv_.getTermByName("ms-ms search").toXMLString(cv_ns) + String(" \n\t\t</SearchType>");
        sip += String("\n\t\t<AdditionalSearchParams>\n");
        writeMetaInfos_(sip, it->getSearchParameters(), 3);
        sip += String(3, '\t') + "<userParam name=\"" + "charges" + "\" unitName=\"" + "xsd:string" + "\" value=\"" + it->getSearchParameters().charges + "\"/>" + "\n";
//        sip += String(3, '\t') + "<userParam name=\"" + "missed_cleavages" + "\" unitName=\"" + "xsd:integer" + "\" value=\"" + String(it->getSearchParameters().missed_cleavages) + "\"/>" + "\n";
        sip += String("\t\t</AdditionalSearchParams>\n");
        writeModParam_(sip, it->getSearchParameters().fixed_modifications, it->getSearchParameters().variable_modifications, 2);
        writeEnzyme_(sip, it->getSearchParameters().digestion_enzyme, it->getSearchParameters().missed_cleavages, 2);
        // TODO MassTable section
        sip += String("\t\t<FragmentTolerance>\n");
        String unit_str = "unitCvRef=\"UO\" unitName=\"dalton\" unitAccession=\"UO:0000221\"";
        if (it->getSearchParameters().fragment_mass_tolerance_ppm)
        {
          unit_str = "unitCvRef=\"UO\" unitName=\"parts per million\" unitAccession=\"UO:0000169\"";
        }
        sip += String(3, '\t') + "<cvParam accession=\"MS:1001412\" name=\"search tolerance plus value\" " + unit_str + " cvRef=\"PSI-MS\" value=\"" + String(it->getSearchParameters().fragment_mass_tolerance) + "\"/>" + "\n";
        sip += String(3, '\t') + "<cvParam accession=\"MS:1001413\" name=\"search tolerance minus value\" " + unit_str + " cvRef=\"PSI-MS\" value=\"" + String(it->getSearchParameters().fragment_mass_tolerance) + "\"/>" + "\n";
        sip += String("\t\t</FragmentTolerance>\n");
        sip += String("\t\t<ParentTolerance>\n");
        unit_str = "unitCvRef=\"UO\" unitName=\"dalton\" unitAccession=\"UO:0000221\"";
        if (it->getSearchParameters().precursor_mass_tolerance_ppm)
        {
          unit_str = "unitCvRef=\"UO\" unitName=\"parts per million\" unitAccession=\"UO:0000169\"";
        }
        sip += String(3, '\t') + "<cvParam accession=\"MS:1001412\" name=\"search tolerance plus value\" " + unit_str + " cvRef=\"PSI-MS\" value=\"" + String(it->getSearchParameters().precursor_tolerance) + "\"/>" + "\n";
        sip += String(3, '\t') + "<cvParam accession=\"MS:1001413\" name=\"search tolerance minus value\" " + unit_str + " cvRef=\"PSI-MS\" value=\"" + String(it->getSearchParameters().precursor_tolerance) + "\"/>" + "\n";
        sip += String("\t\t</ParentTolerance>\n");
        sip += String("\t\t<Threshold>\n\t\t\t") + thcv + "\n";
        sip += String("\t\t</Threshold>\n");
        sip += String("\t</SpectrumIdentificationProtocol>\n");
        sip_set.insert(sip);
        sil_2_date.insert(make_pair(sil_id, String(it->getDateTime().getDate() + "T" + it->getDateTime().getTime())));


        //~ collect SpectraData element for each ProteinIdentification
        String sdat_id;
        StringList sdat_files;
        String sdat_file("UNKNOWN");
        
        if (it->metaValueExists("spectra_data"))
        {
          sdat_files = it->getMetaValue("spectra_data");
          if (!sdat_files.empty() && !sdat_files[0].empty())
          {
            sdat_file = sdat_files[0];
          }
        }

        std::map<String, String>::iterator sdit = sdat_ids.find(sdat_file); //this part is strongly connected to AnalysisCollection write part
        if (sdit == sdat_ids.end())
        {
          sdat_id = "SDAT_" + String(UniqueIdGenerator::getUniqueId());

          //xml
          spectra_data += String("\t\t<SpectraData location=\"") + sdat_file + String("\" id=\"") + sdat_id + String("\">");
          spectra_data += String("\n\t\t\t<FileFormat> \n");
          spectra_data += String(4, '\t') + cv_.getTermByName("mzML format").toXMLString(cv_ns);
          spectra_data += String("\n\t\t\t</FileFormat>\n\t\t\t<SpectrumIDFormat> \n ");
          spectra_data += String(4, '\t') + cv_.getTermByName("multiple peak list nativeID format").toXMLString(cv_ns);
          spectra_data += String("\n\t\t\t</SpectrumIDFormat> \n\t\t</SpectraData>\n");

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
        String sdb_file(it->getSearchParameters().db); //TODO @mths for several IdentificationRuns this must be something else, otherwise for two of the same db just one will be created
        std::map<String, String>::iterator dbit = sdb_ids.find(sdb_file);
        if (dbit == sdb_ids.end())
        {
          sdb_id = "SDB_"+ String(UniqueIdGenerator::getUniqueId());

          search_database += String("\t\t<SearchDatabase ");
          search_database += String("location=\"") + sdb_file + "\" ";
          if (!String(it->getSearchParameters().db_version).empty())
          {
            search_database += String("version=\"") + String(it->getSearchParameters().db_version) + "\" ";
          }
          search_database += String("id=\"") + sdb_id + String("\" > \n\t\t\t<FileFormat> \n ");
          //TODO Searchdb file format type cvParam handling
          search_database += String(4, '\t') + cv_.getTermByName("FASTA format").toXMLString(cv_ns);
          search_database += String("\n\t\t\t</FileFormat>\n\t\t\t<DatabaseName>\n\t\t\t\t<userParam name=\"") + sdb_file + String("\"/>\n\t\t\t</DatabaseName>\n");
          search_database += "\t\t</SearchDatabase> \n";

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
      std::map<String, std::vector<String> > pep_evis; //maps the sequence to the corresponding evidence elements for the next scope
      for (std::vector<PeptideIdentification>::const_iterator it = cpep_id_->begin(); it != cpep_id_->end(); ++it)
      {
        String emz(it->getMZ());
        String ert(it->getRT());
        String sid = it->getMetaValue("spectrum_reference");
        if (sid.empty())
        {
          sid = String(it->getMetaValue("spectrum_id"));
          if (sid.empty())
          {
              if (it->getMZ() != it->getMZ())
            {
              emz = "nan";
              LOG_WARN << "Found no spectrum reference and no mz position of identified spectrum! You are probabliy converting from an old format with insufficient data provision. Setting 'nan' - downstream applications might fail unless you set the references right." << std::endl;
            }
            if (it->getRT() != it->getRT())
            {
              ert = "nan";
              LOG_WARN << "Found no spectrum reference and no RT position of identified spectrum! You are probabliy converting from an old format with insufficient data provision. Setting 'nan' - downstream applications might fail unless you set the references right." << std::endl;
            }
            sid = String("MZ:") + emz + String("@RT:") + ert;
          }
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
          LOG_WARN << "Falling back to referencing first spectrum file given because file or identifier could not be mapped." << std::endl;
        }

        sidres += String("\t\t\t<SpectrumIdentificationResult spectraData_ref=\"")
        //multi identification runs lookup from file_origin here
                + sdr + String("\" spectrumID=\"")
                + sid + String("\" id=\"") + sir + String("\"> \n");
        //map.begin access ok here because make sure at least one "UNKOWN" element is in the sdats_ids map

        for (std::vector<PeptideHit>::const_iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
        {
          String pepid =  "PEP_" + String(UniqueIdGenerator::getUniqueId());
          String pepi = jt->getSequence().toString();
          std::map<String, String>::iterator pit = pep_ids.find(pepi);
          if (pit == pep_ids.end())
          {
            String p;
            //~ TODO simplify mod cv param write
            p += String("\t<Peptide id=\"") + pepid + String("\"> \n\t\t<PeptideSequence>") + jt->getSequence().toUnmodifiedString() + String("</PeptideSequence> \n");
            if (jt->getSequence().isModified())
            {
              ModificationsDB* mod_db = ModificationsDB::getInstance();
              if (!jt->getSequence().getNTerminalModification().empty())
              {
                p += "\t\t<Modification location=\"0\"> \n";
                String mod_str = jt->getSequence().getNTerminalModification();
                std::set<const ResidueModification*> mods;
                mod_db->searchTerminalModifications(mods, mod_str, ResidueModification::N_TERM);
                if (!mods.empty())
                {
                  String acc = (*mods.begin())->getUniModAccession();
                  p += "\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':');
                  p += "\" name=\"" +  mod_str;
                  p += "\" cvRef=\"UNIMOD\"/>";
                }
                else // TODO @mths file issue: as this appears to yield hodgepodge 'id's (sometimes e.g. Gln->pyro-Glu other times UNIMOD accessions) - issue is probably in some idXML writing code or xtandem xml consuming code
                {
                  p += "\t\t\t<cvParam accession=\"NA\" name=\"" +  mod_str + "\" cvRef=\"UNIMOD\"/>";
                }
                p += "\n\t\t</Modification> \n"; // "UNIMOD:" prefix??
              }
              if (!jt->getSequence().getCTerminalModification().empty())
              {
                p += "\t\t<Modification location=\"";
                p += String(jt->getSequence().size());
                p += "\"> \n";
                String mod_str = jt->getSequence().getCTerminalModification();
                std::set<const ResidueModification*> mods;
                mod_db->searchTerminalModifications(mods, mod_str, ResidueModification::C_TERM);
                if (!mods.empty())
                {
                  String acc = (*mods.begin())->getUniModAccession();
                  p += "\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':');
                  p += "\" name=\"" +  mod_str;
                  p += "\" cvRef=\"UNIMOD\"/>";
                }
                else // TODO @mths file issue: as this appears to yield hodgepodge 'id's (sometimes e.g. Gln->pyro-Glu other times UNIMOD accessions) - issue is probably in some idXML writing code or xtandem xml consuming code
                {
                  p += "\t\t\t<cvParam accession=\"NA\" name=\"" +  mod_str + "\" cvRef=\"UNIMOD\"/>";
                }

                p += jt->getSequence().getCTerminalModification(); // "UNIMOD:" prefix??
                p += "\n\t\t</Modification> \n";
              }
              for (Size i = 0; i < jt->getSequence().size(); ++i)
              {
                String mod_str =  jt->getSequence()[i].getModification(); // "UNIMOD:" prefix??
                if (!mod_str.empty())
                {
                  std::set<const ResidueModification*> mods;
                  mod_db->searchModifications(mods, jt->getSequence()[i].getOneLetterCode(), mod_str, ResidueModification::ANYWHERE);
                  if (!mods.empty())
                  {
                    //~ p += jt->getSequence()[i].getModification() + "\t" +  jt->getSequence()[i].getOneLetterCode()  + "\t" +  x +   "\n" ;
                    p += "\t\t<Modification location=\"" + String(i + 1);
                    p += "\" residues=\"" + jt->getSequence()[i].getOneLetterCode();
                    String acc = (*mods.begin())->getUniModAccession();
                    p += "\"> \n\t\t\t<cvParam accession=\"UNIMOD:" + acc.suffix(':'); //TODO @all: do not exclusively use unimod ...
                    p += "\" name=\"" +  mod_str;
                    p += "\" cvRef=\"UNIMOD\"/>";
                    p += "\n\t\t</Modification> \n";
                  }
                }
                /* <psi-pi:SubstitutionModification originalResidue="A" replacementResidue="A"/> */
              }
            }
            p += "\t</Peptide> \n ";
            sen_set.insert(p);
            pep_ids.insert(std::make_pair(pepi, pepid));
          }
          else
          {
            pepid = pit->second;
          }

          std::vector<String> pevid_ids;
          if (pit == pep_ids.end())
          {        
            std::vector<PeptideEvidence> peptide_evidences = jt->getPeptideEvidences();
            // TODO idXML allows peptide hits without protein references! Fails in that case - run PeptideIndexer first
            for (std::vector<PeptideEvidence>::const_iterator pe = peptide_evidences.begin(); pe != peptide_evidences.end(); ++pe)
            {
              String pevid =  "PEV_" + String(UniqueIdGenerator::getUniqueId());
              String dBSequence_ref = String(sen_ids.find(pe->getProteinAccession())->second);
              String idec;
              if (jt->metaValueExists("target_decoy"))
              {
                idec = String(boost::lexical_cast<std::string>((String(jt->getMetaValue("target_decoy"))).hasSubstring("decoy")));
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
                e += " start=\"" + String(pe->getStart()) + "\"";
              }
              else if (jt->metaValueExists("start"))
              {
                e += " start=\"" + String(jt->getMetaValue("start")) + "\"";
              }
              else
              {
                LOG_WARN << "Found no start position of peptide hit in protein sequence." << std::endl;
              }
              if (pe->getEnd() != PeptideEvidence::UNKNOWN_POSITION)
              {
                e += " end=\"" + String(pe->getEnd()) + "\"";
              }
              else if (jt->metaValueExists("end"))
              {
                e += " end=\"" + String(jt->getMetaValue("end")) + "\"";
              }
              else
              {
                LOG_WARN << "Found no end position of peptide hit in protein sequence." << std::endl;
              }
              if (!idec.empty())
              {
                e += " isDecoy=\"" + String(idec)+ "\"";
              }
              e += "/> \n";
              sen_set.insert(e);
              pevid_ids.push_back(pevid);
            }
            pep_evis.insert(make_pair(pepi, pevid_ids));
          }
          else
          {
            pevid_ids =  pep_evis[pepi];
          }

          String cmz((jt->getSequence().getMonoWeight() +  jt->getCharge() * Constants::PROTON_MASS_U) / jt->getCharge()); //calculatedMassToCharge
          String r(jt->getRank()); //rank
          String sc(jt->getScore());
          if (sc.empty())
          {
            sc = "NA";
            LOG_WARN << "No score assigned to this PSM: " /*<< jt->getSequence().toString()*/ << std::endl;
          }
          String c(jt->getCharge()); //charge

          String pte;
          if (jt->metaValueExists("pass_threshold"))
          {
            pte = boost::lexical_cast<std::string>(jt->getMetaValue("pass_threshold"));
          }
          else if (pp_identifier_2_thresh.find(it->getIdentifier())!= pp_identifier_2_thresh.end() && pp_identifier_2_thresh.find(it->getIdentifier())->second != 0.0)
          {
            double th = pp_identifier_2_thresh.find(it->getIdentifier())->second;
            //threshold was 'set' in proteinIdentification (!= default value of member, now check pass
            pte = boost::lexical_cast<std::string>(it->isHigherScoreBetter() ? jt->getScore() > th : jt->getScore() < th); //passThreshold-eval
          }
          else
          {
//            if (pp_identifier_2_thresh.find(it->getIdentifier())== pp_identifier_2_thresh.end())
//            {
//                warning("Fuuuuuu");
//            }
            pte = true;
          }

          //write SpectrumIdentificationItem elements
          String sii = "SII_" + String(UniqueIdGenerator::getUniqueId());
          sidres += String("\t\t\t\t<SpectrumIdentificationItem passThreshold=\"")
                  + pte + String("\" rank=\"") + r + String("\" peptide_ref=\"")
                  + pepid + String("\" calculatedMassToCharge=\"") + cmz
                  + String("\" experimentalMassToCharge=\"") + emz
                  + String("\" chargeState=\"") + c +  String("\" id=\"")
                  + sii + String("\"> \n");

          if (pevid_ids.empty())
          {
            LOG_WARN << "PSM without peptide evidence registered in the given search database found. This will cause an invalid mzIdentML file (which OpenMS can still consume)." << std::endl;
          }
          for (std::vector<String>::const_iterator pevref = pevid_ids.begin(); pevref != pevid_ids.end(); ++pevref)
          {
            sidres += "\t\t\t\t\t<PeptideEvidenceRef peptideEvidence_ref=\"" +  String(*pevref) + "\"/> \n";
          }

          std::set<String> peptide_result_details;
          cv_.getAllChildTerms(peptide_result_details, "MS:1001143"); // search engine specific score for PSMs
          MetaInfoInterface copy_jt = *jt;
          String st(it->getScoreType()); //scoretype

          if (cv_.hasTermWithName(st) && peptide_result_details.find(cv_.getTermByName(st).id) != peptide_result_details.end())
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName(st).toXMLString(cv_ns, sc);
            copy_jt.removeMetaValue(cv_.getTermByName(st).id);
          }
          else if (cv_.exists(st) && peptide_result_details.find(st) != peptide_result_details.end())
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTerm(st).toXMLString(cv_ns, sc);
            copy_jt.removeMetaValue(cv_.getTerm(st).id);
          }
          else if (st == "q-value" || st == "FDR")
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName("PSM-level q-value").toXMLString(cv_ns, sc);
            copy_jt.removeMetaValue(cv_.getTermByName("PSM-level q-value").id);
          }
          else if (st == "Posterior Error Probability")
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName("percolator:PEP").toXMLString(cv_ns, sc); // 'percolaror' was not a typo in the code but in the cv.
            copy_jt.removeMetaValue(cv_.getTermByName("percolator:PEP").id);
          }
          else if (st == "OMSSA")
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName("OMSSA:evalue").toXMLString(cv_ns, sc);
            copy_jt.removeMetaValue(cv_.getTermByName("OMSSA:evalue").id);
          }
          else if (st == "Mascot")
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName("Mascot:score").toXMLString(cv_ns, sc);
            copy_jt.removeMetaValue(cv_.getTermByName("Mascot:score").id);
          }
          else if (st == "XTandem")
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName("X\\!Tandem:hyperscore").toXMLString(cv_ns, sc);
            copy_jt.removeMetaValue(cv_.getTermByName("X\\!Tandem:hyperscore").id);
          }
          else if (st == "SEQUEST")
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName("Sequest:xcorr").toXMLString(cv_ns, sc);
            copy_jt.removeMetaValue(cv_.getTermByName("Sequest:xcorr").id);
          }
          else if (st == "MS-GF+")
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName("MS-GF:RawScore").toXMLString(cv_ns, sc);
            copy_jt.removeMetaValue(cv_.getTermByName("MS-GF:RawScore").id);
          }
          else
          {
            sidres +=  "\t\t\t\t\t" + cv_.getTermByName("search engine specific score for PSMs").toXMLString(cv_ns, sc);
            LOG_WARN << "Converting unknown score type to search engine specific score from PSI controlled vocabulary." << std::endl;
          }
          sidres += "\n";

          copy_jt.removeMetaValue("calcMZ");
          copy_jt.removeMetaValue("target_decoy");
          writeMetaInfos_(sidres, copy_jt, 5);

          //~ sidres += "<cvParam accession=\"MS:1000796\" cvRef=\"PSI-MS\" value=\"55.835.842.3.dta\" name=\"spectrum title\"/>";
          sidres += "\t\t\t\t</SpectrumIdentificationItem>\n";
        }
        if (!ert.empty() && ert != "nan" && ert != "NaN")
        {
          sidres +=  "\t\t\t\t" + cv_.getTermByName("retention time").toXMLString(cv_ns, ert) + "\n";
        }
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
          LOG_ERROR << "encountered a PeptideIdentification which is not linked to any ProteinIdentification" << std::endl;
        }
      }

      //--------------------------------------------------------------------------------------------
      // XML header
      //--------------------------------------------------------------------------------------------
      os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
         << "<MzIdentML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
         << "\txsi:schemaLocation=\"http://psidev.info/psi/pi/mzIdentML/1.1 http://psi-pi.googlecode.com/svn/trunk/schema/mzIdentML1.1.0.xsd\"\n"
         << "\txmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\"\n"
         << "\tversion=\"1.1.0\"\n"
         << "\tid=\"OpenMS_" << String(UniqueIdGenerator::getUniqueId()) << "\"\n"
         << "\tcreationDate=\"" << DateTime::now().getDate() << "T" << DateTime::now().getTime() << "\">\n";

      //--------------------------------------------------------------------------------------------
      // CV list
      //--------------------------------------------------------------------------------------------
      os << "<cvList> \n \t<cv id=\"PSI-MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Vocabularies\"  uri=\"http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"3.15.0\"></cv> \n \t<cv id=\"UNIMOD\" fullName=\"UNIMOD\"        uri=\"http://www.unimod.org/obo/unimod.obo\"></cv> \n \t<cv id=\"UO\"     fullName=\"UNIT-ONTOLOGY\" uri=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\"></cv>\n</cvList>\n";

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
        os << "\t<AnalysisSoftware version=\"OpenMS TOPP v"<< VersionInfo::getVersion() <<"\" name=\"TOPP software\" id=\"" << String("SOF_") << String(UniqueIdGenerator::getUniqueId()) << "\"> \n"
           << "\t\t<SoftwareName> \n\t\t\t" << cv_.getTermByName("TOPP software").toXMLString(cv_ns) << " \n\t\t</SoftwareName> \n\t</AnalysisSoftware> \n";
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
        os << "\t\t<SpectrumIdentificationList id=\"" << sil_it->first << String("\"> \n");
        os << "\t\t\t<FragmentationTable>\n"
           << "\t\t\t\t<Measure id=\"Measure_" << sil_it->first << "\">\n"
              // TODO as soon as fragmentation table is reflectable by our internal structures, this has to be mapped separately from spectrumidentificationlist
           << "\t\t\t\t\t<cvParam accession=\"MS:1001225\" cvRef=\"PSI-MS\" unitCvRef=\"PSI-MS\" unitName=\"m/z\" unitAccession=\"MS:1000040\" name=\"product ion m/z\"/>\n"
           << "\t\t\t\t</Measure>\n"
           << "\t\t\t</FragmentationTable>\n";
        os << sil_it->second;
        os << "\t\t</SpectrumIdentificationList>\n";
      }
      os << "\t</AnalysisData>\n</DataCollection>\n";

      //--------------------------------------------------------------------------------------------
      // close XML header
      //--------------------------------------------------------------------------------------------
      os << "</MzIdentML>";

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

          DataValue d = meta.getMetaValue(keys[i]);
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
          s += "\" value=\"" + (String)(d) + "\"/>" + "\n";
        }
      }
    }

    void MzIdentMLHandler::writeEnzyme_(String& s, Enzyme enzy, UInt miss, UInt indent) const
    {
      String cv_ns = cv_.name();
      s += String(indent, '\t') + "<Enzymes independent=\"false\">" + "\n";
      s += String(indent, '\t') + "\t" + "<Enzyme missedCleavages=\"" + String(miss) + "\" id=\"" + String("ENZ_") + String(UniqueIdGenerator::getUniqueId()) + "\">" + "\n";
      s += String(indent, '\t') + "\t\t" + "<EnzymeName>" + "\n";
      String enzymename = enzy.getName();
      if (cv_.hasTermWithName(enzymename))
      {
        s += String(indent, '\t') + "\t\t\t" + cv_.getTermByName(enzymename).toXMLString(cv_ns) + "\n";
      }
      else if (enzymename == "no cleavage")
      {
        s += String(indent, '\t') + "\t\t\t" + cv_.getTermByName("NoEnzyme").toXMLString(cv_ns) + "\n";
      }
      else
      {
        s += String(indent, '\t') + "\t\t\t" + cv_.getTermByName("cleavage agent details").toXMLString(cv_ns) + "\n";
      }
      s += String(indent, '\t') + "\t\t" + "</EnzymeName>" + "\n";
      s += String(indent, '\t') + '\t' + "</Enzyme>" + "\n";
      s += String(indent, '\t') + "</Enzymes>" + "\n";
    }

    void MzIdentMLHandler::writeModParam_(String& s, const std::vector<String>& fixed, const std::vector<String>& variable, UInt indent) const
    {
      String cv_ns = unimod_.name();
      s += String(indent, '\t') + "<ModificationParams>" + "\n";
      for (std::vector<String>::const_iterator it = fixed.begin(); it != fixed.end(); ++it)
      {
        std::set<const ResidueModification*> mods;
        ModificationsDB::getInstance()->searchModifications(mods, *it, ResidueModification::ANYWHERE);
        if (!mods.empty())
        {
          for (std::set<const ResidueModification*>::const_iterator mt = mods.begin(); mt != mods.end(); ++mt)
          {
            s += String(indent, '\t') + '\t' + "<SearchModification fixedMod=\"true\" massDelta=\"" + String((*mt)->getMonoMass()) + "\" residues=\"" + String((*mt)->getOrigin()) + "\">" + "\n";
            String ac = (*mt)->getUniModAccession();
            if (ac.hasPrefix("UniMod:"))
              ac = "UNIMOD:" + ac.suffix(':');
            s += String(indent, '\t') + "\t\t" + unimod_.getTerm(ac).toXMLString(cv_ns) + "\n";
            s += String(indent, '\t') + '\t' + "</SearchModification>" + "\n";
          }
        }
        else
        {
          LOG_WARN << "Registered fixed modification not writable, unknown or unable to convert to cv parameter." << std::endl;
        }
      }
      for (std::vector<String>::const_iterator it = variable.begin(); it != variable.end(); ++it)
      {
        std::set<const ResidueModification*> mods;
        ModificationsDB::getInstance()->searchModifications(mods, *it, ResidueModification::ANYWHERE);
        if (!mods.empty())
        {
          for (std::set<const ResidueModification*>::const_iterator mt = mods.begin(); mt != mods.end(); ++mt)
          {
            s += String(indent, '\t') + '\t' + "<SearchModification fixedMod=\"false\" massDelta=\"" + String((*mt)->getMonoMass()) + "\" residues=\"" + String((*mt)->getOrigin()) + "\">" + "\n";
            String ac = (*mt)->getUniModAccession();
            if (ac.hasPrefix("UniMod:"))
              ac = "UNIMOD:" + ac.suffix(':');
            s += String(indent, '\t') + "\t\t" + unimod_.getTerm(ac).toXMLString(cv_ns) + "\n";
            s += String(indent, '\t') + '\t' + "</SearchModification>" + "\n";
          }
        }
        else
        {
          LOG_WARN << "Registered variable modification not writable, unknown or unable to convert to cv parameter." << std::endl;
        }
      }
      s += String(indent, '\t') + "</ModificationParams>" + "\n";
    }
  } //namespace Internal
} // namespace OpenMS
