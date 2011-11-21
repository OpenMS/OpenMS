// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

#include <OpenMS/KERNEL/StandardTypes.h>

#include <set>

#include <boost/lexical_cast.hpp>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{

  MzIdentMLHandler::MzIdentMLHandler(const Identification& id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
    	//~ ms_exp_(0),
			id_(0),
			cid_(&id)
  {
  	cv_.loadFromOBO("PSI-MS",File::find("/CV/psi-ms.obo"));
  	unimod_.loadFromOBO("PSI-MS",File::find("/CV/unimod.obo"));
  }

  MzIdentMLHandler::MzIdentMLHandler(Identification& id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			//~ ms_exp_(0),
			id_(&id),
			cid_(0)
  {
  	cv_.loadFromOBO("PSI-MS",File::find("/CV/psi-ms.obo"));
		unimod_.loadFromOBO("PSI-MS",File::find("/CV/unimod.obo"));
  }
	
	MzIdentMLHandler::MzIdentMLHandler(const std::vector<ProteinIdentification>& pro_id, const std::vector<PeptideIdentification>& pep_id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			//~ ms_exp_(0),
			pro_id_(0),
			pep_id_(0),
			cpro_id_(&pro_id),
			cpep_id_(&pep_id)
  {
  	cv_.loadFromOBO("PSI-MS",File::find("/CV/psi-ms.obo"));
		unimod_.loadFromOBO("PSI-MS",File::find("/CV/unimod.obo"));
  }
	
	MzIdentMLHandler::MzIdentMLHandler(std::vector<ProteinIdentification>& pro_id, std::vector<PeptideIdentification>& pep_id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			//~ ms_exp_(0),
			pro_id_(&pro_id),
			pep_id_(&pep_id),
			cpro_id_(0),
			cpep_id_(0)
  {
  	cv_.loadFromOBO("PSI-MS",File::find("/CV/psi-ms.obo"));
		unimod_.loadFromOBO("PSI-MS",File::find("/CV/unimod.obo"));
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
		if (to_ignore.size() == 0)
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
			 parent_tag = *(open_tags_.end()-2);
		}
    String parent_parent_tag;
    if (open_tags_.size() > 2)
		{
			parent_parent_tag = *(open_tags_.end()-3);
		}

		static const XMLCh* s_value = xercesc::XMLString::transcode("value");
    static const XMLCh* s_unit_accession = xercesc::XMLString::transcode("unitAccession");
    static const XMLCh* s_cv_ref = xercesc::XMLString::transcode("cvRef");
    //~ static const XMLCh* s_name = xercesc::XMLString::transcode("name");
    static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");


		if (tag_ == "cvParam")
		{
      String value, unit_accession, cv_ref;
      optionalAttributeAsString_(value, attributes, s_value);
      optionalAttributeAsString_(unit_accession, attributes, s_unit_accession);
      optionalAttributeAsString_(cv_ref, attributes, s_cv_ref);
      handleCVParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_accession), /* attributeAsString_(attributes, s_name), value, */ attributes, cv_ref/*,  unit_accession */);
			return;
		}

		if (tag_ == "MzIdentML")
		{
			// TODO handle version
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
			DoubleReal double_value(0);
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
			// TODO write customizations to Sofware
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
			actual_peptide_ = pep;
			return;
		}

    //error(LOAD, "MzIdentMLHandler::characters: Unkown character section found: '" + tag_ + "', ignoring.");
	}

	void MzIdentMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		static set<String> to_ignore;
		if (to_ignore.size() == 0)
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
		Residue::ResidueType res_type_ = Residue::Full;
		String cv_ns = cv_.name();
		String datacollection_element, analysissoftwarelist_element, analysisprotocolcollection_element, analysiscollection_element ;
		String inputs_element, analysisdata_element;
		std::set<String> sdb_set, sen_set, sof_set, sip_set, spd_set;     
		std::map<String,UInt64> sdb_ids, sen_ids, sof_ids, sip_ids, spd_ids, pep_ids;
		std::map<String,String> pie_ids;
		std::vector<String> peps,pepevis,sidlist; 
		//TODO MS:1001035 (date / time search performed) for sidlist

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
		inputs_element += String("<Inputs>");
		String spectra_data;
		//~ for (std::vector<PeptideIdentification>::const_iterator it = cpep_id_->begin(); it != cpep_id_->end(); ++it)
		//~ {			
			//~ for (std::vector<PeptideHit>::const_iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
			//~ {	
				//TODO get spectra_data location not build in idxml or internal structures yet
			//~ }
		//~ }
		if (spd_set.empty())
		{
				UInt64 spdid;
				spdid = UniqueIdGenerator::getUniqueId();
				spd_set.insert("UNKNOWN");
				spd_ids.insert(std::pair<String,UInt64>("UNKNOWN",spdid));

				spectra_data += String("<SpectraData location=\"") + String("UNKNOWN") + String("\" id=\"") + String(spdid) + String("\">");
				spectra_data += String("<FileFormat> \n ");
				spectra_data +=  cv_.getTermByName("mzML file").toXMLString(cv_ns); 
				spectra_data += String("\n </FileFormat>\n <SpectrumIDFormat> \n ");
				spectra_data +=  cv_.getTermByName("multiple peak list nativeID format").toXMLString(cv_ns);
				spectra_data += String("\n </SpectrumIDFormat> \n </SpectraData>");
		}
		
		/*
		1st: iterate over proteinidentification vector
		*/
		for (std::vector<ProteinIdentification>::const_iterator it = cpro_id_->begin(); it != cpro_id_->end(); ++it)
		{
			UInt64 dbid;
			std::set<String>::iterator dbit = sdb_set.find(String(it->getSearchParameters().db));
			if (dbit==sdb_set.end())
			{
				dbid = UniqueIdGenerator::getUniqueId();
				String dbst(it->getSearchParameters().db);
				
				inputs_element += String("<SearchDatabase "); 
				inputs_element += String("location=\"") + dbst + "\" ";
				//TODO get version db += String("version=\"") + String(it->getSearchParameters().version) + "\" ";
				inputs_element += String("id=\"") + String(dbid) + String("\" > \n <FileFormat> \n ");
				//TODO Searchdb file format type cvParam handling
				inputs_element += cv_.getTermByName("FASTA format").toXMLString(cv_ns);
				inputs_element += String(" </FileFormat> \n <DatabaseName>\n <userParam name=\"") + dbst + String("\"/>\n </DatabaseName>\n");
				inputs_element += "</SearchDatabase> \n";
				
				sdb_ids.insert(std::pair<String,UInt64>(dbst,dbid));
			}
			else
			{
				dbid = sdb_ids.find(*dbit)->second;
			}
			
			//~ get a map from identifier to match OpenMS Protein/PeptideIdentification match string;
			pie_ids.insert(std::pair<String,String>(it->getIdentifier(),it->getSearchEngine())); 	
			
			//~ collect analysissoftware in this loop - does not go into inputelement
			UInt64 swid;
			String swcn = String(it->getSearchEngine());
			std::map<String,UInt64>::iterator soit = sof_ids.find(swcn);
			String osecv;
			if (swcn == "OMSSA")
			{
				osecv = "OMSSA";
			} else if (swcn == "Mascot")
			{
				osecv = "MASCOT";
			} else if (swcn == "XTandem")
			{
				osecv = "xtandem";
			} else if (swcn == "SEQUEST")
			{
				osecv = "Sequest";
			} else
			{
				osecv = "analysis software";
			}
			
			if (soit==sof_ids.end())
			{
				swid = UniqueIdGenerator::getUniqueId();
				//~ TODO consider not only searchengine but also version!
				String sost = String("<AnalysisSoftware version=\"") + String(it->getSearchEngineVersion()) + String("\" name=\"") + swcn +  String("\" id=\"") + String(swid) + String("\"> \n")+ String("<SoftwareName> \n "); 				
				sost += cv_.getTermByName(osecv).toXMLString(cv_ns);  
				sost += String("\n </SoftwareName> \n </AnalysisSoftware> \n");
				sof_set.insert(sost);
				sof_ids.insert(std::pair<String,UInt64>(swcn,swid));
			}
			else
			{
				swid = soit->second;
			}

			//~ collect SpectrumIdentificationProtocol for analysisprotocol in this loop - does not go into inputelement
			std::map<String,UInt64>::iterator spit = sip_ids.find(swcn);
			if (spit==sip_ids.end())
			{
				UInt64 spid = UniqueIdGenerator::getUniqueId();
				String sip = String("<SpectrumIdentificationProtocol id=\"") + String(spid) + String("\" analysisSoftware_ref=\"")  + String(swid) + String("\"> \n <SearchType>\n");
				sip += cv_.getTermByName("ms-ms search").toXMLString(cv_ns);
				sip += String(" \n </SearchType>\n<Threshold>\n");
				sip += cv_.getTermByName("no threshold").toXMLString(cv_ns);
				sip += String("\n</Threshold>\n</SpectrumIdentificationProtocol>\n");
				sip_set.insert(sip);
				sip_ids.insert(std::pair<String,UInt64>(swcn,spid));
			}
			
			for (std::vector<ProteinHit>::const_iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
			{
					UInt64 enid; 
					std::map<String,UInt64>::iterator enit = sen_ids.find(String(jt->getAccession()));
					if (enit==sen_ids.end())
					{
						String entry;
						enid = UniqueIdGenerator::getUniqueId(); //TODO IDs from metadata or where its stored at read in;
						String enst(jt->getAccession());
						
						entry += "<DBSequence accession=\"" + enst + "\" ";
						entry += "searchDatabase_ref=\"" + String(dbid) + "\" ";
						entry += "length=\"" + String(jt->getSequence().length()) + "\" ";
						entry += String("id=\"") + String(enid) + String("\">\n");
						String s = String(jt->getSequence());
						if (!s.empty())
						{
							entry += "<Seq>" + s + "</Seq>\n";
						}
						entry += cv_.getTermByName("protein description").toXMLString(cv_ns, enst);
						entry += "</DBSequence>\n";
						
						sen_ids.insert(std::pair<String,UInt64>(enst,enid));
						sen_set.insert(entry);
						
					}
					else
					{
						enid = enit->second;
					}
			}
			
		}
		inputs_element += spectra_data;
		inputs_element += "</Inputs>\n";
	
		/*
		2nd: iterate over peptideidentification vector
		*/
		for (std::vector<PeptideIdentification>::const_iterator it = cpep_id_->begin(); it != cpep_id_->end(); ++it)
		{			
			String pro_pep_matchstring = it->getIdentifier();//~ TODO getIdentifier() lookup in proteinidentification get search db etc

			for (std::vector<PeptideHit>::const_iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
			{	
				String pepi = jt->getSequence().toString();
				UInt64 pepid =  UniqueIdGenerator::getUniqueId();
				
				std::map<String,UInt64>::iterator pit = pep_ids.find(pepi);
				if (pit==pep_ids.end())
				{
					String p;
					//~ TODO simplify mod cv param write
					p += String("<Peptide id=\"") + String(pepid) + String("\"> \n <PeptideSequence>") +jt->getSequence().toUnmodifiedString() + String("</PeptideSequence> \n");
					if(jt->getSequence().isModified())
					{
						ModificationsDB* mod_db = ModificationsDB::getInstance();
						if(!jt->getSequence().getNTerminalModification().empty())
						{ 
							p += "<Modification location=\"0\"> \n <cvParam accession=\"";
							p += jt->getSequence().getNTerminalModification(); // "UNIMOD:" prefix??
							p += "\" cvRef=\"UNIMOD\"/> \n </Modification> \n";
						}
						if(!jt->getSequence().getCTerminalModification().empty())
						{
							p += "<Modification location=\"";
							p += String(jt->getSequence().size());
							p += "\"> \n <cvParam accession=\"";
							p += jt->getSequence().getCTerminalModification(); // "UNIMOD:" prefix??
							p += "\" cvRef=\"UNIMOD\"/> \n </Modification> \n";
						}
						for (Size i = 0; i < jt->getSequence().size(); ++i)
						{
							String mod_str =  jt->getSequence()[i].getModification(); // "UNIMOD:" prefix??
							if (!mod_str.empty())
							{
								std::set<const ResidueModification*> mods;
								mod_db->searchModifications(mods, jt->getSequence()[i].getOneLetterCode(), mod_str,ResidueModification::ANYWHERE);
								if(!mods.empty()) 
								{								
									//~ p += jt->getSequence()[i].getModification() + "\t" +  jt->getSequence()[i].getOneLetterCode()  + "\t" +  x +   "\n" ;
									p += "<Modification location=\"" + String(i+1);
									p += "\" residues=\"" + jt->getSequence()[i].getOneLetterCode();
									String acc = (*mods.begin())->getUniModAccession();
									p += "\"> \n <cvParam accession=\"UNIMOD:"+ acc.suffix(':'); //TODO repair ResidueModification which gives UniMod ... anyways do not exclusively use unimod ... 
									p += "\" name=\"" +  mod_str;
									p += "\" cvRef=\"UNIMOD\"/>";
									p += "\n </Modification> \n";
								}            
							}
							/* <psi-pi:SubstitutionModification originalResidue="A" replacementResidue="A"/> */
						}
					}
					p += "</Peptide> \n ";
					sen_set.insert(p);
				}
				else
				{
					pepid = pit->second;
				}
				
				std::vector<String> accs = jt->getProteinAccessions(); //TODO idxml allows peptidehits without protein_refs!!! Fails in that case run peptideindexer first
				std::vector<UInt64> pevid_ids;
				for  (std::vector<String>::const_iterator at = accs.begin(); at != accs.end(); ++at)
				{
					UInt64 pevid =  UniqueIdGenerator::getUniqueId();
					String dBSequence_ref = String(sen_ids.find(*at)->second);

					String e;
					String idec(boost::lexical_cast<std::string>((String(jt->getMetaValue("target_decoy"))).hasSubstring("decoy")));
					e += "<PeptideEvidence id=\"" + String(pevid) + "\" peptide_ref=\"" + String(pepid) + "\" dBSequence_ref=\"" + dBSequence_ref;
					
					//~ TODO no '*' allowed!!
					String po = String(jt->getAAAfter());
					if (! po.empty() && po != " " && po != "*")
					{
						e += "\" post=\"" + po;
					}
					String pr = String(jt->getAABefore());
          if (! pr.empty() && pr != " " && pr != "*")
					{
						e += "\" pre=\"" + pr;
					}
					e += "\" isDecoy=\"" + String(idec) + "\"> \n";
										e += "</PeptideEvidence>\n";
					sen_set.insert(e);
					pevid_ids.push_back(pevid);
				}
				
				String cmz(jt->getSequence().getMonoWeight(res_type_, jt->getCharge())); //calculatedMassToCharge
				String emz(String(it->getMetaValue("MZ"))); //precursor MassToCharge
				String ert(String(it->getMetaValue("RT"))); //precursor MassToCharge
				String r(jt->getRank());//rank
				String sc(jt->getScore()); //score TODO what if there is no score?
				String st(it->getScoreType()); //scoretype
				String c(jt->getCharge()); //charge
				String pte(boost::lexical_cast<std::string>(it->isHigherScoreBetter()?jt->getScore()>it->getSignificanceThreshold():jt->getScore()<it->getSignificanceThreshold())); //passThreshold-eval
				
				String sidres;
				UInt64 sir =  UniqueIdGenerator::getUniqueId();
				UInt64 sii =  UniqueIdGenerator::getUniqueId();
				//~ TODO get spectra_data when loading idxml if possible - then add here and in spectra_data section
				sidres += String("<SpectrumIdentificationResult spectraData_ref=\"") + String(spd_ids.begin()->second) + String("\" spectrumID=\"") + String("MZ:")+ emz + String("@RT:") + ert + String("\" id=\"") + String(sir) + String("\"> \n"); //map.begin access ok here because make sure at least one "UNKOWN" element is in the spd_ids map
				sidres += String("<SpectrumIdentificationItem passThreshold=\"") + pte + String("\" rank=\"") + r + String("\" peptide_ref=\"") + String(pepid) + String("\" calculatedMassToCharge=\"") + cmz + String("\" experimentalMassToCharge=\"") + emz + String("\" chargeState=\"") + c +  String("\" id=\"") + String(sii) + String("\"> \n");
				for  (std::vector<UInt64>::const_iterator pevref = pevid_ids.begin(); pevref != pevid_ids.end(); ++pevref)
				{
					sidres += "<PeptideEvidenceRef peptideEvidence_ref=\"" +  String(*pevref) + "\"/> \n";
				}

				//~ TODO nicer cvParam handling!
				if (st == "q-value" || st == "FDR")
				{
					sidres +=  cv_.getTermByName("pep:global FDR").toXMLString(cv_ns,sc);
				}else if (st == "Posterior Error Probability")
				{
					sidres +=  cv_.getTermByName("percolaror:PEP").toXMLString(cv_ns,sc); // 'percolaror' is not a typo.
				}else if (pie_ids[pro_pep_matchstring] == "OMSSA")
				{
					sidres +=  cv_.getTermByName("OMSSA:evalue").toXMLString(cv_ns,sc);
				} else if (pie_ids[pro_pep_matchstring] == "Mascot")
				{
					sidres +=  cv_.getTermByName("MASCOT:score").toXMLString(cv_ns,sc);
				} else if (pie_ids[pro_pep_matchstring] == "XTandem")
				{
					sidres +=  cv_.getTermByName("X!Tandem:expect").toXMLString(cv_ns,sc);
				} else if (pie_ids[pro_pep_matchstring] == "SEQUEST")
				{
					sidres +=  cv_.getTermByName("Sequest:xcorr").toXMLString(cv_ns,sc);
				} else
				{
					sidres +=  cv_.getTermByName("search engine specific score for peptides").toXMLString(cv_ns,sc);
				}
				//~ sidres += "<cvParam accession=\"MS:1000796\" cvRef=\"PSI-MS\" value=\"55.835.842.3.dta\" name=\"spectrum title\"/>";
				sidres += "</SpectrumIdentificationItem>\n </SpectrumIdentificationResult>\n";
				sidlist.push_back(sidres);
				
			}
		}
		
		//--------------------------------------------------------------------------------------------
		// XML header
		//--------------------------------------------------------------------------------------------
		os	<< "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
				<< "<MzIdentML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
				<< "xsi:schemaLocation=\"http://psidev.info/psi/pi/mzIdentML/1.1 http://psi-pi.googlecode.com/svn/trunk/schema/mzIdentML1.1.0.xsd\"\n"
				<< "xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\" \n  id=\"\" \n version=\"1.1.0\"\n"
				<< "creationDate=\"2011-11-11T11:11:11\">\n"; 
				//~ TODO dateTime
		//--------------------------------------------------------------------------------------------
		// CV list
		//--------------------------------------------------------------------------------------------
		os << "<cvList> \n <cv id=\"PSI-MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Vocabularies\"  uri=\"http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"3.15.0\"></cv> \n <cv id=\"UNIMOD\" fullName=\"UNIMOD\"        uri=\"http://www.unimod.org/obo/unimod.obo\"></cv> \n <cv id=\"UO\"     fullName=\"UNIT-ONTOLOGY\" uri=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\"></cv>\n</cvList>\n";

		//--------------------------------------------------------------------------------------------
		// AnalysisSoftwareList
		//--------------------------------------------------------------------------------------------
		os << "<AnalysisSoftwareList>\n";
		for  (std::set<String>::const_iterator sof = sof_set.begin(); sof != sof_set.end(); ++sof)
		{
			os << *sof;
		}
		
		std::map<String,UInt64>::iterator soit = sof_ids.find("TOPP software");
		if (soit==sof_ids.end())
		{
			 os << "<AnalysisSoftware version=\"OpenMS TOPP v1.9\" name=\"TOPP software\" id=\"" << String(UniqueIdGenerator::getUniqueId()) << "\"> \n" 
					 << "<SoftwareName> \n " << cv_.getTermByName("TOPP software").toXMLString(cv_ns) << " \n </SoftwareName> \n </AnalysisSoftware> \n";
		}
		os << "</AnalysisSoftwareList>\n";

		//--------------------------------------------------------------------------------------------
		// SequenceCollection
		//--------------------------------------------------------------------------------------------
		os << "<SequenceCollection>\n";
		for  (std::set<String>::const_iterator sen = sen_set.begin(); sen != sen_set.end(); ++sen)
		{
			os << *sen;
		}
		os << "</SequenceCollection>\n";

		//--------------------------------------------------------------------------------------------
		// AnalysisCollection:
		//+SpectrumIdentification + SpectrumIdentification
		//TODO ProteinDetection
		//--------------------------------------------------------------------------------------------
		//~ TODO for every pair of input vector<Protein/PeptideIdentification> create one SpectrumIdentificationList and fill 
		UInt64 silly  = UniqueIdGenerator::getUniqueId();

		os <<	"<AnalysisCollection>";
		for  (std::map<String,UInt64>::const_iterator sip = sip_ids.begin(); sip != sip_ids.end(); ++sip)
		{
			//~ TODO unsure when to create several lists instead of one SpectrumIdentificationList - for now only one list
			//~ for  (std::set<String>::const_iterator sip = sip_set.begin(); sip != sip_set.end(); ++sip)
			//~ {
				UInt64 ss  = UniqueIdGenerator::getUniqueId();
				String entry = String("<SpectrumIdentification id=\"") + String(ss) + String("\" spectrumIdentificationProtocol_ref=\"") + String(sip->second) + String("\" spectrumIdentificationList_ref=\"") + String(silly) + String("\">") + String("<InputSpectra/>\n<SearchDatabaseRef/>\n </SpectrumIdentification>\n");
				os <<	entry;
			//~ }
		}
		os << "</AnalysisCollection>\n";	
		
		//--------------------------------------------------------------------------------------------
		// AnalysisProtocolCollection
		//+SpectrumIdentificationProtocol + SearchType
		//  																									+ Threshold
		//--------------------------------------------------------------------------------------------
		os << "<AnalysisProtocolCollection>\n";
		for  (std::set<String>::const_iterator sip = sip_set.begin(); sip != sip_set.end(); ++sip)
		{
			os << *sip;
		}
		os << "</AnalysisProtocolCollection>\n";

		//--------------------------------------------------------------------------------------------
		// DataCollection
		//+Inputs 
		//+AnalysisData 
		//--------------------------------------------------------------------------------------------
		os << "<DataCollection>\n" << inputs_element;
		//~ TODO for every pair of input vector<Protein/PeptideIdentification> create one SpectrumIdentificationList and fill 
		os << "<AnalysisData>\n<SpectrumIdentificationList id=\"" 
				<< String(silly) 
				<< String("\"> \n");
		for  (std::vector<String>::const_iterator sid = sidlist.begin(); sid != sidlist.end(); ++sid)
		{
			os << *sid;
		}
		os << "</SpectrumIdentificationList>\n</AnalysisData>\n</DataCollection>\n";

		//--------------------------------------------------------------------------------------------
		// close XML header
		//--------------------------------------------------------------------------------------------
		os << "</MzIdentML>";

	}
	} //namespace Internal
} // namespace OpenMS


