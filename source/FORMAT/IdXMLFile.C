

// -*- Mode: C++; tab-widt: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS 
{

	IdXMLFile::IdXMLFile()
		: XMLHandler(""),
			last_meta_(0)
	{
	  	
	}

  void IdXMLFile::load(const String& filename,  vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids)
  	 throw (Exception::FileNotFound, Exception::ParseError)
  {
  	//Filename for error messages in XMLHandler
  	file_ = filename;
  	
  	protein_ids.clear();
  	peptide_ids.clear();
  	
  	prot_ids_ = &protein_ids;
  	pep_ids_ = &peptide_ids;
  	//try to open file
		if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
		
		// initialize parser
		try 
		{
			xercesc::XMLPlatformUtils::Initialize();
		}
		catch (const xercesc::XMLException& toCatch) 
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
	  }
		xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
		parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
		parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);

		parser->setContentHandler(const_cast<IdXMLFile*>(this));
		parser->setErrorHandler(const_cast<IdXMLFile*>(this));
		
		xercesc::LocalFileInputSource source( xercesc::XMLString::transcode(filename.c_str()) );
		try 
    {
    	parser->parse(source);
    	delete(parser);
    }
    catch (const xercesc::XMLException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
    }
    catch (const xercesc::SAXException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
    }
    
    //reset members
    prot_ids_ = 0;
		pep_ids_ = 0;
		last_meta_ = 0;
		parameters_.clear();
		param_ = ProteinIdentification::SearchParameters();
		String id_ = "";
		prot_id_ = ProteinIdentification();
		pep_id_ = PeptideIdentification();
		prot_hit_ = ProteinHit();
		pep_hit_ = PeptideHit();
		proteinid_to_accession_.clear();
  }
  					 
  void IdXMLFile::store(String filename, const vector<ProteinIdentification>& protein_ids, const vector<PeptideIdentification>& peptide_ids) throw (Exception::UnableToCreateFile)
  {
  		//open stream
		std::ofstream os(filename.c_str());
		if (!os)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		//write header
		os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
		os << "<IdXML xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/IdXML_1_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << endl;
		
		//look up different search parameters
		vector<ProteinIdentification::SearchParameters> params;
		for (vector<ProteinIdentification>::const_iterator it = protein_ids.begin(); it!=protein_ids.end(); ++it)
		{
			if (find(params.begin(), params.end(), it->getSearchParameters())==params.end())
			{
				params.push_back(it->getSearchParameters());
			}
		}
		
		//write search parameters
		for(UInt i=0; i!=params.size();++i)
		{
			os << "  <SearchParameters "
				 << "id=\"SP_" << i << "\" "
				 << "db=\"" << params[i].db << "\" "
				 << "db_version=\"" << params[i].db_version << "\" "
				 << "taxonomy=\"" << params[i].taxonomy << "\" ";
			if (params[i].mass_type == ProteinIdentification::MONOISOTOPIC)
			{ 
				os << "mass_type=\"monoisotopic\" ";
			}
			else if (params[i].mass_type == ProteinIdentification::AVERAGE)
			{ 
				os << "mass_type=\"average\" ";
			}
			os << "charges=\"" << params[i].charges << "\" ";
			if (params[i].enzyme == ProteinIdentification::TRYPSIN)
			{ 
				os << "enzyme=\"trypsin\" ";
			}
			if (params[i].enzyme == ProteinIdentification::PEPSIN_A)
			{ 
				os << "enzyme=\"pepsin_a\" ";
			}
			if (params[i].enzyme == ProteinIdentification::PROTEASE_K)
			{ 
				os << "enzyme=\"protease_k\" ";
			}
			if (params[i].enzyme == ProteinIdentification::CHYMOTRYPSIN)
			{ 
				os << "enzyme=\"chymotrypsin\" ";
			}
			else if (params[i].enzyme == ProteinIdentification::NO_ENZYME)
			{ 
				os << "enzyme=\"no_enzyme\" ";
			}
			else if (params[i].enzyme == ProteinIdentification::UNKNOWN_ENZYME)
			{ 
				os << "enzyme=\"unknown_enzyme\" ";
			}
			os << "missed_cleavages=\"" << params[i].missed_cleavages << "\" "
				 << "precursor_peak_tolerance=\"" << params[i].precursor_tolerance << "\" "
				 << "peak_mass_tolerance=\"" << params[i].peak_mass_tolerance << "\" "
				 << ">" << endl;
			
			//modifications
			for (UInt j=0; j!=params[i].fixed_modifications.size(); ++j)
			{
				os << "    <FixedModification name=\"" << params[i].fixed_modifications[j] << "\" />" << endl;
				//Add MetaInfo, when modifications has it (Andreas)
			}
			for (UInt j=0; j!=params[i].variable_modifications.size(); ++j)
			{
				os << "    <VariableModification name=\"" << params[i].variable_modifications[j] << "\" />" << endl;
				//Add MetaInfo, when modifications has it (Andreas)
			}
			
			writeUserParam_(os, params[i], 4);
			
			os << "  </SearchParameters>" << endl;
		}
		
		UInt prot_count = 0;
		map<String,UInt> accession_to_id;
		
		//Identifiers of protein identifications that are already written
		vector<String> done_identifiers;
		
		//write ProteinIdentification Runs
		for (UInt i=0; i<protein_ids.size(); ++i)
		{
			done_identifiers.push_back(protein_ids[i].getIdentifier());
			
			os << "  <IdentificationRun ";
			String time, date;
			protein_ids[i].getDateTime().getDate(date);
			protein_ids[i].getDateTime().getTime(time);
			os << "date=\"" << date << "T" << time << "\" ";
			os << "search_engine=\"" << protein_ids[i].getSearchEngine() << "\" ";
			os << "search_engine_version=\"" << protein_ids[i].getSearchEngineVersion() << "\" ";
			//identifier
			for(UInt j=0; j!=params.size();++j)
			{
				if (params[j]==protein_ids[i].getSearchParameters())
				{
					os << "search_parameters_ref=\"SP_" << j << "\" ";
					break;
				}
			}
			os << ">" << endl;
			os << "    <ProteinIdentification ";
			os << "score_type=\"" << protein_ids[i].getScoreType() << "\" ";
			if (protein_ids[i].isHigherScoreBetter())
			{
				os << "higher_score_better=\"true\" ";
			}
			else
			{
				os << "higher_score_better=\"false\" ";	
			}
			os << "significance_threshold=\"" << protein_ids[i].getSignificanceThreshold() << "\" >" << endl;
			
			//write protein hits
			for(UInt j=0; j<protein_ids[i].getHits().size(); ++j)
			{
				os << "      <ProteinHit ";
				os << "id=\"PH_" << prot_count << "\" ";
				accession_to_id[protein_ids[i].getHits()[j].getAccession()] = prot_count++;
				os << "accession=\"" << protein_ids[i].getHits()[j].getAccession() << "\" ";
				os << "score=\"" << protein_ids[i].getHits()[j].getScore() << "\" ";
				os << "sequence=\"" << protein_ids[i].getHits()[j].getSequence() << "\" >" << endl;
				writeUserParam_(os, protein_ids[i].getHits()[j], 4);
				os << "      </ProteinHit>" << endl;
			}
			
			writeUserParam_(os, protein_ids[i], 3);
			os << "    </ProteinIdentification>" << endl;

			//write PeptideIdentifications
			for (UInt l=0; l<peptide_ids.size(); ++l)
			{
				if (peptide_ids[l].getIdentifier()==protein_ids[i].getIdentifier())
				{
					os << "    <PeptideIdentification ";
					os << "score_type=\"" << peptide_ids[l].getScoreType() << "\" ";
					if (peptide_ids[l].isHigherScoreBetter())
					{
						os << "higher_score_better=\"true\" ";
					}
					else
					{
						os << "higher_score_better=\"false\" ";	
					}
					os << "significance_threshold=\"" << peptide_ids[l].getSignificanceThreshold() << "\" ";
					//mz
					DataValue dv = peptide_ids[l].getMetaValue("MZ");
					if (dv!=DataValue::EMPTY)
					{
						os << "MZ=\"" << dv.toString() << "\" ";
					}
					//rt
					dv = peptide_ids[l].getMetaValue("RT");
					if (dv!=DataValue::EMPTY)
					{
						os << "RT=\"" << dv.toString() << "\" ";
					}
					//spectrum_reference
					dv = peptide_ids[l].getMetaValue("spectrum_reference");
					if (dv!=DataValue::EMPTY)
					{
						os << "spectrum_reference=\"" << dv.toString() << "\" ";
					}
					os << ">" << endl;
					
					//write peptide hits
					for(UInt j=0; j<peptide_ids[l].getHits().size(); ++j)
					{
						os << "      <PeptideHit ";
						os << "score=\"" << peptide_ids[l].getHits()[j].getScore() << "\" ";
						os << "sequence=\"" << peptide_ids[l].getHits()[j].getSequence() << "\" ";
						os << "charge=\"" << peptide_ids[l].getHits()[j].getCharge() << "\" ";
						if (peptide_ids[l].getHits()[j].getAABefore()!=' ')
						{
							os << "aa_before=\"" << peptide_ids[l].getHits()[j].getAABefore() << "\" ";
						}
						if (peptide_ids[l].getHits()[j].getAAAfter()!=' ')
						{
							os << "aa_after=\"" << peptide_ids[l].getHits()[j].getAAAfter() << "\" ";
						}	
						if(peptide_ids[l].getHits()[j].getProteinAccessions().size()!=0)
						{
							String accs = "";
							for (UInt m=0; m<peptide_ids[l].getHits()[j].getProteinAccessions().size(); ++m)
							{
								if (accs!="")
								{
									accs = accs + " ";
								}
								accs = accs + "PH_" + accession_to_id[peptide_ids[l].getHits()[j].getProteinAccessions()[m]];
							}
							os << "protein_refs=\"" << accs << "\" ";
						}
						os << ">" << endl;
						writeUserParam_(os, peptide_ids[l].getHits()[j], 4);
						os << "      </PeptideHit>" << endl;
					}
					
					writeUserParam_(os, peptide_ids[l], 3);
					os << "    </PeptideIdentification>" << endl;
				}
			}

			os << "  </IdentificationRun>" << endl;
		}
		
		for (UInt i=0; i<peptide_ids.size(); ++i)
		{
			if (find(done_identifiers.begin(), done_identifiers.end(), peptide_ids[i].getIdentifier())==done_identifiers.end())
			{
				cerr << "Warning (IdXMLFile): Omitting peptide identification because of missing ProteinIdentification with identifier '" << peptide_ids[i].getIdentifier() << "'!" << endl;
			}
		}
		//write footer
		os << "</IdXML>" << endl;
		
		//close stream
		os.close();
  
  	//reset members
    prot_ids_ = 0;
		pep_ids_ = 0;
		last_meta_ = 0;
		parameters_.clear();
		param_ = ProteinIdentification::SearchParameters();
		String id_ = "";
		prot_id_ = ProteinIdentification();
		pep_id_ = PeptideIdentification();
		prot_hit_ = ProteinHit();
		pep_hit_ = PeptideHit();
		proteinid_to_accession_.clear();
  }

  
	void IdXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{		
		String element = xercesc::XMLString::transcode(qname);
		
		//cout << "Start: " << element << endl;
		
		//START
		if (element == "IdXML")
		{
	
		}
		
		//SEARCH PARAMETERS
		else if (element == "SearchParameters")
		{
			//store id
			id_ =  attributeAsString(attributes,"id");
			
			//reset parameters
			param_ = ProteinIdentification::SearchParameters();
			
			//load parameters
			param_.db = attributeAsString(attributes,"db");
			param_.db_version = attributeAsString(attributes,"db_version");

			optionalAttributeAsString(param_.taxonomy, attributes,"taxonomy");
			param_.charges = attributeAsString(attributes,"charges");
			optionalAttributeAsUInt(param_.missed_cleavages, attributes,"missed_cleavages");
			param_.peak_mass_tolerance = attributeAsDouble(attributes,"peak_mass_tolerance");
			param_.precursor_tolerance = attributeAsDouble(attributes,"precursor_peak_tolerance");
			//mass type
			const XMLCh* mass_type = attributes.getValue(xercesc::XMLString::transcode("mass_type"));
			if (xercesc::XMLString::equals(mass_type,xercesc::XMLString::transcode("monoisotopic")))
			{
				param_.mass_type = ProteinIdentification::MONOISOTOPIC;
			}
			else if (xercesc::XMLString::equals(mass_type,xercesc::XMLString::transcode("average")))
			{
				param_.mass_type = ProteinIdentification::AVERAGE;
			}
			//enzyme
			const XMLCh* enzyme = attributes.getValue(xercesc::XMLString::transcode("enzyme"));
			if (enzyme!=0)
			{
				if (xercesc::XMLString::equals(enzyme,xercesc::XMLString::transcode("trypsin")))
				{
					param_.enzyme = ProteinIdentification::TRYPSIN;
				}
				else if (xercesc::XMLString::equals(enzyme,xercesc::XMLString::transcode("pepsin_a")))
				{
					param_.enzyme = ProteinIdentification::PEPSIN_A;
				}
				else if (xercesc::XMLString::equals(enzyme,xercesc::XMLString::transcode("protease_k")))
				{
					param_.enzyme = ProteinIdentification::PROTEASE_K;
				}
				else if (xercesc::XMLString::equals(enzyme,xercesc::XMLString::transcode("chymotrypsin")))
				{
					param_.enzyme = ProteinIdentification::CHYMOTRYPSIN;
				}			 
				else if (xercesc::XMLString::equals(enzyme,xercesc::XMLString::transcode("no_enzyme")))
				{
					param_.enzyme = ProteinIdentification::NO_ENZYME;
				}
				else if (xercesc::XMLString::equals(enzyme,xercesc::XMLString::transcode("unknown_enzyme")))
				{
					param_.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
				}
			}
			last_meta_ = &param_;	
		}
		else if (element == "FixedModification")
		{
			param_.fixed_modifications.push_back(attributeAsString(attributes,"name"));
			//change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
			last_meta_ = 0;
		}
		else if (element == "VariableModification")
		{
			param_.variable_modifications.push_back(attributeAsString(attributes,"name"));
			//change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
			last_meta_ = 0;
		}
		
		// RUN
		else if (element == "IdentificationRun")
		{
			pep_id_ = PeptideIdentification();
			prot_id_ = ProteinIdentification();

			prot_id_.setSearchEngine(attributeAsString(attributes,"search_engine"));
			prot_id_.setSearchEngineVersion(attributeAsString(attributes,"search_engine_version"));

			//search parameters
			String ref = attributeAsString(attributes,"search_parameters_ref");
			if (parameters_.find(ref)==parameters_.end())
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Invalid search parameters reference '") + ref + "'" );
			}
			prot_id_.setSearchParameters(parameters_[ref]);
			
			//date
			prot_id_.setDateTime(DateTime::fromString(attributeAsString(attributes,"date"),"yyyy-MM-ddThh:mm:ss"));
			
			//set identifier
			prot_id_.setIdentifier(prot_id_.getSearchEngine() + '_' + attributeAsString(attributes,"date"));
		}
		
		//PROTE ProteinIdentification
		else if (element == "ProteinIdentification")
		{
			prot_id_.setScoreType(attributeAsString(attributes,"score_type"));
			
			//optional significance threshold
			DoubleReal tmp=0.0;
			optionalAttributeAsDouble(tmp,attributes,"significance_threshold");
			if (tmp!=0.0)
			{
				prot_id_.setSignificanceThreshold(tmp);
			}
			
			//score orientation
			const XMLCh* higher_score_better = attributes.getValue(xercesc::XMLString::transcode("higher_score_better"));
			if (xercesc::XMLString::equals(higher_score_better,xercesc::XMLString::transcode("true")))
			{
				prot_id_.setHigherScoreBetter(true);	
			}
			else if (xercesc::XMLString::equals(higher_score_better,xercesc::XMLString::transcode("false")))
			{
				prot_id_.setHigherScoreBetter(false);					
			}
			else
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", "Invalid value for 'higher_score_better '");				
			}
			last_meta_ = &prot_id_;
		}
		else if (element == "ProteinHit")
		{
			prot_hit_ = ProteinHit();
			String accession = attributeAsString(attributes,"accession");
			prot_hit_.setAccession(accession);
			prot_hit_.setScore(attributeAsDouble(attributes,"score"));
			
			//sequence
			String tmp="";
			optionalAttributeAsString(tmp,attributes,"sequence");
			prot_hit_.setSequence(tmp);
			
			last_meta_ = &prot_hit_;			
			
			//insert id and accession to map
			proteinid_to_accession_[attributeAsString(attributes,"id")] = accession;
		}
		
		//PEPTIDES
		else if (element == "PeptideIdentification")
		{
			
			//set identifier
			pep_id_.setIdentifier(prot_ids_->back().getIdentifier());
			
			pep_id_.setScoreType(attributeAsString(attributes,"score_type"));
			
			//optional significance threshold
			DoubleReal tmp=0.0;
			optionalAttributeAsDouble(tmp,attributes,"significance_threshold");
			if (tmp!=0.0)
			{
				pep_id_.setSignificanceThreshold(tmp);
			}

			//score orientation
			const XMLCh* higher_score_better = attributes.getValue(xercesc::XMLString::transcode("higher_score_better"));
			if (xercesc::XMLString::equals(higher_score_better,xercesc::XMLString::transcode("true")))
			{
				pep_id_.setHigherScoreBetter(true);	
			}
			else if (xercesc::XMLString::equals(higher_score_better,xercesc::XMLString::transcode("false")))
			{
				pep_id_.setHigherScoreBetter(false);					
			}
			else
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", "Invalid value for 'higher_score_better '");				
			}
			
			//MZ
			DoubleReal tmp2=-numeric_limits<DoubleReal>::max();
			optionalAttributeAsDouble(tmp2, attributes,"MZ");
			if (tmp2!=-numeric_limits<DoubleReal>::max())
			{
				pep_id_.setMetaValue("MZ", tmp2);
			}
			//RT
			tmp2=-numeric_limits<DoubleReal>::max();
			optionalAttributeAsDouble(tmp2, attributes,"RT");
			if (tmp2!=-numeric_limits<DoubleReal>::max())
			{
				pep_id_.setMetaValue("RT", tmp2);
			}
			Int tmp3=-numeric_limits<Int>::max();
			optionalAttributeAsInt(tmp3, attributes,"spectrum_reference");
			if (tmp3!=-numeric_limits<Int>::max())
			{
				pep_id_.setMetaValue("spectrum_reference", tmp3);				
			}
			
			last_meta_ = &pep_id_;
		}
		else if (element == "PeptideHit")
		{
			pep_hit_ = PeptideHit();
			
			pep_hit_.setCharge(attributeAsInt(attributes,"charge"));
			pep_hit_.setScore(attributeAsDouble(attributes,"score"));
			pep_hit_.setSequence(attributeAsString(attributes,"sequence"));
			
			//aa_before
			String tmp="";
			optionalAttributeAsString(tmp,attributes,"aa_before");
			if (!tmp.empty())
			{
				pep_hit_.setAABefore(tmp[0]);
			}
			//aa_after
			tmp="";
			optionalAttributeAsString(tmp,attributes,"aa_after");
			if (!tmp.empty())
			{
				pep_hit_.setAAAfter(tmp[0]);
			}
			
			//parse optional protein ids to determine accessions
			const XMLCh* refs = attributes.getValue(xercesc::XMLString::transcode("protein_refs"));
			if (refs!=0)
			{
				String accession_string = xercesc::XMLString::transcode(refs);
				accession_string.trim();
				vector<String> accessions;
				accession_string.split(' ', accessions);
				if (accession_string!="" && accessions.size()==0)
				{
					accessions.push_back(accession_string);
				}
				for(vector<String>::const_iterator it = accessions.begin(); it!=accessions.end(); ++it)
				{
					map<String,String>::const_iterator it2 = proteinid_to_accession_.find(*it);
					if (it2!=proteinid_to_accession_.end())
					{
						pep_hit_.addProteinAccession(it2->second);
					}
					else
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Invalid protein reference '") + *it + "'" );
					}
				}
			}
			last_meta_ = &pep_hit_;
		}
		
		//USERPARAM
		else if (element == "UserParam")
		{
			if (last_meta_ == 0)
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", "UserParam unexpected!" );
			}

			static XMLCh* s_name = xercesc::XMLString::transcode("name");
			static XMLCh* s_value = xercesc::XMLString::transcode("value");
			static XMLCh* s_type = xercesc::XMLString::transcode("type");
			static XMLCh* s_int = xercesc::XMLString::transcode("int");
			static XMLCh* s_float = xercesc::XMLString::transcode("float");
			static XMLCh* s_string = xercesc::XMLString::transcode("string");	
			
			const XMLCh* value = attributes.getValue(s_value);
			const XMLCh* type = attributes.getValue(s_type);
			
			//register name
			String name = xercesc::XMLString::transcode(attributes.getValue(s_name));
			last_meta_->metaRegistry().registerName(name,"","");
			
			if(*type==*s_int)
			{
				last_meta_->setMetaValue(name, xercesc::XMLString::parseInt(value));
			}
			else if (*type==*s_float)
			{
				last_meta_->setMetaValue(name, atof(xercesc::XMLString::transcode(value)) );
			}
			else if (*type==*s_string)
			{
				last_meta_->setMetaValue(name, (String)xercesc::XMLString::transcode(value));
			}
			else
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Invlid UserParam type '") + xercesc::XMLString::transcode(type) + "'" );
			}
		}
	}
	
	void IdXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		String element = xercesc::XMLString::transcode(qname);
		
		//cout << "End: " << element << endl;
		
		//START
		if (element == "IdXML")
		{
			
		}

		///SEARCH PARAMETERS
		else if (element == "SearchParameters")
		{
			last_meta_ = 0;
			parameters_[id_] = param_;
		}		
		else if (element == "FixedModification")
		{
			last_meta_ = &param_;
		}
		else if (element == "VariableModification")
		{
			
			last_meta_ = &param_;
		}
		
		// RUN
		else if (element == "IdentificationRun")
		{

		}
		
		//PROTE IDENTIFICATIONS
		else if (element == "ProteinIdentification")
		{
			prot_ids_->push_back(prot_id_);
			prot_id_ = ProteinIdentification();
			last_meta_  = 0;		
		}
		else if (element == "ProteinHit")
		{
			prot_id_.insertHit(prot_hit_);
			last_meta_ = &prot_id_;
		}
		
		//PEPTIDES
		else if (element == "PeptideIdentification")
		{
			pep_ids_->push_back(pep_id_);
			pep_id_ = PeptideIdentification();
			last_meta_  = 0;
		}
		else if (element == "PeptideHit")
		{
			pep_id_.insertHit(pep_hit_);
			last_meta_ = &pep_id_;
		}
	}

	void IdXMLFile::writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent) const
	{
		std::vector<String> keys;
		meta.getKeys(keys);
		
		for (UInt i = 0; i!=keys.size();++i)
		{
			if (keys[i]!="MZ" && keys[i]!="RT" && keys[i]!="spectrum_reference")
			{
				os << String(2*indent,' ') << "<UserParam type=\"";
				
				DataValue d = meta.getMetaValue(keys[i]);
				//determine type
				if (d.valueType()==DataValue::STRVALUE)
				{
					os << "string\" name=\"" << keys[i] << "\" value=\"" << (String)(d) << "\"/>" << endl;
				}
				if (d.valueType()==DataValue::INTVALUE || d.valueType()==DataValue::SHOVALUE || d.valueType()==DataValue::LONVALUE)
				{
					os << "int\" name=\"" << keys[i] << "\" value=\"" << (String)(d) << "\"/>" << endl;
				}
				if (d.valueType()==DataValue::DOUVALUE || d.valueType()==DataValue::FLOVALUE)
				{
					os << "double\" name=\"" << keys[i] << "\" value=\"" << (String)(d) << "\"/>" << endl;
				}
			}
		}
	}


} // namespace OpenMS
