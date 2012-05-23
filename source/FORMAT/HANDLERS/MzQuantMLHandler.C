// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzQuantMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <set>
#include <vector>
#include <map>

using namespace std;

// TODO getter/setter for MSQuantifications
// TODO Software DefaultTag for each file: OpenMS

namespace OpenMS
{
	namespace Internal
	{
		MzQuantMLHandler::MzQuantMLHandler(const MSQuantifications& msq, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			msq_(0),
			cmsq_(&msq)
		{
				cv_.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
		}

		MzQuantMLHandler::MzQuantMLHandler(MSQuantifications& msq, /* FeatureMap& feature_map, */ const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
			logger_(logger),
			msq_(&msq),
			cmsq_(0)
		{
				cv_.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
		}

		MzQuantMLHandler::~MzQuantMLHandler()
		{
		}
		
		void MzQuantMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
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

			if (tag_ == "MzQuantML")
			{
				// TODO handle version and experiment type
				return;
			}

			//~ if (tag_ == "SourceFile")
			//~ {
				//~ // start new
				//~ actual_inputfile_ = MSQuantifications::Inputfile();
				
				//~ // location attribute
				//~ actual_inputfile_.location = attributeAsString_(attributes, "location");

				//~ // name attribute (opt)
				//~ String name;
				//~ if (optionalAttributeAsString_(name, attributes, "name"))
				//~ {
					//~ actual_inputfile_.name = name;
				//~ }
				
				// ext doc uri element in cahracter tag ExternalFormatDocumentation
				// FileFormat is own tag -> CVparam
				//~ return;
			//~ }
			error(LOAD, "MzIdentMLHandler::startElement: Unkown element found: '" + tag_ + "' in tag '" + parent_tag + "', ignoring.");
		}

		void MzQuantMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
		{
			//~ if (tag_ == "ExternalFormatDocumentation")
			//~ {
				//~ actual_inputfile_.doc_uri = sm_.convert(chars);

				//~ return;
			//~ }

			//error(LOAD, "MzIdentMLHandler::characters: Unkown character section found: '" + tag_ + "', ignoring.");
		}

		void MzQuantMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			//~ static set<String> to_ignore;
			//~ if (to_ignore.empty())
			//~ {
				//~ to_ignore.insert("mzIdentML");
				//~ to_ignore.insert("cvParam");
			//~ }

			//~ tag_ = sm_.convert(qname);
			//~ open_tags_.pop_back();

			//~ if (to_ignore.find(tag_) != to_ignore.end())
			//~ {
				//~ return;
			//~ }

			//~ if (tag_ == "SourceFile")
			//~ {
				//~ inputfiles_.addHit(actual_inputfile_);
				//~ return;
			//~ }
		}

		void MzQuantMLHandler::handleCVParam_(const String& /* parent_parent_tag*/, const String& parent_tag, const String& accession, /* const String& name, */ /* const String& value, */ const xercesc::Attributes& attributes, const String& cv_ref /* , const String& unit_accession */)
		{
				//~ TODO
		}


		void MzQuantMLHandler::writeTo(std::ostream&  os)
		{
			//~ TODO logger_.startProgress(0,exp.size(),"storing mzML file");
			String line; //everyone walk the line!!!
			//header
			//~ TODO CreationDate
			os << "<mzQuantML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psidev.info/psi/pi/mzQuantML/1.0.0-rc2 ../../schema/mzQuantML_1_0_0-rc2.xsd\" xmlns=\"http://psidev.info/psi/pi/mzQuantML/1.0.0-rc2\"" << " version=\"1.0.0-rc2\"" << ">\n";
			
			//CVList
			os << "<CvList>\n \t<Cv id=\"PSI-MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Vocabularies\"  uri=\"http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"2.25.0\"/>\n\t<Cv id=\"UO\" fullName=\"Unit Ontology\" uri=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\"/>\n</CvList>\n";
			
			// Inputfiles TODO
			//~ os << "\t<InputFiles>\n";
			//~ os << "\t\t<RawFilesGroup id=\"ig_"+String(UniqueIdGenerator::getUniqueId())+"\" >\n";
			//~ for (std::vector<InputFile>::const_iterator it = cmsq_.getInputFiles().begin(); it != cmsq_.getInputFiles().end(); ++it)
			//~ {
				//~ String id, location;
				//~ if (*it).id != ""
				//~ {
					//~ id = (*it).id
				//~ }					
				//~ else 
				//~ {
					//~ id = String(UniqueIdGenerator::getUniqueId());
				//~ }
				//~ line = "\t\t\t<RawFile id=\"" + String((*it).second) + "\" location=\"" + String((*it).location) + "\"/>\n";
				//~ os << line;
			//~ }
			//~ os << "\t\t</RawFilesGroup>\n";
				// TODO decide on Inputfiles structure in MSQuantifications
				//~ if conversion from consensusXML:
				//~ inputfile_tag += "\t\t\t<SourceFile id=\"is_" + String(UniqueIdGenerator::getUniqueId()) + "\" location=\"" + cmsq_->getLoadedFilePath() + "\"/>\n";

				//~ std::map<UInt64,String> idid_to_idfilenames;
				//~ for ( Size i = 0; i < cmsq_->getDataProcessing().size(); ++i )
				//~ {
						//~ for (std::set<DataProcessing::ProcessingAction>::const_iterator it = cmsq_->getDataProcessing()[i].getProcessingActions().begin(); it != cmsq_->getDataProcessing()[i].getProcessingActions().end(); ++it)
						//~ {
								//~ if (cmsq_->getDataProcessing()[i].metaValueExists("parameter: id"))
								//~ {
										//~ idid_to_idfilenames.insert(std::pair<UInt64,String>(UniqueIdGenerator::getUniqueId(),String(cmsq_->getDataProcessing()[i].getMetaValue ("parameter: id"))));
								//~ }
						//~ }
				//~ }
				
				//~ inputfile_tag += "\t<IdentificationFiles>\n";
				//~ for (std::map<UInt64,String>::const_iterator it=idid_to_idfilenames.begin(); it != idid_to_idfilenames.end(); ++it)
				//~ {
						//~ inputfile_tag += "\t\t\t<IdentificationFile id=\"ir_" + String(it->first) + "\" location=\"" +  (it->second) + "\"/>\n";
				//~ }
				//~ inputfile_tag += "\t\t</IdentificationFiles>\n";
			//~ os << "\t</InputFiles>\n";
				
			// SoftwareList +	DataProcessing
			String softwarelist_tag;
			softwarelist_tag += "\t<SoftwareList>\n"	;
				
			String dataprocessinglist_tag;
			dataprocessinglist_tag += "\t<DataProcessingList>\n";
			Size order_d = 0;
			for (std::vector<DataProcessing>::const_iterator dit = cmsq_->getDataProcessingList().begin(); dit != cmsq_->getDataProcessingList().end(); ++dit)
			{
				String sw_ref;
				sw_ref = "sw_" + String(UniqueIdGenerator::getUniqueId());
				softwarelist_tag += "\t\t<Software id=\"" +  sw_ref + "\" version=\"" + String(dit->getSoftware().getVersion()) + "\"/>\n";
				writeCVParams_(softwarelist_tag, dit->getSoftware().getCVTerms(), UInt(3));
				writeUserParams_(softwarelist_tag, dit->getSoftware(), UInt(3));
				// TODO if no CV == name write userParam with Software.getName()!!!
				softwarelist_tag += "\t\t</Software>\n";
				++order_d;		
				dataprocessinglist_tag += "\t\t<DataProcessing id=\"dp_" + String(UniqueIdGenerator::getUniqueId()) + "\" software_Ref=\"" + sw_ref + "\" order=\"" + String(order_d) + "\">\n";
				Size order_c = 0;
				for (std::set<DataProcessing::ProcessingAction>::const_iterator pit = dit->getProcessingActions().begin(); pit != dit->getProcessingActions().end(); ++pit)
				{
					//~ TODO rewrite OpenMS::DataProcessing
					++order_c;
					dataprocessinglist_tag += "\t\t\t<ProcessingMethod order=\"" + String(order_c) + "\">\n";
					//~ TODO add CVTermList/MetaInfoInterfaceObject and Order to "ProcessingAction" 
					//~ writeUserParam_(dataprocessinglist_tag, pit->getUserParams(), UInt(4));  //writeUserParam_(String& s, const MetaInfoInterface& meta, UInt indent)
					//~ writeCVParams_(dataprocessinglist_tag, (pit->getCVParams.getCVTerms(), UInt(4));  //writeCVParams_(String& s, const Map< String, std::vector < CVTerm > > & , UInt indent)
					dataprocessinglist_tag += "\t\t\t\t<userParam name=\"" + DataProcessing::NamesOfProcessingAction[*pit] + "\" value=\"" + dit->getSoftware().getName() + "\" />\n";
					dataprocessinglist_tag += "\t\t\t</ProcessingMethod>\n";
				}
				dataprocessinglist_tag += "\t\t</DataProcessing>\n";
			}
				
			dataprocessinglist_tag += "\t</DataProcessingList>\n";
				
			softwarelist_tag += "\t</SoftwareList>\n";
				
			os << dataprocessinglist_tag << softwarelist_tag;

			// Assay & StudyVariables tags

			//  Assay: each  "channel" gets its assay
			String assay_xml("\t<AssayList>\n"), study_xml("\t<StudyVariableList>\n"), inputfiles_xml("\t<InputFiles>\n");
			
			for (std::vector<MSQuantifications::Assay>::const_iterator ait = cmsq_->getAssays().begin(); ait != cmsq_->getAssays().end(); ++ait)
			{
				String rfgr,ar,vr;
				rfgr = String(UniqueIdGenerator::getUniqueId());
				ar = String(UniqueIdGenerator::getUniqueId());
				vr = String(UniqueIdGenerator::getUniqueId());
				
				for (std::vector<ExperimentalSettings>::const_iterator iit = ait->raw_files_.begin(); iit != ait->raw_files_.end(); ++iit)
				{
					inputfiles_xml += "\t\t<RawFilesGroup id=\"rfg_" + rfgr + "\">\n";
					// TODO do better grouping in MSQuantifications class' Assay
					for (std::vector<SourceFile>::const_iterator sit = iit->getSourceFiles().begin(); sit != iit->getSourceFiles().end(); ++sit)
					{
						inputfiles_xml += "\t\t\t<RawFile id=\"" + String(UniqueIdGenerator::getUniqueId()) + "\" location=\"" + sit->getPathToFile() + sit->getNameOfFile() + "\"/>\n";
					}
					inputfiles_xml += "\t\t</RawFilesGroup>\n";
				}
				
				assay_xml += "\t\t<Assay id=\"a_" + ar + "\" RawFilesGroup_refs=\"" + rfgr + "\">\n";
				assay_xml += "\t\t\t<Label>\n";
				switch(ait->label_type_)
				{ //NONE=0, SILAC_LIGHTER, SILAC_LIGHT, SILAC_MEDIUM, SILAC_HEAVY, ITRAQ114, ITRAQ115,ITRAQ116, ITRAQ117, SIZE_OF_LABEL_TYPES
					case 2:
					// TODO read the mod masses somewhere from params - no hardcode!   
						//~ TODO CVTerms
						//~ [term]
						//~ id: MS:1001819
						//~ name: two sample run
						assay_xml += "\t\t\t\t<Modification massDelta=\"0\" residues=\"L\">\n";
						assay_xml += "\t\t\t\t\t<cvParam cvRef="" accession="" name=\"None\" value=\"0\"/>\n";
					break;
					case 4:
						assay_xml += "\t\t\t\t<Modification massDelta=\"8.0141988132\" residues=\"L\">\n";
						assay_xml += "\t\t\t\t\t<cvParam name=\"Lys8\" value=\"8.0141988132\"/>\n";
					break;
					case 5:
					case 6:
					case 7:
					case 8:
					{
						assay_xml += "\t\t\t\t<Modification massDelta=\"145\" residues=\"N-term\">\n";
						assay_xml += "\t\t\t\t\t<cvParam name =\"itraq label\"/>\n";
						break;
					}
					default:
						assay_xml += "\t\t\t\t<Modification massDelta=\"0\" residues=\"X\">\n";
						assay_xml += "\t\t\t\t\t<cvParam name =\"no label\"/>\n";
				}
				assay_xml += "\t\t\t\t</Modification>\n";
				assay_xml += "\t\t\t</Label>\n";
				assay_xml += "\t\t</Assay>\n";

				// for SILACAnalyzer/iTRAQAnalyzer one assay is one studyvariable, this may change!!! TODO for iTRAQ
				study_xml += "\t<StudyVariable id=\"v_" + vr + "\" name=\"String\">\n";
				study_xml += "\t\t\t<AssayRef assay_ref=\"a_" + ar + "\"/>\n";
				study_xml += "\t</StudyVariable>\n";
			}
			assay_xml += "</AssayList>\n";
			inputfiles_xml += "</InputFiles>\n";
			study_xml += "</StudyVariableList>\n";
			os << inputfiles_xml << assay_xml << study_xml;

			// Features and QuantLayers
			std::vector<UInt64> fid;
			std::vector<Real> fin, fwi /*, fqu */;
			std::vector< std::vector< std::vector<UInt64> >  >cid; //per consensusmap - per consensus - per feature (first entry is consensus idref)
			std::vector< std::vector<Real> > f2i;
			String feature_xml = "";
			feature_xml += "\t<FeatureList>\n"; // TODO spectrum_refs
			for (std::vector<ConsensusMap>::const_iterator mit = cmsq_->getConsensusMaps().begin(); mit != cmsq_->getConsensusMaps().end(); ++mit)
			{
				for (ConsensusMap::const_iterator cit = mit->begin(); cit != mit->end(); ++cit)
				{
					std::vector< std::vector<UInt64> > cvec;
					const std::set< FeatureHandle,FeatureHandle::IndexLess>& feature_handles = cit->getFeatures();
					switch (cmsq_->getAnalysisSummary().quant_type_) //enum QUANT_TYPES {MS1LABEL=0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES}; // derived from processing applied
					{
						case 0: //ms1label
						{
							std::vector<UInt64> idvec;
							idvec.push_back(UniqueIdGenerator::getUniqueId());
							for (std::set< FeatureHandle,FeatureHandle::IndexLess>::const_iterator fit = feature_handles.begin(); fit != feature_handles.end(); ++fit)
							{
								fid.push_back(UniqueIdGenerator::getUniqueId());
								idvec.push_back(fid.back());
								fin.push_back(fit->getIntensity());
								fwi.push_back(fit->getWidth());
								//~ fqu.push_back(jt->getQuality());
								feature_xml += "\t\t<Feature id=\"f_" + String(fid.back()) + "\" RT=\"" + String(fit->getRT()) + "\" MZ=\"" + String(fit->getMZ()) + "\" charge=\"" + String(fit->getCharge()) + "\"/>\n";
								// TODO as soon as SILACanalyzer incorporate convex hulls read from the featuremap
								//~ writeUserParam_(os, *jt, UInt(2)); // FeatureHandle has no MetaInfoInterface!!!
							}
							cvec.push_back(idvec);
						}break;
						case 1: //ms2label
						{
							std::vector<Real> fi;
							fid.push_back(UniqueIdGenerator::getUniqueId());
							feature_xml += "\t\t<Feature id=\"f_" + String(fid.back()) + "\" RT=\"" + String(cit->getRT()) + "\" MZ=\"" + String(cit->getMZ()) + "\" charge=\"" + String(cit->getCharge()) + "\"/>\n";
							//~ std::vector<UInt64> cidvec;
							//~ cidvec.push_back(fid.back());
							for (std::set< FeatureHandle,FeatureHandle::IndexLess>::const_iterator fit = feature_handles.begin(); fit != feature_handles.end(); ++fit)
							{
								fi.push_back(fit->getIntensity());
							}
							f2i.push_back(fi);
						}break;
					}
					cid.push_back(cvec);
				}
			}
			os << feature_xml;
				
			switch (cmsq_->getAnalysisSummary().quant_type_) //enum QUANT_TYPES {MS1LABEL=0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES}; // derived from processing applied
			{
				case 0: //ms1label
				{
					os << "\t\t<FeatureQuantLayer id=\"" << "q_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n\t\t\t<ColumnDefinition>\n";

					os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"intensity\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>";
					os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"width\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>";
					//~ os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"quality\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>"; // getQuality erst ab BaseFeature - nicht in FeatureHandle
					os << "</ColumnDefinition>\t\t\t\t\n<DataMatrix>\n";
					for (Size i=0; i < fid.size(); ++i)
					{
						os <<"\t\t\t\t\t<Row object_ref=\"f_" << String(fid[i]) << "\">";
						os << fin[i] << " " << fwi[i] /* << " " << fiq[i] */;
						os << "</Row>\n";
					}
					os << "\t\t\t</DataMatrix>\n";
					os << "\t\t</FeatureQuantLayer>\n";
				}
				break;
				case 1: //ms2label
				{
					os << "\t\t<MS2AssayQuantLayer id=\""+ String(UniqueIdGenerator::getUniqueId()) +"\">\n\t\t\t<DataType>\n\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"intensity\"/>\n\t\t\t</DataType>\n\t\t\t<ColumnDefinition>\n\t\t\t\t\t";
					//~ for (Size i = 0; i < file_descriptions.size(); ++i)
					//~ {
						//~ os << "TODO";
					//~ }
					os <<  "\n\t\t\t\t</ColumnDefinition>\t\t\t<DataMatrix>\n";
					for (Size i = 0; i < fid.size(); ++i)
					{
						os << "\t\t\t\t\t<Row object_ref=\"f_" + String(fid[i]) + "\">";
						for (Size j = 0; j < f2i[i].size(); ++j)
						{
							os << String(f2i[i][j]) << " ";
						}
						os << "</Row>\n";
					}
					os << "\t\t\t</DataMatrix>\n\t\t</MS2AssayQuantLayer>\n";
				}
				break;
			}
			
			os << "\t</FeatureList>\n";

			// Peptides
			os << "\t<PeptideList  id=\"" << "m_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n";
			
			for (Size k = 0; k < cid.size(); ++k)
			{
				switch (cmsq_->getAnalysisSummary().quant_type_) //enum QUANT_TYPES {MS1LABEL=0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES}; // derived from processing applied
				{
						case 0: // ms1label - iterate consensusmap?
						{
							for (Size i = 0; i < cid[k].size(); ++i)
							{
								//~ TODO also report consensus without ids
								if (!(*cmsq_).getConsensusMaps()[k][i].getPeptideIdentifications().empty())
								{
									os << "\t\t<PeptideConsensus id=\"" << "c_" << String(cid[k][i].front()) << "\" charge=\""+ String((*cmsq_).getConsensusMaps()[k][i].getCharge()) +"\">\n";
									os << "\t\t\t<feature_refs>";
									for (Size j=1; j < cid[k][i].size(); ++j)
									{
										os << "f_" << cid[k][i][j]<< " ";
									}
									os << "\t\t\t</feature_refs>\n";
									//~ os << "\t\t\t<IdentificationRef id_ref=\"";
									//~ os << (*cmsq_).getConsensusMaps()[k][i].getPeptideIdentifications().front().getIdentifier() << "\" IdentificationFile_ref=\"";
									//~ os << idid_to_idfilenames.begin()->first  << "\"/>\n";
									os << "\t\t</PeptideConsensus>\n";
								}
							}
								// QuantLayers
								os << "\t\t<RatioQuantLayer id=\"" << "q_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n\t\t\t<ColumnIndex>\n";

								os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"ratio\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>";
								//~ os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"quality\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>"; // getQuality erst ab BaseFeature - nicht in FeatureHandle
								os << "</ColumnIndex>\t\t\t\t\n<DataMatrix>\n";
								for (Size i=0; i < cid[k].size(); ++i)
								{
										os <<"\t\t\t\t\t<Row object_ref=\"c_" << String(cid[k][i].front()) << "\">";
										os << String((*cmsq_).getConsensusMaps()[k][i].getIntensity());
										os << "</Row>\n";
								}
								os << "\t\t\t</DataMatrix>\n";
								os << "\t\t</RatioQuantLayer>\n";
						}
						case 1: // ms2label - iterate featuremap
						{
								for (Size i = 0; i < fid.size(); ++i)
								{
										if (!(*cmsq_).getConsensusMaps()[k][i].getPeptideIdentifications().empty())
										{
												os << "\t\t<PeptideConsensus id=\"" << "c_" << String(fid[i]) << "\" charge=\""+ String((*cmsq_).getConsensusMaps()[k][i].getCharge()) +"\">\n";
												os << "\t\t\t<feature_refs>";
												for (Size j = 0; j < f2i[i].size(); ++j)
												{
														os << String(f2i[i][j]) << " ";
												}
												os << "\t\t\t</feature_refs>\n";
												//~ os << "\t\t\t<IdentificationRef id_ref=\"";
												//~ os << (*cmsq_).getConsensusMaps()[k][i].getPeptideIdentifications().front().getIdentifier() << "\" IdentificationFile_ref=\"";
												//~ os << idid_to_idfilenames.begin()->first  << "\"/>\n";
												os << "\t\t</PeptideConsensus>\n";
										}
								}
								//~ TODO ratios, when available (not yet for the iTRAQ tuples of iTRAQAnalyzer)
						}
				}
				os << "\t</PeptideList>\n";
			}


				//--------------------------------------------------------------------------------------------
				// Ratio
				//--------------------------------------------------------------------------------------------
				//~ if (experiment_type == 0)
				//~ {
						//~ os << "\t<RatioList>\n";
						//~ for (Size i = 0; i < assayid_to_studyvarid.size(); ++i)
						//~ {
								//~ os << "\t\t<Ratio id=\"" << "r_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n";
								//~ //todo paarungen
								//~ os << "\t\t</Ratio>\n";
						//~ }
						//~ os << "\t</RatioList>\n";
				//~ }

				//--------------------------------------------------------------------------------------------
				// Proteins and Proteingroups
				//--------------------------------------------------------------------------------------------
				// TODO - omitted as there are no ids yet

				os << "</mzQuantML>\n";
		}

		void MzQuantMLHandler::writeCVParams_(String& s, const Map< String, std::vector < CVTerm > > & cvl, UInt indent)
		{
			String inden ((size_t)indent, '\t');
			for (std::map<String,std::vector<CVTerm> >::const_iterator jt = cvl.begin(); jt != cvl.end(); ++jt)
			{
				for (std::vector<CVTerm>::const_iterator kt =  (*jt).second.begin(); kt !=  (*jt).second.end(); ++kt)
				{
					s += inden;
					s += "<cvParam cvRef=\"" + kt->getCVIdentifierRef()+ "\" accession=\"" + (*jt).first + "\" name=\"" + kt->getName() ;
					if ( kt->hasValue() )
					{
						s += "\" value=\"" + kt->getValue().toString() + "\"/>\n"; // value is OpenMS::DataValue
					}
					else
					{
						s +=	 "\"/>\n";
					}
				}
			}
		}

		void MzQuantMLHandler::writeUserParams_(std::ostream& os, const MetaInfoInterface& meta, UInt indent)
		{
			String h;
			writeUserParams_(h, meta, indent);
			os << h;
		}
		
		void MzQuantMLHandler::writeUserParams_(String& s, const MetaInfoInterface& meta, UInt indent)
		{
			std::vector<String> keys;
			meta.getKeys(keys);

			for (Size i = 0; i!=keys.size();++i)
			{
				s += String(indent,'\t') + "<userParam name=\"" + keys[i] + "\" unitName=\"";

				DataValue d = meta.getMetaValue(keys[i]);
				//determine type
				if (d.valueType()==DataValue::INT_VALUE)
				{
					s += "xsd:integer";
				}
				else if (d.valueType()==DataValue::DOUBLE_VALUE)
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

		void MzQuantMLHandler::writeFeature_(ostream& os, const String& identifier_prefix, UInt64 identifier, UInt indentation_level)
		{
				//~ String indent = String(indentation_level,'\t');

				//~ os << indent << "\t\t<Feature id=\"" << identifier_prefix << identifier ;
						//~ os << " RT=" << feat.getRT() << " MZ" << feat.getMZ();
						//~ os << " charge=" << feat.getCharge() << "\">\n";

				//~ // TODO dataprocessing_ref & assay_ref

				//~ os << "\t\t\t<Masstrace>";
				//~ vector<ConvexHull2D> hulls = feat.getConvexHulls();
				//~ vector<ConvexHull2D>::iterator citer = hulls.begin();
				//~ Size hulls_count = hulls.size();

				//~ for (Size i = 0;i < hulls_count; i++)
				//~ {
						//~ ConvexHull2D current_hull = hulls[i];
						//~ current_hull.compress();
						//~ Size hull_size	= current_hull.getHullPoints().size();

						//~ for (Size j=0;j<hull_size;j++)
						//~ {
								//~ DPosition<2> pos = current_hull.getHullPoints()[j];
								//~ /*Size pos_size = pos.size();
								//~ os << indent << "\t\t\t\t<hullpoint>\n";
								//~ for (Size k=0; k<pos_size; k++)
								//~ {
										//~ os << indent << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << precisionWrapper(pos[k]) << "</hposition>\n";
								//~ }
								//~ os << indent << "\t\t\t\t</hullpoint>\n";*/
								//~ os << indent << precisionWrapper(pos[0]) << " " << precisionWrapper(pos[1]) << " ";
						//~ }
						//~ os << "\n";
				//~ }
				//~ os << "</Masstrace>\n";

				//~ os << indent << "<userParam name=\"overallquality\" value=\"" << precisionWrapper(feat.getOverallQuality()) << "\" unitName=\"quality\""/>

				//~ /*TODO write model description as a userParam
				//~ ModelDescription<2> desc = feat.getModelDescription();
				//~ if (!desc.getName().empty() || !desc.getParam().empty())
				//~ {
					//~ os << indent << "\t\t\t<model name=\"" << desc.getName() << "\">\n";
					//~ Param modelp = desc.getParam();
					//~ Param::ParamIterator piter = modelp.begin();
					//~ while (piter != modelp.end())
					//~ {
						//~ os << indent << "\t\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\"/>\n";
						//~ piter++;
					//~ }
					//~ os << indent << "\t\t\t</model>\n";
				//~ } */
				//~ writeUserParam_(os, feat, indentation_level + 3);

				//~ os << indent << "\t\t</feature>\n";

 }




	} //namespace Internal
} // namespace OpenMS
