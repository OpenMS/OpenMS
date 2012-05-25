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
				// handle version and experiment type
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
				//~ ...
		}


		void MzQuantMLHandler::writeTo(std::ostream&  os)
		{
			//~ TODO logger_.startProgress(0,exp.size(),"storing mzML file");
			String line; //everyone walk the line!!!
			std::vector<UInt64> rid;
 			
			//header
			//~ TODO CreationDate
			os << "<mzQuantML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psidev.info/psi/pi/mzQuantML/1.0.0-rc2 ../../schema/mzQuantML_1_0_0-rc2.xsd\" xmlns=\"http://psidev.info/psi/pi/mzQuantML/1.0.0-rc2\"" << " version=\"1.0.0-rc2\"" << ">\n";
			
			//CVList
			os << "<CvList>\n \t<Cv id=\"PSI-MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Vocabularies\"  uri=\"http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"2.25.0\"/>\n\t<Cv id=\"UO\" fullName=\"Unit Ontology\" uri=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\"/>\n</CvList>\n";
			
			// TODO write analysissummary
			
			//Software & DataProcessing
			String softwarelist_tag;
			softwarelist_tag += "\t<SoftwareList>\n"	;
				
			String dataprocessinglist_tag;
			dataprocessinglist_tag += "\t<DataProcessingList>\n";
			// TODO Software DefaultTag for each file: OpenMS
			Size order_d = 0;

			std::vector<DataProcessing> pl = cmsq_->getDataProcessingList();
			for (std::vector<DataProcessing>::const_iterator dit = pl.begin(); dit != pl.end(); ++dit)
			//~ for (std::vector<DataProcessing>::const_iterator dit = cmsq_->getDataProcessingList().begin(); dit != cmsq_->getDataProcessingList().end(); ++dit) // soome wierd bug is making this impossible resulting in segfault - to tired to work this one out right now
			{
				String sw_ref;
				sw_ref = "sw_" + String(UniqueIdGenerator::getUniqueId());
				softwarelist_tag += "\t\t<Software id=\"" +  sw_ref + "\" version=\"" + String(dit->getSoftware().getVersion()) + "\"/>\n";
				writeCVParams_(softwarelist_tag, dit->getSoftware().getCVTerms(), UInt(3));
				if (dit->getSoftware().getCVTerms().empty())
				{
					softwarelist_tag += "\t\t\t<userParam name=\""+ String(dit->getSoftware().getName()) +"\"/>\n";
				}
				softwarelist_tag += "\t\t</Software>\n";
				++order_d;		
				dataprocessinglist_tag += "\t\t<DataProcessing id=\"dp_" + String(UniqueIdGenerator::getUniqueId()) + "\" software_Ref=\"" + sw_ref + "\" order=\"" + String(order_d) + "\">\n";
				Size order_c = 0;
				for (std::set<DataProcessing::ProcessingAction>::const_iterator pit = dit->getProcessingActions().begin(); pit != dit->getProcessingActions().end(); ++pit)
				{
					//~ TODO rewrite OpenMS::DataProcessing
					//~ TODO add CVTermList/MetaInfoInterfaceObject to DataProcessing and ParamGroup/Order to "ProcessingAction" or document implicit ordering
					++order_c;
					dataprocessinglist_tag += "\t\t\t<ProcessingMethod order=\"" + String(order_c) + "\">\n";
					//~ writeUserParam_(dataprocessinglist_tag, pit->getUserParams(), UInt(4));  //writeUserParam_(String& s, const MetaInfoInterface& meta, UInt indent)
					//~ writeCVParams_(dataprocessinglist_tag, (pit->getCVParams.getCVTerms(), UInt(4));  //writeCVParams_(String& s, const Map< String, std::vector < CVTerm > > & , UInt indent)
					dataprocessinglist_tag += "\t\t\t\t<userParam name=\"" + String(DataProcessing::NamesOfProcessingAction[*pit]) + "\" value=\"" + String(dit->getSoftware().getName()) + "\" />\n";
					dataprocessinglist_tag += "\t\t\t</ProcessingMethod>\n";
				}
				dataprocessinglist_tag += "\t\t</DataProcessing>\n";
			}
				
			dataprocessinglist_tag += "\t</DataProcessingList>\n";
				
			softwarelist_tag += "\t</SoftwareList>\n";
				
			os << dataprocessinglist_tag << softwarelist_tag;
			
			// Ratios tag
			String ratio_xml("\t<RatioList>\n");
			for (std::vector<MSQuantifications::Assay>::const_iterator ait = cmsq_->getAssays().begin()+1; ait != cmsq_->getAssays().end(); ++ait)
			{
				UInt64 ri = UniqueIdGenerator::getUniqueId();
				rid.push_back(ri);
				ratio_xml += "\t\t<Ratio id=\"r_" + String(ri) +" numerator_ref=a_\""+ String(cmsq_->getAssays().begin()->uid_) +"\" denominator_ref=a_\""+ String(ait->uid_) +"\" >\n";
				ratio_xml += "\t\t\t<RatioCalculation>\n\t\t\t\t<userParam name=\"Simple ratio calc\"/>\n\t\t\t\t<userParam name=\"light to medium/.../heavy\"/>\n\t\t\t</RatioCalculation>\n\t\t</Ratio>\n";
			}
			ratio_xml += "\t</RatioList>\n";
			
			// Assay & StudyVariables: each  "channel" gets its assay
			String assay_xml("\t<AssayList>\n"), study_xml("\t<StudyVariableList>\n"), inputfiles_xml("\t<InputFiles>\n");
			std::map<String,String> files;
			for (std::vector<MSQuantifications::Assay>::const_iterator ait = cmsq_->getAssays().begin(); ait != cmsq_->getAssays().end(); ++ait)
			{
				String rfgr,ar,vr;
				rfgr = String(UniqueIdGenerator::getUniqueId());
				vr = String(UniqueIdGenerator::getUniqueId());
				//TODO regroup at Rawfilesgroup level
				String rgs;
				bool group_exists = true;
				rgs += "\t\t<RawFilesGroup id=\"rfg_" + rfgr + "\">\n";
				for (std::vector<ExperimentalSettings>::const_iterator iit = ait->raw_files_.begin(); iit != ait->raw_files_.end(); ++iit)
				{
					if (files.find(iit->getLoadedFilePath()) == files.end())
					{
						group_exists = false;
						UInt64 rid = UniqueIdGenerator::getUniqueId();
						files.insert(std::make_pair<String,String>(iit->getLoadedFilePath(),rfgr));
						rgs += "\t\t\t<RawFile id=\"r_" +String(rid)  + "\" location=\"" + iit->getLoadedFilePath() + "\"/>\n";
						// TODO write proteowizards sourcefiles (if there is any mentioning of that in the mzml) into OpenMS::ExperimentalSettings of the exp
					}
					else
					{
						rfgr = "rfg_" + String(files.find(iit->getLoadedFilePath())->second);
					}
				}
				rgs += "\t\t</RawFilesGroup>\n";
				
				if(!group_exists)
				{
					inputfiles_xml += rgs;
				}
				
				assay_xml += "\t\t<Assay id=\"a_" + String(ait->uid_)  + "\" RawFilesGroup_refs=\"" + rfgr + "\">\n";
				assay_xml += "\t\t\t<Label>\n";
				switch(cmsq_->getAnalysisSummary().quant_type_)
				{ //enum QUANT_TYPES {MS1LABEL=0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES}; // derived from processing applied
					case 0:
						for (std::vector< std::pair<String, DoubleReal> >::const_iterator lit = ait->mods_.begin(); lit != ait->mods_.end(); ++lit)
						{ 	//~ TODO CVTerms
							assay_xml += "\t\t\t\t<Modification massDelta=\""+String(lit->second)+"\" residues=\"TODO\">\n";
							assay_xml += "\t\t\t\t\t<cvParam cvRef=\"TODO\" accession=\"TODO\" name=\""+String(lit->first)+"\" value=\"0\"/>\n";							
						}
					break;
					
					case 1:
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
			os << inputfiles_xml << assay_xml << study_xml << ratio_xml;

			// Features and QuantLayers
			std::vector<UInt64> fid;
			std::vector<Real> fin, fwi /*, fqu */;
			std::vector< std::vector< std::vector<UInt64> >  >cid; //per consensusmap - per consensus - per feature (first entry is consensus idref)
			std::vector< std::vector<Real> > f2i;
			String feature_xml = "";
			feature_xml += "\t<FeatureList>\n"; // TODO spectrum_refs
			for (std::vector<ConsensusMap>::const_iterator mit = cmsq_->getConsensusMaps().begin(); mit != cmsq_->getConsensusMaps().end(); ++mit)
			{
				std::vector< std::vector<UInt64> > cmid;
				for (ConsensusMap::const_iterator cit = mit->begin(); cit != mit->end(); ++cit)
				{
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
							cmid.push_back(idvec);
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
				}
				cid.push_back(cmid);
			}
			os << feature_xml;
				
			switch (cmsq_->getAnalysisSummary().quant_type_) //enum QUANT_TYPES {MS1LABEL=0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES}; // derived from processing applied
			{
				case 0: //ms1label
				{
					os << "\t\t<FeatureQuantLayer id=\"" << "q_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n\t\t\t<ColumnDefinition>\n";

					os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"intensity\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>";
					os << "\t\t\t\t<Column index=\"1\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"width\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>";
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
								os << "\t\t<PeptideConsensus id=\"" << "c_" << String(cid[k][i].front()) << "\" charge=\""+ String((*cmsq_).getConsensusMaps()[k][i].getCharge()) +"\">\n";
								os << "\t\t\t<feature_refs>";
								for (Size j=1; j < cid[k][i].size(); ++j)
								{
									os << "f_" << cid[k][i][j]<< " ";
								}
								os << "\t\t\t</feature_refs>\n";
								if (!(*cmsq_).getConsensusMaps()[k][i].getPeptideIdentifications().empty())
								{									
									//~ os << "\t\t\t<IdentificationRef id_ref=\"";
									//~ os << (*cmsq_).getConsensusMaps()[k][i].getPeptideIdentifications().front().getIdentifier() << "\" IdentificationFile_ref=\"";
									//~ os << idid_to_idfilenames.begin()->first  << "\"/>\n";
								}
								os << "\t\t</PeptideConsensus>\n";
							}
							
							// QuantLayers
							os << "\t\t<RatioQuantLayer id=\"" << "q_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n";
							os << "\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001132\" name=\"peptide ratio\"/>\n\t\t\t\t\t</DataType>\n";
							os << "\t\t\t\t<Column index>";
							//~ todo <ColumnIndex>ratio_L_M ... ratio_L_H</ColumnIndex>
							os << "</ColumnIndex>\n\t\t\t<DataMatrix>\n";

							//~ calculate ratios
							for (Size i=0; i < cid[k].size(); ++i)
							{
								const std::set< FeatureHandle,FeatureHandle::IndexLess>& feature_handles = (*cmsq_).getConsensusMaps()[k][i].getFeatures();
								String dis;
								if (feature_handles.size() > 1)
								{
									os <<"\t\t\t\t<Row object_ref=\"c_" << String(cid[k][i].front()) << "\">";
									std::set< FeatureHandle,FeatureHandle::IndexLess>::const_iterator fit = feature_handles.begin(); // this is unlabeled
									fit++;
									for (fit; fit != feature_handles.end(); ++fit)
									{
										dis += String( feature_handles.begin()->getIntensity() / fit->getIntensity() );
										dis += " ";
									}
									os << dis.trim() << "</Row>\n";
								}
							}
							os << "\t\t\t</DataMatrix>\n";
							os << "\t\t</RatioQuantLayer>\n";
						}
						break;
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
						break;
				}
				os << "\t</PeptideList>\n";
			}

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
			if (meta.isMetaEmpty())
			{
				return;
			}				
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
