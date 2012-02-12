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
#include <set>
#include <vector>
#include <map>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{

		MzQuantMLHandler::MzQuantMLHandler(const ConsensusMap& consensus_map, /* const FeatureMap& feature_map, */ const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			cm_(0),
			ccm_(&consensus_map)
		{
				cv_.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
		}

		MzQuantMLHandler::MzQuantMLHandler(ConsensusMap& consensus_map, /* FeatureMap& feature_map, */ const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
			logger_(logger),
			cm_(&consensus_map),
			ccm_(0)
		{
				cv_.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
		}

		MzQuantMLHandler::~MzQuantMLHandler()
		{
		}

		void MzQuantMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
		{
				//~ TODO all
		}

		void MzQuantMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
		{
				//~ TODO all
		}

		void MzQuantMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
				//~ TODO all
		}

		void MzQuantMLHandler::handleCVParam_(const String& /* parent_parent_tag*/, const String& parent_tag, const String& accession, /* const String& name, */ /* const String& value, */ const xercesc::Attributes& attributes, const String& cv_ref /* , const String& unit_accession */)
		{
				//~ TODO
		}


		void MzQuantMLHandler::writeTo(std::ostream&  os)
		{
				//--------------------------------------------------------------------------------------------
				// Header & id generation
				//--------------------------------------------------------------------------------------------
				os << "<!-- This is a converted  example (from OpenMS consensusXML) -->\n";
				os << "<mzQuantML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psidev.info/psi/pi/mzQuantML/0.1 ../../schema/mzQuantML_0_1_7.xsd\" xmlns=\"http://psidev.info/psi/pi/mzQuantML/0.1\">\n";
				os << "<CvList>\n \t<Cv id=\"PSI-MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Vocabularies\"  uri=\"http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"2.25.0\"/>\n\t<Cv id=\"UO\" fullName=\"Unit Ontology\" uri=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\"/>\n</CvList>\n";

				const ConsensusMap::FileDescriptions& file_descriptions = ccm_->getFileDescriptions();
				std::map<String,UInt64> name_to_rawid;
				std::map<UInt64,UInt64> assayid_to_rawid, assayid_to_studyvarid;
				std::map<UInt64,Size> assayid_to_filedescriptionsindex;
				std::map<Size,UInt64> filedescriptionsindex_to_assayid;
				for ( Size i = 0; i < file_descriptions.size(); ++i)
				{
						std::pair<map<String,UInt64>::iterator,bool> ret;
						ret = name_to_rawid.insert(std::pair<String,UInt64> (file_descriptions[i].filename,UniqueIdGenerator::getUniqueId()));
						UInt64 ind(UniqueIdGenerator::getUniqueId());
						assayid_to_rawid.insert(std::pair<UInt64,UInt64> (ind, (*ret.first).second));
						assayid_to_filedescriptionsindex.insert(std::pair<UInt64,Size> (ind, i));
						filedescriptionsindex_to_assayid.insert(std::pair<Size, UInt64> (i,ind));
						// for SILACAnalyzer/iTRAQAnalyzer one assay is one studyvariable, this may change!!!
						assayid_to_studyvarid.insert(std::pair<UInt64,UInt64> (ind, UniqueIdGenerator::getUniqueId()));
				}

				//--------------------------------------------------------------------------------------------
				// Inputfiles
				//--------------------------------------------------------------------------------------------
				String inputfile_tag("<InputFiles>\n");
				inputfile_tag += "\t<RawFilesGroup id=\"ig_"+String(UniqueIdGenerator::getUniqueId())+"\" >\n";
				for (std::map<String,UInt64>::const_iterator it=name_to_rawid.begin(); it != name_to_rawid.end(); ++it)
				{
						inputfile_tag += "\t\t\t<RawFile id=\"ir_" + String((*it).second) + "\" location=\"" + String((*it).first) + "\"/>\n";
				}
				inputfile_tag += "\t\t</RawFilesGroup>\n";

				//~ if conversion from consensusXML:
				inputfile_tag += "\t\t\t<SourceFile id=\"is_" + String(UniqueIdGenerator::getUniqueId()) + "\" location=\"" + ccm_->getLoadedFilePath() + "\"/>\n";

				std::map<UInt64,String> idid_to_idfilenames;
				for ( Size i = 0; i < ccm_->getDataProcessing().size(); ++i )
				{
						for (std::set<DataProcessing::ProcessingAction>::const_iterator it = ccm_->getDataProcessing()[i].getProcessingActions().begin(); it != ccm_->getDataProcessing()[i].getProcessingActions().end(); ++it)
						{
								if (ccm_->getDataProcessing()[i].metaValueExists("parameter: id"))
								{
										idid_to_idfilenames.insert(std::pair<UInt64,String>(UniqueIdGenerator::getUniqueId(),String(ccm_->getDataProcessing()[i].getMetaValue ("parameter: id"))));
								}
						}
				}
				inputfile_tag += "\t<IdentificationFiles>\n";
				for (std::map<UInt64,String>::const_iterator it=idid_to_idfilenames.begin(); it != idid_to_idfilenames.end(); ++it)
				{
						inputfile_tag += "\t\t\t<IdentificationFile id=\"ir_" + String(it->first) + "\" location=\"" +  (it->second) + "\"/>\n";
				}
				inputfile_tag += "\t\t</IdentificationFiles>\n";

				inputfile_tag += "\t</InputFiles>\n";
				os << inputfile_tag;

				QUANT_TYPES experiment_type; // derived from processing applied
				//--------------------------------------------------------------------------------------------
				// SoftwareList and DataProcessing
				//--------------------------------------------------------------------------------------------
				String softwarelist_tag("<SoftwareList>\n");
				String dataprocessinglist_tag("<DataProcessingList>\n");
				for ( Size i = 0; i < ccm_->getDataProcessing().size(); ++i )
				{
						const DataProcessing& processing = ccm_->getDataProcessing()[i];
						UInt64 suid(UniqueIdGenerator::getUniqueId());
						softwarelist_tag += "\t<Software id=\"s_" + String(suid) + "\" version=\"" + String(processing.getSoftware().getVersion()) + "\" />\n";
						softwarelist_tag += "\t\t<cvParam cvRef=\"\" accession=\"\" name=\"OpenMS_1.9\"/>\n";
						Size order_c = 0;
						for (std::set<DataProcessing::ProcessingAction>::const_iterator it = processing.getProcessingActions().begin(); it != processing.getProcessingActions().end(); ++it)
						{
								// TODO data not yet linked in openms core formats
								if(DataProcessing::NamesOfProcessingAction[*it] == String("Quantitation"))
										{
										if (processing.getSoftware().getName()==String("SILACAnalyzer"))
										{
												experiment_type = MS1LABEL;
										}
										else if (processing.getSoftware().getName()==String("ITRAQAnalyzer"))
										{
												experiment_type = MS2LABEL;
										}
										else
										{
												experiment_type = LABELFREE;
										}
								}
								dataprocessinglist_tag += "\t<dataProcessing id=\"d_" + String(UniqueIdGenerator::getUniqueId()) + "\" software_Ref=\"s_" + String(suid) + "\">\n";
								dataprocessinglist_tag += "\t\t<ProcessingMethod order=\"" + String(order_c) + "\">\n";
								dataprocessinglist_tag += "\t\t\t<userParam name=\"" + DataProcessing::NamesOfProcessingAction[*it] + "\" value=\"" + processing.getSoftware().getName() + "\" />\n";
								writeUserParam_(dataprocessinglist_tag, processing, UInt(3));
								dataprocessinglist_tag += "\t\t</ProcessingMethod>\n";
								dataprocessinglist_tag += "\t</dataProcessing>\n";
								++order_c;
						}
				}
				softwarelist_tag += "\t</SoftwareList>\n";
				dataprocessinglist_tag += "\t</DataProcessingList>\n";
				os << softwarelist_tag << dataprocessinglist_tag;

				//--------------------------------------------------------------------------------------------
				// Assay & StudyVariables
				//--------------------------------------------------------------------------------------------
				//  unique id dispatch only here - String(UniqueIdGenerator::getUniqueId())
				String assay_xml("<AssayList>\n"), study_xml("<StudyVariableList>\n");
				for (std::map<UInt64,UInt64>::const_iterator it= assayid_to_rawid.begin(); it != assayid_to_rawid.end(); ++it)
				{
						assay_xml += "\t<Assay id=\"a_" + String((*it).first) + "\">\n";
						assay_xml +="\t\t<RawFilesGroup_refs>ir_" + String((*it).second) + "</RawFilesGroup_refs>\n";
						assay_xml +="\t\t<Label>\n";
						switch(experiment_type)
						{
								case 0:
										// TODO read from userparams or so which aminoacids were used in light and heavy - for the time being the most common lysine
										if (String(file_descriptions[ assayid_to_filedescriptionsindex[(*it).first] ].label).compare(String("light"))==0)
										{
												assay_xml += "\t\t\t<Modification massDelta=\"0\" residues=\"L\">\n";
												assay_xml += "\t\t\t\t<cvParam cvRef="" accession="" name=\"None\" value=\"0\"/>\n";
										}
										else
										{
												assay_xml += "\t\t\t<Modification massDelta=\"8.0141988132\" residues=\"L\">\n";
												assay_xml += "\t\t\t\t<cvParam name=\"Lys8\" value=\"8.0141988132\"/>\n";

										}
								break;
								case 1:
												//~ TODO extract value from userparam
												assay_xml += "\t\t\t<Modification massDelta=\"290\" residues=\"L\">\n";
												assay_xml += "\t\t\t\t<cvParam name =\"" +file_descriptions[ assayid_to_filedescriptionsindex[(*it).first] ].label + "\"/>\n";
								break;
								//~ case 2: has no label!
						}
						assay_xml += "\t\t\t</Modification>\n";
						assay_xml +="\t\t</Label>\n";
						assay_xml += "\t</Assay>\n";

						// for SILACAnalyzer/iTRAQAnalyzer one assay is one studyvariable, this may change!!!
						study_xml += "\t<StudyVariable id=\"v_" + String(assayid_to_studyvarid[(*it).first]) + "\" name=\"String\">\n";
						study_xml += "\t\t\t<AssayRef assay_ref=\"a_" + String((*it).first) + "\"/>\n";
						study_xml += "\t</StudyVariable>\n";

				}
				assay_xml += "</AssayList>\n";
				study_xml += "</StudyVariableList>\n";
				os << assay_xml << study_xml;

				//--------------------------------------------------------------------------------------------
				// Features and layers
				//--------------------------------------------------------------------------------------------
				std::vector<UInt64> fid;
				std::vector<Real> fin, fwi /*, fqu */;
				std::vector< std::vector<UInt64> > cid;
				std::vector< std::vector<Real> > f2i;
				String feature_xml = "";
				feature_xml += "\t<FeatureList>\n"; // TODO rawfilegroup_refs
				for (ConsensusMap::const_iterator it = ccm_->begin(); it != ccm_->end(); ++it)
				{
						const std::set< FeatureHandle,FeatureHandle::IndexLess>& feature_handles = it->getFeatures();
						switch (experiment_type)
						{
								case 0:
								{
										std::vector<UInt64> cidvec;
										cidvec.push_back(UniqueIdGenerator::getUniqueId());
										for (std::set< FeatureHandle,FeatureHandle::IndexLess>::const_iterator jt = feature_handles.begin(); jt != feature_handles.end(); ++jt)
										{
												fid.push_back(UniqueIdGenerator::getUniqueId());
												cidvec.push_back(fid.back());
												fin.push_back(jt->getIntensity());
												fwi.push_back(jt->getWidth());
												//~ fqu.push_back(jt->getQuality());
												feature_xml += "\t\t<Feature id=\"f_" + String(fid.back()) + "\" RT=\"" + String(jt->getRT()) + "\" MZ=\"" + String(jt->getMZ()) + "\" charge=\"" + String(jt->getCharge()) + "\">\n\t\t</Feature>\n";
												// TODO as soon as SILACanalyzer incorporate convex hulls read from the featuremap
												//~ writeUserParam_(os, *jt, UInt(2)); // FeatureHandle has no MetaInfoInterface!!!
												// TODO as soon as SILACanalyzer and so on can incorporate ids
										}
										cid.push_back(cidvec);
								}break;
								case 1:
								{
										std::vector<Real> fi;
										fid.push_back(UniqueIdGenerator::getUniqueId());
										feature_xml += "\t\t<Feature id=\"f_" + String(fid.back()) + "\" RT=\"" + String(it->getRT()) + "\" MZ=\"" + String(it->getMZ()) + "\" charge=\"" + String(it->getCharge()) + "\"/>\n";
										//~ std::vector<UInt64> cidvec;
										//~ cidvec.push_back(fid.back());
										for (std::set< FeatureHandle,FeatureHandle::IndexLess>::const_iterator jt = feature_handles.begin(); jt != feature_handles.end(); ++jt)
										{
												fi.push_back(jt->getIntensity());
										}
										f2i.push_back(fi);
								}break;
						}
				}

				os << feature_xml;
				switch (experiment_type)
				{
						case 0:
								os << "\t\t<FeatureQuantLayer id=\"" << "q_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n\t\t\t<ColumnIndex>\n";

								os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"intensity\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>";
								os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"width\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>";
								//~ os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"quality\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>"; // getQuality erst ab BaseFeature - nicht in FeatureHandle
								os << "</ColumnIndex>\t\t\t\t\n<DataMatrix>\n";
								for (Size i=0; i < fid.size(); ++i)
								{
										os <<"\t\t\t\t\t<Row object_ref=\"f_" << String(fid[i]) << "\">";
										os << fin[i] << " " << fwi[i] /* << " " << fiq[i] */;
										os << "</Row>\n";
								}
								os << "\t\t\t</DataMatrix>\n";
								os << "\t\t</FeatureQuantLayer>\n";
						break;
						case 1:
								os << "\t\t<MS2AssayQuantLayer id=\""+ String(UniqueIdGenerator::getUniqueId()) +"\">\n\t\t\t<DataType>\n\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"intensity\"/>\n\t\t\t</DataType>\n\t\t\t<Column index>\n\t\t\t\t\t";
								for (Size i = 0; i < file_descriptions.size(); ++i)
								{
										os << filedescriptionsindex_to_assayid[i];
								}
								os <<  "\n\t\t\t\t</ColumnIndex>\t\t\t<DataMatrix>\n";
								for (Size i = 0; i < fid.size(); ++i)
								{
										os << "\t\t\t\t\t<Row object_ref=\"f_" + String(fid[i]) + "\">";
										for (Size j = 0; j < f2i[i].size(); ++j)
										{
												os << String(f2i[i][j]) << " ";
										}
										os << "</Row>\n";
								}
								os << "\t\t\t</DataMatrix>\n\t\t</MS2StudyVariableQuantLayer>\n";
						break;
				}
				os << "\t</FeatureList>\n";


				//--------------------------------------------------------------------------------------------
				// Peptides
				//--------------------------------------------------------------------------------------------
				os << "\t<PeptideList  id=\"" << "m_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n";
				switch (experiment_type)
				{
						case 0: //iterate consensusmap
						{
								for (Size i = 0; i < cid.size(); ++i)
								{
										if (!(*ccm_)[i].getPeptideIdentifications().empty())
										{
												os << "\t\t<PeptideConsensus id=\"" << "c_" << String(cid[i].front()) << "\" charge=\""+ String((*ccm_)[i].getCharge()) +"\">\n";
												os << "\t\t\t<feature_refs>";
												for (Size j=1; j < cid[i].size(); ++j)
												{
														os << "f_" << cid[i][j]<< " ";
												}
												os << "\t\t\t</feature_refs>\n";
												os << "\t\t\t<IdentificationRef id_ref=\"";
												os << (*ccm_)[i].getPeptideIdentifications().front().getIdentifier() << "\" IdentificationFile_ref=\"";
												os << idid_to_idfilenames.begin()->first  << "\"/>\n";
												os << "\t\t</PeptideConsensus>\n";
										}
								}
								// QuantLayers
								os << "\t\t<RatioQuantLayer id=\"" << "q_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n\t\t\t<ColumnIndex>\n";

								os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"ratio\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>";
								//~ os << "\t\t\t\t<Column index=\"0\">\n\t\t\t\t\t<DataType>\n\t\t\t\t\t\t<cvParam cvRef=\"PSI-MS\" accession=\"TODO\" name=\"quality\"/>\n\t\t\t\t\t</DataType>\n\t\t\t\t</Column>"; // getQuality erst ab BaseFeature - nicht in FeatureHandle
								os << "</ColumnIndex>\t\t\t\t\n<DataMatrix>\n";
								for (Size i=0; i < cid.size(); ++i)
								{
										os <<"\t\t\t\t\t<Row object_ref=\"c_" << String(cid[i].front()) << "\">";
										os << String((*ccm_)[i].getIntensity());
										os << "</Row>\n";
								}
								os << "\t\t\t</DataMatrix>\n";
								os << "\t\t</RatioQuantLayer>\n";
						}
						case 1: //iterate featuremap
						{
								for (Size i = 0; i < fid.size(); ++i)
								{
										if (!(*ccm_)[i].getPeptideIdentifications().empty())
										{
												os << "\t\t<PeptideConsensus id=\"" << "c_" << String(fid[i]) << "\" charge=\""+ String((*ccm_)[i].getCharge()) +"\">\n";
												os << "\t\t\t<feature_refs>";
												for (Size j = 0; j < f2i[i].size(); ++j)
												{
														os << String(f2i[i][j]) << " ";
												}
												os << "\t\t\t</feature_refs>\n";
												os << "\t\t\t<IdentificationRef id_ref=\"";
												os << (*ccm_)[i].getPeptideIdentifications().front().getIdentifier() << "\" IdentificationFile_ref=\"";
												os << idid_to_idfilenames.begin()->first  << "\"/>\n";
												os << "\t\t</PeptideConsensus>\n";
										}
								}
								//~ TODO ratios, when available (not yet for the iTRAQ tuples of iTRAQAnalyzer)
						}
				}



				os << "\t</PeptideList>\n";


				//--------------------------------------------------------------------------------------------
				// Ratio
				//--------------------------------------------------------------------------------------------
				if (experiment_type == 0)
				{
						os << "\t<RatioList>\n";
						for (Size i = 0; i < assayid_to_studyvarid.size(); ++i)
						{
								os << "\t\t<Ratio id=\"" << "r_" << String(UniqueIdGenerator::getUniqueId()) << "\">\n";
								//todo paarungen
								os << "\t\t</Ratio>\n";
						}
						os << "\t</RatioList>\n";
				}

				//--------------------------------------------------------------------------------------------
				// Proteins and Proteingroups
				//--------------------------------------------------------------------------------------------
				// TODO - omitted as there are no ids yet

				os << "\t</mzQuantML>\n";
		}

		void MzQuantMLHandler::writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent)
		{
			String h;
			writeUserParam_(h, meta, indent);
			os << h;
		}

		void MzQuantMLHandler::writeUserParam_(String& s, const MetaInfoInterface& meta, UInt indent)
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
