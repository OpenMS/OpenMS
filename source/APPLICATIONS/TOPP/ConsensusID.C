// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDFeatureMapper.h>
#include <OpenMS/ANALYSIS/ID/ConsensusID.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page ConsensusID_TOPP ConsensusID
	
	@brief Computes a consensus identification from peptide identification engines.
	
	For a detailed description of the algorithms and parameters see the documentation of
	the ConsensusID class.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

//Helper class
struct IDData
{
	DoubleReal mz;
	DoubleReal rt;
	vector<PeptideIdentification> ids;
};

class TOPPConsensusID
	: public TOPPBase
{
	public:
		TOPPConsensusID()
			: TOPPBase("ConsensusID","Computes a consensus identification from peptide identifications of several identification engines.")
		{
			
		}
	
	protected:

		Param getSubsectionDefaults_(const String& /*section*/) const
		{
			return ConsensusID().getDefaults();
		}

		void registerOptionsAndFlags_()
		{
			registerStringOption_("ids","<file>","","one or more IdXML files separated by comma (without blanks)");
			registerOutputFile_("out","<file>","","output file ");
			setValidFormats_("out",StringList::create("IdXML"));
			registerInputFile_("features","<file>","","feature input file. If this file is given, all identifications\n"
																								"are mapped to features and the consensus is made for features.\n",false);
			setValidFormats_("features",StringList::create("FeatureXML"));
			registerOutputFile_("features_out","<file>","","Features that have identifications are stored in this file.\n"
			                                               "This file is created only if the 'features' file is given.\n",false);
			setValidFormats_("features_out",StringList::create("FeatureXML"));
			registerSubsection_("algorithm","Consensus algorithm section");
		}
	
		ExitCodes main_(int , const char**)
		{

			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
			
			vector<String> ids;
			String tmp = getStringOption_("ids");
			tmp.split(',', ids);
			if (ids.size() == 0)
			{
				ids.push_back(tmp);
			}
			for(UInt i = 0; i < ids.size(); ++i)
			{
				inputFileReadable_(ids[i]);
			}
	
			String out = getStringOption_("out");

			String feature_file = getStringOption_("features");
			String feature_out_file = "";
			bool feature_mode = false;
			if (feature_file!="")
			{
				feature_mode = true;
				feature_out_file = getStringOption_("features_out");
			}
			
			//-------------------------------------------------------------
			// loading input + calculations
			//-------------------------------------------------------------
			
			//set up ConsensusID
			ConsensusID consensus;
			Param alg_param = getParam_().copy("algorithm:",true);
			if (alg_param.empty())
			{
				writeLog_("No parameters for ConsensusID given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			alg_param.setValue("number_of_runs",(UInt)ids.size());
			writeDebug_("Parameters passed to ConsensusID", alg_param, 3);
			consensus.setParameters(alg_param);

			IdXMLFile ax_file;
			vector<ProteinIdentification> prot_ids;
			vector<PeptideIdentification> all_ids;

			//feature mode
			if (feature_mode)
			{
				//load features
				FeatureMap<> features;
				FeatureXMLFile feat_file;
				feat_file.load(feature_file,features);
				
				//map ids to features
				IDFeatureMapper mapper;
				for(UInt i = 0; i < ids.size(); ++i)
				{
					writeDebug_(String("Mapping ids: ") + ids[i], 2);
					ax_file.load(ids[i],prot_ids, all_ids);
					mapper.annotate(features,all_ids,prot_ids);
				}
				
				//do consensus
				for (UInt i = 0; i < features.size(); ++i)
				{
					//cout << "ConsensusID -- Feature " << features[i].getRT() << " / " << features[i].getMZ() << " ("<< i+1 << ")" << endl;
					consensus.apply(features[i].getPeptideIdentifications());
				}
				
				// writing output
				all_ids.clear();
				for (UInt i = 0; i < features.size(); ++i)
				{

					for (vector<PeptideIdentification>::const_iterator id_it = features[i].getPeptideIdentifications().begin(); 
							 id_it != features[i].getPeptideIdentifications().end(); 
							 ++id_it)
					{
						all_ids.push_back(*id_it);
						all_ids.back().setMetaValue("RT",features[i].getRT());
						all_ids.back().setMetaValue("MZ",features[i].getMZ());
					}
				}
				vector< ProteinIdentification > prot_id_out(1);
				DateTime date;
				date.now();
				prot_id_out[0].setDateTime(date);
				prot_id_out[0].setSearchEngine("OpenMS/ConsensusID");
				prot_id_out[0].setSearchEngineVersion(VersionInfo::getVersion());
				ax_file.store(out,prot_id_out,all_ids);
				
				//write output features (those with IDs)
				if (feature_out_file!="")
				{
					FeatureMap<> features_out;
					for (UInt i = 0; i < features.size(); ++i)
					{
						if (features[i].getPeptideIdentifications().size()!=0 && features[i].getPeptideIdentifications()[0].getHits().size()!=0)
						{
							features_out.push_back(features[i]);
						}
					}
					feat_file.store(feature_out_file,features_out);
				}
			}
			else //non-feature mode
			{
				//Identifications merged by precursor position
				vector<IDData> prec_data;
				
				//load and merge ids (by precursor position)
				for(UInt i = 0; i < ids.size(); ++i)
				{
					writeDebug_(String("Mapping ids: ") + ids[i], 2);
					ax_file.load(ids[i],prot_ids, all_ids);
					// Insert peptide IDs
					for (vector<PeptideIdentification>::iterator ins = all_ids.begin(); ins != all_ids.end(); ++ins)
					{
						DoubleReal rt = (DoubleReal)(ins->getMetaValue("RT"));
						DoubleReal mz = (DoubleReal)(ins->getMetaValue("MZ"));
						writeDebug_(String("  ID: ") + rt + " / " + mz, 4);
						vector<IDData>::iterator pos = prec_data.begin();
						while (pos != prec_data.end())
						{
							if (fabs(pos->rt - rt) < 0.1 && fabs(pos->mz - mz) < 0.1 )
							{
								break;
							}
							++pos;
						}
						//right position was found => append ids
						if (pos != prec_data.end())
						{
							writeDebug_(String("    Appending IDs to precursor: ") + pos->rt + " / " + pos->mz, 4);
							pos->ids.push_back(*ins);
						}
						//insert new entry
						else
						{
							IDData tmp;
							tmp.mz = mz;
							tmp.rt = rt;
							tmp.ids.push_back(*ins);
							prec_data.push_back(tmp);
							writeDebug_(String("    Inserting new precursor: ") + tmp.rt + " / " + tmp.mz, 4);
						}
					}
				}
				
				//do consensus
				for (vector<IDData>::iterator it = prec_data.begin(); it!=prec_data.end(); ++it)
				{
					writeDebug_(String("Calculating consensus for : ") + it->rt + " / " + it->mz + " #peptide ids: " + it->ids.size(), 4);
					//cout << "ConsensusID -- Precursor " << it->rt << " / " << it->mz  << endl;
					/*
					for (UInt i = 0; i != it->ids.size(); ++i)
					{
						vector<PeptideHit> hits = it->ids[i].getHits();
						for (vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2)
						{
							String seq(it2->getSequence());
							replace(seq.begin(), seq.end(), 'I', 'L');
							it2->setSequence(seq);
						}
						it->ids[i].setHits(hits);
					}
					*/
					consensus.apply(it->ids);
				}
				
				// writing output
				all_ids.clear();
				PeptideIdentification id;
				DateTime date;
				date.now();
				String date_str;
				date.get(date_str);
				for (vector<IDData>::iterator it = prec_data.begin(); it!=prec_data.end(); ++it)
				{
					all_ids.push_back(it->ids[0]);
					all_ids.back().setMetaValue("RT",it->rt);
					all_ids.back().setMetaValue("MZ",it->mz);
					all_ids.back().setIdentifier("Consensus_" + date_str);
				}
				
				//store consensus
				vector< ProteinIdentification > prot_id_out(1);
				prot_id_out[0].setDateTime(date);
				prot_id_out[0].setSearchEngine("OpenMS/ConsensusID");
				prot_id_out[0].setSearchEngineVersion(VersionInfo::getVersion());
				prot_id_out[0].setIdentifier("Consensus_" + date_str);
				ax_file.store(out,prot_id_out,all_ids);
			}
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPConsensusID tool;
	return tool.main(argc,argv);
}

/// @endcond
