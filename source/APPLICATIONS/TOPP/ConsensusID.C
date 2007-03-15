// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/FORMAT/FeatureMapFile.h>
#include <OpenMS/ANALYSIS/ID/IDFeatureMapper.h>
#include <OpenMS/ANALYSIS/ID/ConsensusID.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page ConsensusID ConsensusID
	
	@brief Combines results of Identification engines.

	@todo write test (Marc)
	@todo document (Marc)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

//Helper class
struct IDData
{
	DoubleReal mz;
	DoubleReal rt;
	vector<Identification> ids;
};

class TOPPConsensusID
	: public TOPPBase
{
	public:
		TOPPConsensusID()
			: TOPPBase("ConsensusID","Combines results of Identification engines.")
		{
			
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerStringOption_("ids","<file>","","one or more analysisXML files separated by comma (without blanks)");
			registerStringOption_("out","<file>","","output file in AnalysisXML format");
			registerStringOption_("features","<file>","","input feature file. If this file is given, all identifications\n"
																									 "are mapped to features and the consensus is made for feaatures.",false);
			registerStringOption_("features_out","<file>","","Features that have identifications are stored in this file."
			                                                 "Only available when 'features' file is given!",false);
			registerSubsection_("algorithm");
		}
	
		ExitCodes main_(int , char**)
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
			outputFileWritable_(out);

			String feature_file = getStringOption_("features");
			String feature_out_file = "";
			bool feature_mode = false;
			if (feature_file!="")
			{
				inputFileReadable_(feature_file);
				feature_mode = true;
				feature_out_file = getStringOption_("features_out");
				if (feature_out_file!="")
				{
					outputFileWritable_(feature_out_file);
				}
			}
			
			//-------------------------------------------------------------
			// loading input + calculations
			//-------------------------------------------------------------
			
			//set up ConsensusID
			ConsensusID consensus;
			const Param& alg_param = getParam_().copy("algorithm:",true);
			writeDebug_("Parameters passed to ConsensusID", alg_param, 3);
			if (alg_param.empty())
			{
				writeLog_("No parameters for ConsensusID given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			consensus.setParameters(alg_param);

			AnalysisXMLFile ax_file;
			vector<ProteinIdentification> prot_ids;
			vector<IdentificationData> all_ids;

			//feature mode
			if (feature_mode)
			{
				//load features
				FeatureMap<> features;
				FeatureMapFile feat_file;
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
					cout << "ConsensusID -- Feature " << features[i].getRT() << " / " << features[i].getMZ() << endl;
					consensus.apply(features[i].getIdentifications());
				}
				
				// writing output
				all_ids.clear();
				IdentificationData id;
				for (UInt i = 0; i < features.size(); ++i)
				{
					id.rt = features[i].getRT();
					id.mz = features[i].getMZ();
					for (vector<Identification>::const_iterator id_it = features[i].getIdentifications().begin(); 
							 id_it != features[i].getIdentifications().end(); 
							 ++id_it)
					{
						id.id = *id_it;
						all_ids.push_back(id);
					}
				}
				ax_file.store(out,features.getProteinIdentifications(),all_ids);
				
				//write output features (those with IDs)
				if (feature_out_file!="")
				{
					FeatureMap<> features_out;
					for (UInt i = 0; i < features.size(); ++i)
					{
						if (features[i].getIdentifications().size()!=0 && features[i].getIdentifications()[0].getPeptideHits().size()!=0)
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
				//Storage for protein identifications
				vector< ProteinIdentification > prot_id_out;
				
				//load and merge ids (by precursor position)
				for(UInt i = 0; i < ids.size(); ++i)
				{
					writeDebug_(String("Mapping ids: ") + ids[i], 2);
					ax_file.load(ids[i],prot_ids, all_ids);
					// Append protein IDs
					prot_id_out.insert(prot_id_out.end(), prot_ids.begin(), prot_ids.end());
					// Insert peptide IDs
					for (vector<IdentificationData>::iterator ins = all_ids.begin(); ins != all_ids.end(); ++ins)
					{
						vector<IDData>::iterator pos = prec_data.begin();
						while (pos != prec_data.end())
						{
							if (fabs(pos->rt - ins->rt) < 0.1 && fabs(pos->mz - ins->mz) < 0.1 )
							{
								break;
							}
						}
						//right position was found => append ids
						if (pos != prec_data.end())
						{
							pos->ids.push_back(ins->id);
						}
						//insert new entry
						else
						{
							IDData tmp;
							tmp.mz = ins->mz;
							tmp.rt = ins->rt;
							tmp.ids.push_back(ins->id);
							prec_data.push_back(tmp);
						}
					}
				}
				
				//do consensus
				for (vector<IDData>::iterator it = prec_data.begin(); it!=prec_data.end(); ++it)
				{
					cout << "ConsensusID -- Precursor " << it->rt << " / " << it->mz  << endl;
					consensus.apply(it->ids);
				}
				
				// writing output
				all_ids.clear();
				IdentificationData id;
				for (vector<IDData>::iterator it = prec_data.begin(); it!=prec_data.end(); ++it)
				{
					id.rt = it->rt;
					id.mz = it->mz;
					id.id = it->ids[0];
					all_ids.push_back(id);
				}
				ax_file.store(out,prot_id_out,all_ids);
			}
			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPConsensusID tool;
	return tool.main(argc,argv);
}

/// @endcond
