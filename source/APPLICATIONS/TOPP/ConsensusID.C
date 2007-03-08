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
	

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

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
			registerStringOption_("in","<file>","","input feature file");
			registerStringOption_("ids","<file>","","one or more analysisXML files separated by comma (without blanks)");
			registerStringOption_("out","<file>","","output file in AnalysisXML format");
			
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

			String in = getStringOption_("in");
			inputFileReadable_(in);
	
			String out = getStringOption_("out");
			outputFileWritable_(out);
			
			//-------------------------------------------------------------
			// loading input + calculations
			//-------------------------------------------------------------
			
			//load features
			FeatureMap<> features;
			FeatureMapFile().load(in,features);
			
			//map ids to features
			IDFeatureMapper mapper;
			AnalysisXMLFile ax_file;
			vector<ProteinIdentification> prot_ids;
			vector<IdentificationData> all_ids;
			for(UInt i = 0; i < ids.size(); ++i)
			{
				ax_file.load(ids[i],prot_ids, all_ids);
				mapper.annotate(features,all_ids,prot_ids);
			}
			
			//do merge
			ConsensusID consensus;
			const Param& alg_param = getParam_().copy("algorithm:",true);
	
			writeDebug_("Parameters passed to ConsensusID", alg_param, 3);
			
			if (alg_param.empty())
			{
				writeLog_("No parameters for ConsensusID given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			
			consensus.setParameters(alg_param);
			for (UInt i = 0; i < features.size(); ++i)
			{
				consensus.apply(features[i]);
			}
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
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

			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPConsensusID tool;
	return tool.main(argc,argv);
}

/// @endcond
