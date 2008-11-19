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
// $Maintainer: Lars Nilse, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page QuantIDMerger QuantIDMerger
	
	@todo use IDMapper and write output to file (Marc, Chris)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPQuantIDMerger
	: public TOPPBase
{

	public:

		TOPPQuantIDMerger()
			: TOPPBase("QuantIDMerger", "")
		{
		}
		
		struct QuantData
		{
			Int id;
			DoubleReal rt;
			DoubleReal mz;
			
			QuantData()
				: id(-1),
					rt(0.0),
					mz(0.0)
			{
			}
		};
	
	protected:
	
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in_quant", "<file>", "", "Qunatitation input file.");
			setValidFormats_("in_quant",StringList::create("featureXML,consensusXML"));
			registerInputFile_("in_id", "<file>", "", "Identification input file.");
			setValidFormats_("in_id",StringList::create("idXML"));
			registerOutputFile_("out", "<file>", "", "Output file in text format");			
			
			addEmptyLine_();
			registerDoubleOption_("rt_cutoff","<value>",0.0, "Maximum allowed RT deviation between identification and quantitation.");
			setMinFloat_("rt_cutoff",0.0);
			registerDoubleOption_("mz_cutoff","<value>",0.0, "Maximum allowed m/z deviation between identification and quantitation.");
			setMinFloat_("mz_cutoff",0.0);

		}
	
		ExitCodes main_(int , const char**)
		{
			//----------------------------------------------------------------
			//load quant data
			String in_quant = getStringOption_("in_quant");
			FileHandler::Type in_type = FileHandler::getType(in_quant);
			
			vector<QuantData> quant_data;
			ConsensusMap consensus_map;
			FeatureMap<> feature_map;
			if (in_type == FileHandler::CONSENSUSXML)
			{
				ConsensusXMLFile().load(in_quant,consensus_map);
				for (UInt i=0; i<consensus_map.size(); ++i)
				{
					const ConsensusFeature::HandleSetType& handles = consensus_map[i].getFeatures();
					for (ConsensusFeature::HandleSetType::const_iterator it = handles.begin(); it!=handles.end(); ++it)
					{
						QuantData tmp;
						tmp.id = i;
						tmp.mz = it->getMZ();
						tmp.rt = it->getRT();
						quant_data.push_back(tmp);
					}
				}
			}
			else if (in_type == FileHandler::FEATUREXML)
			{
				FeatureXMLFile().load(in_quant,feature_map);
				for (UInt i=0; i<feature_map.size(); ++i)
				{
					QuantData tmp;
					tmp.id = i;
					tmp.mz = feature_map[i].getMZ();
					tmp.rt = feature_map[i].getRT();
					quant_data.push_back(tmp);
				}
			}
			
			//----------------------------------------------------------------
			//load id data
			vector<ProteinIdentification> protein_ids;
			vector<PeptideIdentification> peptide_ids;
			String in_id = getStringOption_("in_id");
			IdXMLFile().load(in_id,protein_ids,peptide_ids);
			
			//----------------------------------------------------------------
			//parse identification engine data
			Map<String, String> engines;
			for (UInt i=0; i<protein_ids.size(); ++i)
			{
				engines[protein_ids[i].getIdentifier()] = protein_ids[i].getSearchEngine()
				                                         + " "
				                                         + protein_ids[i].getSearchEngineVersion()
				                                         + " "
				                                         + protein_ids[i].getDateTime().get();
			}
			
			//----------------------------------------------------------------
			//map id to quant
			DoubleReal mz_cutoff = getDoubleOption_("mz_cutoff");
			DoubleReal rt_cutoff = getDoubleOption_("rt_cutoff");
			Map<UInt,vector<UInt> > quant_to_ids;
			for (UInt i=0; i<peptide_ids.size(); ++i)
			{
				//determine hits
				vector<UInt> matches; 
				for (UInt q=0; q<quant_data.size(); ++q)
				{
					if (fabs((DoubleReal)(peptide_ids[i].getMetaValue("MZ"))-quant_data[q].mz)<=mz_cutoff
						  &&
						  fabs((DoubleReal)(peptide_ids[i].getMetaValue("RT"))-quant_data[q].rt)<=rt_cutoff
						 )
					{
						matches.push_back(q);
					}
				}
				//one hit => report the one hit
				if (matches.size()==1)
				{
					quant_to_ids[matches[0]].push_back(i);
				}
				//several hits => assign to closest data points
				else
				{
					DoubleReal min_dist = 1.0;
					Int min_index = -1;
					for (UInt m=0; m<matches.size(); ++m)
					{
						DoubleReal dist_rt = fabs((DoubleReal)(peptide_ids[i].getMetaValue("RT"))-quant_data[m].rt)/rt_cutoff;
						DoubleReal dist_mz = fabs((DoubleReal)(peptide_ids[i].getMetaValue("MZ"))-quant_data[m].mz)/mz_cutoff;
						DoubleReal dist = sqrt(dist_rt*dist_rt + dist_mz*dist_mz);
						if (dist < min_dist)
						{
							min_dist = dist;
							min_index = matches[m];
						}
					}
					quant_to_ids[min_index].push_back(i);
				}
			}
			//----------------------------------------------------------------
			//output
			for (UInt q=0; q<quant_data.size(); ++q)
			{
				//print quantitation data
				if (in_type == FileHandler::CONSENSUSXML)
				{
					
					// TODO: this is wrong!! ConsensusMap does not have quant_data.size() entries (it has only that many featureHandles!)
					cout << "Consensus feature " << q << ":" << endl;
					cout << "- rt        : " << consensus_map[q].getRT() << endl;
					cout << "- mz        : " << consensus_map[q].getMZ() << endl;
					cout << "- ratio     : " << consensus_map[q].getIntensity() << endl;
					cout << "- charge    : " << consensus_map[q].getCharge() << endl;
				}
				else if (in_type == FileHandler::FEATUREXML)
				{
					cout << "Feature " << q << ":" << endl;
					cout << "- rt        : " << feature_map[q].getRT() << endl;
					cout << "- mz        : " << feature_map[q].getMZ() << endl;
					cout << "- intensity : " << feature_map[q].getIntensity() << endl;
					cout << "- charge    : " << feature_map[q].getCharge() << endl;
				}
				//print hits
				if (quant_to_ids.has(q))
				{
					for (UInt i=0; i<quant_to_ids[q].size();++i)
					{
						const PeptideIdentification& pep_id = peptide_ids[quant_to_ids[q][i]];
						cout << "- Peptide Identification - rt: " << pep_id.getMetaValue("RT") 
						                             << " - mz: " << pep_id.getMetaValue("MZ") 
						                             << " - engine: " << engines[pep_id.getIdentifier()]
						                             << endl;
						for (UInt h=0; h<pep_id.getHits().size(); ++h)
						{
							const PeptideHit& hit = pep_id.getHits()[h];
							cout << "  - hit: " << hit.getSequence() << " (" << hit.getScore() << ")" << endl; 
						}
					}
				}
				cout << endl;
			}
			
			
			//----------------------------------------------------------------
			//write output
			
			return EXECUTION_OK;
		}

};


int main( int argc, const char** argv )
{
	TOPPQuantIDMerger tool;
	return tool.main(argc, argv);
}

/// @endcond
