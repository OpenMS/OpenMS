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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/ANALYSIS/ID/ConsensusID.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page ConsensusID_TOPP ConsensusID
	
	@brief Computes a consensus identification from peptide identification engines.
	
	The input file can contain several searches, e.g. from several identification engines.
	
	@see IDMerger
	
	For a detailed description of the algorithms and parameters see the documentation of
	the %OpenMS ConsensusID class.
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
			registerInputFile_("in","<file>","","input file");
			setValidFormats_("in",StringList::create("IdXML,featureXML,consensusXML"));
			registerOutputFile_("out","<file>","","output file");
			setValidFormats_("out",StringList::create("IdXML,featureXML,consensusXML"));

			addEmptyLine_();
			registerDoubleOption_("rt_delta","<value>",0.1, "Maximum allowed precursor RT deviation between identifications.", false);
			setMinFloat_("rt_delta",0.0);
			registerDoubleOption_("mz_delta","<value>",0.1, "Maximum allowed precursor m/z deviation between identifications.", false);
			setMinFloat_("mz_delta",0.0);
			
			registerSubsection_("algorithm","Consensus algorithm section");
		}
	
		ExitCodes main_(int , const char**)
		{
			String in = getStringOption_("in");
			FileHandler::Type in_type = FileHandler::getType(in);
			String out = getStringOption_("out");
			
			DoubleReal rt_delta = getDoubleOption_("rt_delta");
			DoubleReal mz_delta = getDoubleOption_("mz_delta");
			
			//----------------------------------------------------------------
			//set up ConsensusID
			//----------------------------------------------------------------
			ConsensusID consensus;
			Param alg_param = getParam_().copy("algorithm:",true);
			if (alg_param.empty())
			{
				writeLog_("No parameters for ConsensusID given. Aborting!");
				return ILLEGAL_PARAMETERS;
			}
			writeDebug_("Parameters passed to ConsensusID (without number of runs)", alg_param, 3);

			//----------------------------------------------------------------
			// idXML
			//----------------------------------------------------------------
			if (in_type == FileHandler::IDXML)
			{
				vector<ProteinIdentification> prot_ids;
				vector<PeptideIdentification> pep_ids;
				IdXMLFile().load(in,prot_ids, pep_ids);
			
				//merge peptide ids by precursor position
				vector<IDData> prec_data;
				for (vector<PeptideIdentification>::iterator pep_id_it = pep_ids.begin(); pep_id_it != pep_ids.end(); ++pep_id_it)
				{
					DoubleReal rt = (DoubleReal)(pep_id_it->getMetaValue("RT"));
					DoubleReal mz = (DoubleReal)(pep_id_it->getMetaValue("MZ"));
					writeDebug_(String("  ID: ") + rt + " / " + mz, 4);
					vector<IDData>::iterator pos = prec_data.begin();
					while (pos != prec_data.end())
					{
						if (fabs(pos->rt - rt) < rt_delta && fabs(pos->mz - mz) < mz_delta )
						{
							break;
						}
						++pos;
					}
					//right position was found => append ids
					if (pos != prec_data.end())
					{
						writeDebug_(String("    Appending IDs to precursor: ") + pos->rt + " / " + pos->mz, 4);
						pos->ids.push_back(*pep_id_it);
					}
					//insert new entry
					else
					{
						IDData tmp;
						tmp.mz = mz;
						tmp.rt = rt;
						tmp.ids.push_back(*pep_id_it);
						prec_data.push_back(tmp);
						writeDebug_(String("    Inserting new precursor: ") + tmp.rt + " / " + tmp.mz, 4);
					}
				}
				
				//compute consensus
				alg_param.setValue("number_of_runs",(UInt)prot_ids.size());
				consensus.setParameters(alg_param);
				for (vector<IDData>::iterator it = prec_data.begin(); it!=prec_data.end(); ++it)
				{
					writeDebug_(String("Calculating consensus for : ") + it->rt + " / " + it->mz + " #peptide ids: " + it->ids.size(), 4);
					consensus.apply(it->ids);
				}
				
				// writing output
				pep_ids.clear();
				DateTime date;
				date.now();
				for (vector<IDData>::iterator it = prec_data.begin(); it!=prec_data.end(); ++it)
				{
					pep_ids.push_back(it->ids[0]);
					pep_ids.back().setMetaValue("RT",it->rt);
					pep_ids.back().setMetaValue("MZ",it->mz);
					pep_ids.back().setIdentifier("Consensus_" + date.get());
				}
				
				//store consensus
				vector< ProteinIdentification > prot_id_out(1);
				prot_id_out[0].setDateTime(date);
				prot_id_out[0].setSearchEngine("OpenMS/ConsensusID");
				prot_id_out[0].setSearchEngineVersion(VersionInfo::getVersion());
				prot_id_out[0].setIdentifier("Consensus_" + date.get());
				IdXMLFile().store(out,prot_id_out,pep_ids);
			}
			
			//----------------------------------------------------------------
			// featureXML
			//----------------------------------------------------------------
//			if (in_type == FileHandler::FEATUREXML)
//			{
//				//load features
//				FeatureMap<> features;
//				FeatureXMLFile feat_file;
//				feat_file.load(feature_file,features);
//				
//				//map ids to features
//				IDMapper mapper;
//				for(UInt i = 0; i < ids.size(); ++i)
//				{
//					writeDebug_(String("Mapping ids: ") + ids[i], 2);
//					IdXMLFile().load(ids[i],prot_ids, all_ids);
//					mapper.annotate(features,all_ids,prot_ids);
//				}
//				
//				//do consensus
//				for (UInt i = 0; i < features.size(); ++i)
//				{
//					//cout << "ConsensusID -- Feature " << features[i].getRT() << " / " << features[i].getMZ() << " ("<< i+1 << ")" << endl;
//					consensus.apply(features[i].getPeptideIdentifications());
//				}
//				
//				// writing output
//				all_ids.clear();
//				for (UInt i = 0; i < features.size(); ++i)
//				{
//
//					for (vector<PeptideIdentification>::const_iterator id_it = features[i].getPeptideIdentifications().begin(); 
//							 id_it != features[i].getPeptideIdentifications().end(); 
//							 ++id_it)
//					{
//						all_ids.push_back(*id_it);
//						all_ids.back().setMetaValue("RT",features[i].getRT());
//						all_ids.back().setMetaValue("MZ",features[i].getMZ());
//					}
//				}
//				vector< ProteinIdentification > prot_id_out(1);
//				DateTime date;
//				date.now();
//				prot_id_out[0].setDateTime(date);
//				prot_id_out[0].setSearchEngine("OpenMS/ConsensusID");
//				prot_id_out[0].setSearchEngineVersion(VersionInfo::getVersion());
//				IdXMLFile().store(out,prot_id_out,all_ids);
//				
//				//write output features (those with IDs)
//				if (feature_out_file!="")
//				{
//					FeatureMap<> features_out;
//					for (UInt i = 0; i < features.size(); ++i)
//					{
//						if (features[i].getPeptideIdentifications().size()!=0 && features[i].getPeptideIdentifications()[0].getHits().size()!=0)
//						{
//							features_out.push_back(features[i]);
//						}
//					}
//					feat_file.store(feature_out_file,features_out);
//				}
//			}

			//----------------------------------------------------------------
			// consensusXML
			//----------------------------------------------------------------
			if (in_type == FileHandler::CONSENSUSXML)
			{
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
