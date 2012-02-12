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
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen and others$
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/ANALYSIS/ID/ConsensusID.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ConsensusID ConsensusID

	@brief Computes a consensus identification from peptide identification engines.

	<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ ConsensusID \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_FalseDiscoveryRate </td>
		</tr>
		<tr>
		  <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDPosteriorErrorProbability </td>
		</tr>
		<tr>
		  <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDMapper </td>
		</tr>
	</table>
	</CENTER>

	The input file can contain several searches, e.g. from several identification engines. In order
	to use the PEPMatrix or the PEPIons algorithm, posterior
	error probabilities (PEPs) need to be calculated using the @ref TOPP_IDPosteriorErrorProbability tool
	for all individual search engines. After PEP calculation, the different search engine results
	have to be combined using @ref TOPP_IDMerger. Identification runs can be mapped
	to featureXML and consensusXML with the @ref TOPP_IDMapper tool. The merged file can now be fed into
	into the @ref TOPP_ConsensusID tool. For the statistical assessment of the results it is recommended
	to use target-decoy databases for peptide identifications. The false discovery rates (FDRs) can be
	calcultated using the @ref TOPP_FalseDiscoveryRate tool.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ConsensusID.cli

	For the parameters of the algorithm section see the algorithms documentation: @n
		@ref OpenMS::ConsensusID "Consensus algorithm" @n
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

//Helper class
struct IDData
{
	DoubleReal mz;
	DoubleReal rt;
	String sourcefile;
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
			setValidFormats_("in",StringList::create("idXML,featureXML,consensusXML"));
			registerOutputFile_("out","<file>","","output file");
			setValidFormats_("out",StringList::create("idXML,featureXML,consensusXML"));

			addEmptyLine_();
			registerDoubleOption_("rt_delta","<value>",0.1, "Maximum allowed precursor RT deviation between identifications.", false);
			setMinFloat_("rt_delta",0.0);
			registerDoubleOption_("mz_delta","<value>",0.1, "Maximum allowed precursor m/z deviation between identifications.", false);
			setMinFloat_("mz_delta",0.0);
			registerIntOption_("min_length","<value>",6, "Minimum of length of peptides for final consensus list", false);
			setMinInt_("min_length",1);
			registerFlag_("use_all_hits","If 'true' not only the first hit, but all are used (peptides only)");

			registerSubsection_("algorithm","Consensus algorithm section");
		}

		ExitCodes main_(int , const char**)
		{
			String in = getStringOption_("in");
			FileTypes::Type in_type = FileHandler::getType(in);
			String out = getStringOption_("out");
			bool use_all_hits(getFlag_("use_all_hits"));

			DoubleReal rt_delta = getDoubleOption_("rt_delta");
			DoubleReal mz_delta = getDoubleOption_("mz_delta");
			UInt min_length = getIntOption_("min_length");

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
			if (in_type == FileTypes::IDXML)
			{
				vector<ProteinIdentification> prot_ids;
				vector<PeptideIdentification> pep_ids;
				String document_id;
				IdXMLFile().load(in,prot_ids, pep_ids, document_id);

				//merge peptide ids by precursor position Sven:Ideally one should merge all peptide hits form the the different peptide identifications and keep the the information on the identification runs as a meta value
				vector<IDData> prec_data,final;
				for (vector<PeptideIdentification>::iterator pep_id_it = pep_ids.begin(); pep_id_it != pep_ids.end(); ++pep_id_it)
				{
					PeptideIdentification t;
					t=*pep_id_it;
					String file_origin = (String)pep_id_it->getMetaValue("file_origin");
					String scoring = (String)pep_id_it->getIdentifier();
					DoubleReal rt = (DoubleReal)(pep_id_it->getMetaValue("RT"));
					DoubleReal mz = (DoubleReal)(pep_id_it->getMetaValue("MZ"));
					writeDebug_(String("  ID: ") + rt + " / " + mz, 4);
					vector<IDData>::iterator pos = prec_data.begin();
					while (pos != prec_data.end())
					{
						if (fabs(pos->rt - rt) < rt_delta && fabs(pos->mz - mz) < mz_delta && pos->sourcefile == file_origin)
						{
							break;
						}
						++pos;
					}
					//right position was found => append ids
					if (pos != prec_data.end())
					{
						writeDebug_(String("    Appending IDs to precursor: ") + pos->rt + " / " + pos->mz, 4);
						//write information on search engine
						vector<PeptideHit> hits;
						for (vector<PeptideHit>::const_iterator pit = t.getHits().begin(); pit != t.getHits().end();++pit)
							{
								PeptideHit hit = *pit;
								if (hit.getSequence().size()>=min_length)
								{
									if (hit.metaValueExists("scoring"))
									{
										String meta_value = (String)hit.getMetaValue("scoring");
									}
									hit.setMetaValue("scoring", pep_id_it->getIdentifier());
									hits.push_back(hit);
									if (!use_all_hits || pit->getScore() > 0.98)
									{
										break;
									}
								}
								}
						t.setHits(hits);
						pos->sourcefile=file_origin;
						pos->ids.push_back(t);
					}
					//insert new entry
					else
					{
						IDData tmp;
						tmp.mz = mz;
						tmp.rt = rt;
						tmp.sourcefile = file_origin;
						vector<PeptideHit> hits;
						for (vector<PeptideHit>::const_iterator pit = t.getHits().begin(); pit != t.getHits().end();++pit)
						{
							PeptideHit hit = *pit;
							if (hit.getSequence().size()>=min_length)
							{
								if (hit.metaValueExists("scoring"))
								{
									String meta_value = (String)hit.getMetaValue("scoring");
								}
								hit.setMetaValue("scoring", pep_id_it->getIdentifier());
								hits.push_back(hit);
								if (!use_all_hits || pit->getScore() > 0.98)
								{
									break;
								}
							}
							cout<<pep_id_it->getIdentifier()<<endl;
     				}
						t.setHits(hits);
						tmp.ids.push_back(t);
						prec_data.push_back(tmp);
						writeDebug_(String("    Inserting new precursor: ") + tmp.rt + " / " + tmp.mz, 4);
					}
				}
				//iterate over prec_data and write to final only one peptide identification per rt mz

				for (vector<IDData>::iterator fin = prec_data.begin(); fin != prec_data.end(); ++fin)
				{
					IDData tmp;
					tmp.mz=fin->mz;
					tmp.rt=fin->rt;
					tmp.sourcefile = fin->sourcefile;
					PeptideIdentification t;
					vector<PeptideHit> P;

					for (vector<PeptideIdentification>::iterator tt = fin->ids.begin(); tt != fin->ids.end(); ++tt)
					{
						for (vector<PeptideHit>::const_iterator pit = tt->getHits().begin(); pit != tt->getHits().end();++pit)
							{
								P.push_back(*pit);
							}
					}
					t.setHits(P);
					tmp.ids.push_back(t);
					final.push_back(tmp);
				}


				///iterate over prec_data and write to final only one peptide identification per rt mz

				//compute consensus
				alg_param.setValue("number_of_runs",(UInt)prot_ids.size());
				consensus.setParameters(alg_param);
				for (vector<IDData>::iterator it = final.begin(); it!=final.end(); ++it)
				{
					writeDebug_(String("Calculating consensus for : ") + it->rt + " / " + it->mz + " #peptide ids: " + it->ids.size(), 4);
					consensus.apply(it->ids);
				}

				// writing output
				pep_ids.clear();
				for (vector<IDData>::iterator it = final.begin(); it!=final.end(); ++it)
				{
					pep_ids.push_back(it->ids[0]);
					pep_ids.back().setMetaValue("RT",it->rt);
					pep_ids.back().setMetaValue("MZ",it->mz);
					pep_ids.back().setMetaValue("file_origin",it->sourcefile);
				}

				//create new identification run
				vector< ProteinIdentification > prot_id_out(1);
				prot_id_out[0].setDateTime(DateTime::now());
				prot_id_out[0].setSearchEngine("OpenMS/ConsensusID");
				prot_id_out[0].setSearchEngineVersion(VersionInfo::getVersion());

				//store consensus
				IdXMLFile().store(out,prot_id_out,pep_ids);
			}

			//----------------------------------------------------------------
			// featureXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::FEATUREXML)
			{
				//load map
				FeatureMap<> map;
				FeatureXMLFile().load(in,map);

				//compute consensus
				alg_param.setValue("number_of_runs",(UInt)map.getProteinIdentifications().size());
				consensus.setParameters(alg_param);
				for (Size i = 0; i < map.size(); ++i)
				{
					consensus.apply(map[i].getPeptideIdentifications());
				}

				//create new identification run
				map.getProteinIdentifications().clear();
				map.getProteinIdentifications().resize(1);
				map.getProteinIdentifications()[0].setDateTime(DateTime::now());
				map.getProteinIdentifications()[0].setSearchEngine("OpenMS/ConsensusID");
				map.getProteinIdentifications()[0].setSearchEngineVersion(VersionInfo::getVersion());

				//store consensus
				FeatureXMLFile().store(out,map);
			}

			//----------------------------------------------------------------
			// consensusXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::CONSENSUSXML)
			{
				//load map
				ConsensusMap map;
				ConsensusXMLFile().load(in,map);

				//compute consensus
				alg_param.setValue("number_of_runs",(UInt)map.getProteinIdentifications().size());
				consensus.setParameters(alg_param);
				for (Size i = 0; i < map.size(); ++i)
				{
					consensus.apply(map[i].getPeptideIdentifications());
				}

				//create new identification run
				map.getProteinIdentifications().clear();
				map.getProteinIdentifications().resize(1);
				map.getProteinIdentifications()[0].setDateTime(DateTime::now());
				map.getProteinIdentifications()[0].setSearchEngine("OpenMS/ConsensusID");
				map.getProteinIdentifications()[0].setSearchEngineVersion(VersionInfo::getVersion());

				//store consensus
				ConsensusXMLFile().store(out,map);
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
