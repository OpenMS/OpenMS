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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/AScore.h>


using namespace OpenMS;
using namespace std;

/**
	@page TOPP_PhosphoScoring PhosphoScoring

	@brief Tool to score phosphorylation sites of a peptide.
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ PhosphoScoring \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_??? (or other ID engines) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_??? </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_???? </td>
		</tr>
	</table>
</CENTER>

	TODO: docu

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PhosphoScoring.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPhosphoScoring
	: public TOPPBase
{
	public:
		TOPPPhosphoScoring()
			: TOPPBase("PhosphoScoring", "Scores potential phosphorylation sites and therby tries to localize the most probable sites.")
		{
		}

	protected:


		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "mzML", false);
			setValidFormats_("in", StringList::create("mzML"));
			registerInputFile_("id", "<file>", "", "Identification input file which contains a search against a concatenated sequence database", false);
			setValidFormats_("id", StringList::create("idXML"));
			registerOutputFile_("out", "<file>", "", "Identification output with annotated phosphorylation scores");
			registerDoubleOption_("fragment_mass_tolerance","<tolerance>",0.5,"Fragment mass error",false);
		//	registerStringOption_("type", "If set, the FDR of the proteins only is calculated");
		//	setValidFormats_("type",StringList::create("ascore"));
			//registerSubsection_("algorithm","Parameter section for the FDR calculation algorithm");

			addEmptyLine_();
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
			//Param alg_param = getParam_().copy("algorithm:", true);
			
			String in(getStringOption_("in"));
			String id(getStringOption_("id"));
			String out(getStringOption_("out"));
			DoubleReal fragment_mass_tolerance(getDoubleOption_("fragment_mass_tolerance"));
			//String type(getStringOption_("type"));
      AScore scoring_function;
	    MzMLFile f;
	    f.setLogType(log_type_);
     	RichPeakMap exp;
     	//-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------
				vector<PeptideIdentification> pep_ids;
				vector<ProteinIdentification> prot_ids;
				vector<PeptideIdentification> pep_out;
				IdXMLFile().load(id, prot_ids, pep_ids);
				MzMLFile().load(in,exp);				
      //-------------------------------------------------------------
     	// calculations
     	//-------------------------------------------------------------
			// map the ids to the spectra
			IDMapper id_mapper;
			id_mapper.annotate(exp, pep_ids, prot_ids); 
			for(RichPeakMap::iterator it = exp.begin(); it != exp.end(); ++it)
			{
				if(it->getPeptideIdentifications().size() == 0)
				{
					continue;
				}

				for(vector<PeptideIdentification>::iterator hits = it->getPeptideIdentifications().begin(); hits < it->getPeptideIdentifications().end(); ++hits)
				{
					vector<PeptideHit> scored_peptides;
					for(vector<PeptideHit>::const_iterator hit = hits->getHits().begin(); hit < hits->getHits().end(); ++hit)
					{
						PeptideHit scored_hit = *hit;
						RichPeakSpectrum& temp  = *it;
						//compute number of possible phosphorylation sites
						/*String without_phospho_str(scored_hit.getSequence().toString());
						size_t found = without_phospho_str.find("(Phospho)");
						while(found != string::npos)
						{
							without_phospho_str.erase(found,String("(Phospho)").size());
							found = without_phospho_str.find("(Phospho)");
						}
						AASequence without_phospho(without_phospho_str);
						DoubleReal prec = hits->getMetaValue("MZ");
						DoubleReal prec_mz =prec * scored_hit.getCharge();
						prec_mz -= scored_hit.getCharge();
						DoubleReal mono_weight = without_phospho.getMonoWeight();
						DoubleReal ha = prec_mz - mono_weight;
						DoubleReal nps = ha / 79.966331 ;*/

						//if(fabs( 1- nps)<= 0.1 || nps >= 1 )
						//{
							scoring_function.compute(scored_hit, temp, fragment_mass_tolerance);
						//}
						scored_peptides.push_back(scored_hit);
					}

					PeptideIdentification new_hits(*hits);
					new_hits.setHits(scored_peptides);
					pep_out.push_back(new_hits);
				}
			}
			//-------------------------------------------------------------
     	// writing output
     	//-------------------------------------------------------------
			
			IdXMLFile().store(out,prot_ids,pep_out);
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPPhosphoScoring tool;
	return tool.main(argc,argv);
	//Teste den hier es müssten 3 verschiedene theoretische Modelle entstehen = GCTSERGRYRIL
}

