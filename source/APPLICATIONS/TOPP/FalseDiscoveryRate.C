// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_FalseDiscoveryRate FalseDiscoveryRate
	
	@brief Tool to estimate the false discovery rate on peptide and protein level

	This TOPP tool can calulate the false discovery rate (FDR) given a forward and
	backward search. Most useful is this on protein level, however, it also can be 
	applied to peptides.

	The false discovery rate is defined as the number of false discoveries (the hits
	in the reversed search) over the number of false and correct discoveries (the hits 
	in both databases) given a score.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FalseDiscoveryRate.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFalseDiscoveryRate
	: public TOPPBase
{
	public:
		TOPPFalseDiscoveryRate()
			: TOPPBase("FalseDiscoveryRate", "Estimates the false discovery rate on peptide and protein level using decoy searches.")
		{
		}
	
	protected:

		Param getSubsectionDefaults_(const String& /*section*/) const
		{
			return FalseDiscoveryRate().getDefaults();
		}

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "Identification input file which contains a search against a concatenated sequence databse", false);
			setValidFormats_("in", StringList::create("idXML"));
			registerInputFile_("fwd_in", "<file>", "", "Identification input to estimate FDR, forward", false);
			setValidFormats_("fwd_in", StringList::create("idXML"));
			registerInputFile_("rev_in", "<file>", "", "Identification input to estimate FDR, decoy run", false);
			setValidFormats_("rev_in", StringList::create("idXML"));
			registerOutputFile_("out", "<file>", "", "Identification output with annotated FDR");
			registerFlag_("proteins_only", "if set, the FDR of the proteins only is calculated");
			registerFlag_("peptides_only", "if set, the FDR of the peptides only is caluclated");
			registerFlag_("q_value", "if set, the q-values will be calculated instead of the FDRs");
		
			addEmptyLine_();		
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			Param alg_param = getParam_().copy("q_value",false);
			FalseDiscoveryRate fdr;

			if (!alg_param.empty())
			{
				fdr.setParameters(alg_param);
				writeDebug_("Parameters passed to FalseDiscoveryRate", alg_param, 3);
			}

			//input/output files
			// either fwd_in and rev_in must be given or just the in which contains results of a search against a concatenated target decoy sequence db
			String fwd_in(getStringOption_("fwd_in")), rev_in(getStringOption_("rev_in")), in(getStringOption_("in"));
			bool combined(false);
			if (fwd_in != "" && rev_in != "")
			{
				if (in != "")
				{
					writeLog_("Error, either 'fwd_in' and 'rev_in' must be given or 'in', but not both");
					return ILLEGAL_PARAMETERS;
				}
			}
			else
			{
				if (in != "")
				{
					combined = true;
				}
				else
				{
					writeLog_("Error, at least 'fwd_in' and 'rev_in' or 'in' must be given");
					return ILLEGAL_PARAMETERS;
				}
			}
			String out(getStringOption_("out"));
			bool proteins_only(getFlag_("proteins_only"));
			bool peptides_only(getFlag_("peptides_only"));

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

			if (!combined)
			{
				vector<PeptideIdentification> fwd_pep, rev_pep;
				vector<ProteinIdentification> fwd_prot, rev_prot;
				String document_id;
				IdXMLFile().load(fwd_in, fwd_prot, fwd_pep, document_id);
				IdXMLFile().load(rev_in, rev_prot, rev_pep, document_id);
			
      	//-------------------------------------------------------------
      	// calculations
      	//-------------------------------------------------------------
			
				writeDebug_("Starting calculations", 1);

				if (!proteins_only)
				{
					fdr.apply(fwd_pep, rev_pep);
				}
				if (!peptides_only)
				{
					fdr.apply(fwd_prot, rev_prot);
				}

				//-------------------------------------------------------------
      	// writing output
      	//-------------------------------------------------------------

      	IdXMLFile().store(out, fwd_prot, fwd_pep);
			}
			else
			{
				vector<PeptideIdentification> pep_ids;
				vector<ProteinIdentification> prot_ids;
				IdXMLFile().load(in, prot_ids, pep_ids);
				if (!proteins_only)
				{
					fdr.apply(pep_ids);
				}
				if (!peptides_only)
				{
					fdr.apply(prot_ids);
				}

				IdXMLFile().store(out, prot_ids, pep_ids);
			}

			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPFalseDiscoveryRate tool;
	return tool.main(argc,argv);
}

