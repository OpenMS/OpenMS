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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

using namespace OpenMS;
using namespace std;

/**
	@page IDDecoyProbability IDDecoyProbability
	
	@brief Tool to estimate probability of peptide hits

	So far an estimation of the false score distribution with a gamma distribution
	and the correct score distribution with a gaussian distribution is performed. 
	The probabilities are calculated using bayes law, similar to PeptideProphet.
	This implementation is much simpler than that of PeptideProphet, however it is
	more robust and can deal with more search engines!
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDDecoyProbability
	: public TOPPBase
{
	public:
		TOPPIDDecoyProbability()
			: TOPPBase("IDDecoyProbability", "Estimates peptide probabilities using a decoy search strategy.")
		{
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("fwd_in", "<file>", "", "Identification input of forward run");
			registerInputFile_("rev_in", "<file>", "", "Identification input of decoy run");
			registerOutputFile_("out", "<file>", "", "Identification output with forward scores converted to probabilities");

			registerSubsection_("decoy_algorithm", "Algorithm parameter subsection");
			addEmptyLine_();		
		}

		Param getSubsectionDefaults_(const String& /*section*/) const
		{
			IDDecoyProbability decoy_prob;
			return decoy_prob.getParameters();
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input/output files
			String fwd_in(getStringOption_("fwd_in")), rev_in(getStringOption_("rev_in"));
			String out(getStringOption_("out"));

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

			vector<PeptideIdentification> fwd_pep, rev_pep, out_pep;
			vector<ProteinIdentification> fwd_prot, rev_prot;
			IdXMLFile().load(fwd_in, fwd_prot, fwd_pep);
			IdXMLFile().load(rev_in, rev_prot, rev_pep);
			
      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
			
			writeDebug_("Starting calculations", 1);

			IDDecoyProbability decoy_prob;
			Param decoy_param = getParam_().copy("decoy_algorithm:",true);
			decoy_prob.setParameters(decoy_param);
			decoy_prob.apply(out_pep, fwd_pep, rev_pep);
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
	

			IdXMLFile().store(out, fwd_prot, out_pep);
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPIDDecoyProbability tool;
	return tool.main(argc,argv);
}

