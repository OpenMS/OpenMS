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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ProteinInference ProteinInference

	@brief Computes a protein identification based on the FDRs of peptides.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ProteinInference.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPProteinInference
	: public TOPPBase
{
	public:
		TOPPProteinInference()
			: TOPPBase("ProteinInference","Protein inference based on FDRs of peptides.")
		{
		}

	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file");
			setValidFormats_("in",StringList::create("idXML"));
			registerOutputFile_("out","<file>","","output file");
			setValidFormats_("out",StringList::create("idXML"));

			addEmptyLine_();
			registerIntOption_("min_peptides_per_protein", "<num>", 2, "Minimal number of peptides needed for a protein identification", false);
			setMinInt_("min_peptides_per_protein", 1);

			registerSubsection_("algorithm","Consensus algorithm section");
		}

		ExitCodes main_(int , const char**)
		{
			String in = getStringOption_("in");
			String out = getStringOption_("out");

			vector<ProteinIdentification> prot_ids;
			vector<PeptideIdentification> pep_ids;
			IdXMLFile().load(in, prot_ids, pep_ids);

			


			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPProteinInference tool;
	return tool.main(argc,argv);
}

/// @endcond
