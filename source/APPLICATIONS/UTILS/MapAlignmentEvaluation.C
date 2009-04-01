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
// $Maintainer: Katharina Albers, Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page MapAlignmentEvaluation MapAlignmentEvaluation
		
	@brief Evaluate alignment results against a ground truth
	
	This tool implements the evaluation measures from our paper 
	"Critical assessment of alignment procedures for LC-MS proteomics and metabolomics measurements",
	Eva Lange, Ralf Tautenhahn, Steffen Neumann, Clemens Groepl. BMC Bioinformatics 2008, 9:375.
	doi:10.1186/1471-2105-9-375.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_MapAlignmentEvaluation.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignmentEvaluation
  : public TOPPBase
{

public:
	TOPPMapAlignmentEvaluation()
		: TOPPBase("MapAlignmentEvaluation","Evaluate alignment results against a ground truth.", false)
	{
	}

protected: 
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file: tool",true);
		setValidFormats_("in",StringList::create("consensusXML"));
		registerInputFile_("gt","<file>","","input file: ground truth",true);
		setValidFormats_("gt",StringList::create("consensusXML"));
		//registerOutputFile_("out","<file>","","output file ");
		registerStringOption_("type","<name>","","Caap Evaluation type",true);
		setValidStrings_("type",Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts());

		/*addEmptyLine_();
		addText_("This tool implements the evaluation measures from our paper:\n"
						 "\"Critical assessment of alignment procedures for LC-MS proteomics and metabolomics measurements\"\n"
						 "Eva Lange, Ralf Tautenhahn, Steffen Neumann, Clemens Groepl\n"
						 "BMC Bioinformatics 2008, 9:375.\n"
						 "doi:10.1186/1471-2105-9-375\n"
						); */

	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------

		String in   = getStringOption_("in");
		String gt   = getStringOption_("gt");
		String type = getStringOption_("type");

		DoubleReal out = 0;

		//-------------------------------------------------------------
		// check for valid input
		//-------------------------------------------------------------
		//check if both input files have the correct type
	/*	if (FileHandler::getType(in)!=FileHandler::CONSENSUSXML)
		{
			writeLog_("Error: The input file must be of type ConsensusXML!");
			return ILLEGAL_PARAMETERS;
		}

		if (FileHandler::getType(gt)!=FileHandler::CONSENSUSXML)
		{
			writeLog_("Error: The groundtruth file must be of type ConsensusXML!");
			return ILLEGAL_PARAMETERS;
		}	*/
		//fehler: /home/kate/OpenMS/source/APPLICATIONS/TOPP/MapAlignmentEvaluation.C: In member function ‘virtual OpenMS::TOPPBase::ExitCodes TOPPMapAlignmentEvaluation::main_(int, const char**)’:
		//home/kate/OpenMS/source/APPLICATIONS/TOPP/MapAlignmentEvaluation.C:114: error: ‘CONSENSUSXML’ is not a member of ‘OpenMS::FileHandler’
		//home/kate/OpenMS/source/APPLICATIONS/TOPP/MapAlignmentEvaluation.C:120: error: ‘CONSENSUSXML’ is not a member of ‘OpenMS::FileHandler’

		//-------------------------------------------------------------
		// set up algorithm
		//-------------------------------------------------------------
		MapAlignmentEvaluationAlgorithm* algorithm = Factory<MapAlignmentEvaluationAlgorithm>::create(type);

		
		//-------------------------------------------------------------
		// read input files
		//-------------------------------------------------------------

		// reader
		ConsensusXMLFile consensus_xml_file_in;
		consensus_xml_file_in.setLogType(log_type_);

		// tool -> consensus_map_in
		ConsensusMap consensus_map_in;
		consensus_xml_file_in.load( in, consensus_map_in );
		
		// reader
		ConsensusXMLFile consensus_xml_file_gt;
		consensus_xml_file_gt.setLogType(log_type_); //richtig??

		// gt -> consensus_map_gt
		ConsensusMap consensus_map_gt;
		consensus_xml_file_gt.load( gt, consensus_map_gt );
		

		//evaluate
		algorithm->evaluate(consensus_map_in, consensus_map_gt, out);

		//write output
		cout << type << ": " << out << "\n";

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPMapAlignmentEvaluation tool;
  return tool.main(argc,argv);
}

/// @endcond
