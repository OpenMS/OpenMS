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
// $Maintainer: Clemens Groepl, Chris Bielow $
// $Authors: Katharina Albers $
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
	@page UTILS_MapAlignmentEvaluation MapAlignmentEvaluation
		
	@brief Evaluates alignment results against a ground truth.
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MapAlignmentEvaluationAlgorithm \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (text output)</td>
		</tr>
  </table>
</CENTER>

	This tool implements the evaluation measures published in\n
	"Critical assessment of alignment procedures for LC-MS proteomics and metabolomics measurements",\n
	Eva Lange, Ralf Tautenhahn, Steffen Neumann, Clemens Groepl. BMC Bioinformatics 2008, 9:375.\n
	doi:10.1186/1471-2105-9-375.\n

	Input is a ground truth file as described on the CAAP web page\n
	Output is a recall- or a precision-value.\n
	
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
		: TOPPBase("MapAlignmentEvaluation","Evaluates alignment results against a ground truth.", false)
	{
	}

protected: 
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file: tool",true);
		setValidFormats_("in",StringList::create("consensusXML"));
		registerInputFile_("gt","<file>","","input file: ground truth",true);
		setValidFormats_("gt",StringList::create("consensusXML"));
		registerStringOption_("type","<name>","","Caap Evaluation type",true);
		StringList types = Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts();
		types.push_back("F1");
		setValidStrings_("type", types);
		registerDoubleOption_("rt_dev","<double>",0.1,"Maximum allowed deviation of the retention time", false);
		registerDoubleOption_("mz_dev","<double>",0.1,"Maximum allowed deviation of m/z", false);
		registerDoubleOption_("int_dev","<double>",100,"Maximum allowed deviation of Intensity", false);
		registerFlag_("use_charge", "Use charge criterion when assesing if two features are identical.", false);

		addEmptyLine_();
		addText_("This tool implements the evaluation measures published in:\n"
			"\"Critical assessment of alignment procedures for LC-MS proteomics and metabolomics measurements\"\n"
			"Eva Lange, Ralf Tautenhahn, Steffen Neumann, Clemens Groepl\n"
			"BMC Bioinformatics 2008, 9:375.\n"
			"doi:10.1186/1471-2105-9-375\n"
			);

	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------

		String in   = getStringOption_("in");
		String gt   = getStringOption_("gt");
		String type = getStringOption_("type");

		DoubleReal rt_dev = getDoubleOption_("rt_dev");
		DoubleReal mz_dev = getDoubleOption_("mz_dev");
		DoubleReal int_dev = getDoubleOption_("int_dev");

		bool use_charge = getFlag_("use_charge");

		//-------------------------------------------------------------
		// check for valid input
		//-------------------------------------------------------------
		//check if both input files have the correct type
		if (FileHandler::getType(in) != FileTypes::CONSENSUSXML)
		{
			writeLog_("Error: The input file must be of type ConsensusXML!");
			return ILLEGAL_PARAMETERS;
		}

		if (FileHandler::getType(gt) != FileTypes::CONSENSUSXML)
		{
			writeLog_("Error: The groundtruth file must be of type ConsensusXML!");
			return ILLEGAL_PARAMETERS;
		}
		
		//-------------------------------------------------------------
		// read input files
		//-------------------------------------------------------------

		// reader
		ConsensusXMLFile consensus_xml_file_in;
		consensus_xml_file_in.setLogType(log_type_);

		// tool -> consensus_map_in
		ConsensusMap consensus_map_in;
		consensus_xml_file_in.load( in, consensus_map_in );
		
		// gt -> consensus_map_gt
		ConsensusMap consensus_map_gt;
		consensus_xml_file_in.load( gt, consensus_map_gt );
		

		//-------------------------------------------------------------
		// set up algorithm
		//-------------------------------------------------------------
		if (type == "F1")
		{
      MapAlignmentEvaluationAlgorithm* algorithm_p = Factory<MapAlignmentEvaluationAlgorithm>::create("precision");
			MapAlignmentEvaluationAlgorithm* algorithm_r = Factory<MapAlignmentEvaluationAlgorithm>::create("recall");

      DoubleReal precision = 0;
      DoubleReal recall = 0;

      //evaluate
      algorithm_p->evaluate(consensus_map_in, consensus_map_gt, rt_dev, mz_dev, int_dev, use_charge, precision);
      algorithm_r->evaluate(consensus_map_in, consensus_map_gt, rt_dev, mz_dev, int_dev, use_charge, recall);

			//write output
      cout << "precision" << ": " << precision << "\n";
      cout << "   recall" << ": " << recall << "\n";
      cout << "-->    F1" << ": " << (2*precision*recall)/(precision+recall) << " (2*precision*recall)/(precision+recall)\n";
			
		}
		else
		{
			MapAlignmentEvaluationAlgorithm* algorithm = Factory<MapAlignmentEvaluationAlgorithm>::create(type);

      DoubleReal result = 0;

      //evaluate
      algorithm->evaluate(consensus_map_in, consensus_map_gt, rt_dev, mz_dev, int_dev, use_charge, result);

			//write output
      cout << type << ": " << result << "\n";
		}

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPMapAlignmentEvaluation tool;
  return tool.main(argc,argv);
}

/// @endcond
