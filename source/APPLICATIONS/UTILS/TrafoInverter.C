// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_TrafoInverter TrafoInverter
	
	@brief Computes an (approximate) inverse of a retention time transformation.

	Inverting RT transformations becomes necessary in a certain data analysis pipeline for label-free quantification: After a first pass of @ref TOPP_FeatureFinder "feature detection", @ref TOPP_MapAligner "alignment" and @ref TOPP_FeatureLinker "feature grouping", seed lists generated from the consensus map (cf. @ref TOPP_SeedListGenerator "SeedListGenerator") should be used for a second pass of feature detection, in order to fill "holes" in the consensus map. Since the consensus map was produced from transformed (aligned) feature maps, the RT scales of the seed lists do not match those of the original LC-MS maps (mzML files). However, applying the RT transformations computed from the feature maps to the LC-MS maps is not possible, because it would impair the feature detection. Instead, the transformations have to be inverted, and the inverse transformations applied to the seed lists, before those can be used for targeted feature detection.
  	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_TrafoInverter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPTrafoInverter
	: public TOPPBase
{
public:

	TOPPTrafoInverter()
		: TOPPBase("TrafoInverter", "Computes an (approximate) inverse of a retention time transformation", false)
	{
	}


protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "Input file");
		setValidFormats_("in", StringList::create("trafoXML"));
		registerOutputFile_("out", "<file>", "", "Output file");
		setValidFormats_("out", StringList::create("trafoXML"));
	}


	ExitCodes main_(int , const char**)
	{
		String in = getStringOption_("in"), out = getStringOption_("out");
	
		TransformationDescription trafo;
		TransformationXMLFile trafo_file;
		trafo_file.load(in, trafo);
		trafo.getInverse(trafo);
		trafo_file.store(out, trafo);

		return EXECUTION_OK;
	}
};


int main(int argc, const char** argv)
{
	TOPPTrafoInverter tool;
	return tool.main(argc, argv);
}
  
/// @endcond
