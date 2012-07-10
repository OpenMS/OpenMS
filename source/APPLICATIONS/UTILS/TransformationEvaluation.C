// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Hendrik Weisser
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_TransformationEvaluation TransformationEvaluation
	
	@brief Applies a transformation to a range of values and records the results.

	This is useful for plotting transformations for quality assessment etc.

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_TransformationEvaluation.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude UTILS_TransformationEvaluation.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPTransformationEvaluation
	: public TOPPBase
{
public:

	TOPPTransformationEvaluation()
		: TOPPBase("TransformationEvaluation", "Applies a transformation to a range of values", false)
	{
	}

protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "Input file containing the transformation description");
		setValidFormats_("in", StringList::create("trafoXML"));
    registerOutputFile_("out", "<file>", "", "Output file containing original and transformed values; if empty, output is written to the screen", false);
		registerDoubleOption_("min", "<value>", 0.0, "Minimum value to transform", false);
		registerDoubleOption_("max", "<value>", 0.0, "Maximum value to transform (if at or below 'min', select a suitable maximum based on the transformation description)", false);
		registerDoubleOption_("step", "<value>", 1.0, "Step size between 'min' and 'max'", false);
		setMinFloat_("step", 0.001);
	}


	ExitCodes main_(int, const char**)
	{
		String in = getStringOption_("in"), out = getStringOption_("out");

		TransformationDescription trafo_in;
		TransformationXMLFile().load(in, trafo_in);
		TransformationDescription::DataPoints data;

		DoubleReal min = getDoubleOption_("min"), max = getDoubleOption_("max"),
			step = getDoubleOption_("step");
		if (max <= min)
		{
			data = trafo_in.getDataPoints();
			sort(data.begin(), data.end());
			max = data.back().first;
			DoubleReal magnitude = floor(log10(max));
			max = Math::ceilDecimal(max, magnitude - 1);
			if (max <= min) 
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'min' must be lower than 'max'");
			}
		}

		data.clear();
		for (DoubleReal value = min; value <= max; value += step)
		{
			DoubleReal transformed = trafo_in.apply(value);
			if (out.empty()) 
			{
				cout << value << '\t' << transformed << endl;
			}
			else data.push_back(make_pair(value, transformed));
		}

		if (!out.empty())
		{
			TransformationDescription trafo_out(trafo_in);
			trafo_out.setDataPoints(data);
			TransformationXMLFile().store(out, trafo_out);
		}

		return EXECUTION_OK;
	}
};


int main(int argc, const char** argv)
{
	TOPPTransformationEvaluation tool;
	return tool.main(argc, argv);
}
  
/// @endcond
