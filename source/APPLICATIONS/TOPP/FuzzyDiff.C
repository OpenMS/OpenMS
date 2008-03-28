// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FuzzyDiff FuzzyDiff
	
	@brief Compares two files.  Numeric differences are tolerated up to a given relative or absolute error.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFuzzyDiff
	: public TOPPBase
{
 public:
	TOPPFuzzyDiff()
		: TOPPBase("FuzzyDiff","Compares two files, tolerating numeric differences.")
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in1","<file>","","first input file",true);
		registerInputFile_("in2","<file>","","second input file",true);
		registerDoubleOption_("ratio","<double>",1,"acceptable relative error",false);
		setMinFloat_("ratio",1);
		registerDoubleOption_("absdiff","<double>",0,"acceptable absolute difference",false);
		setMinFloat_("absdiff",0);
		registerIntOption_("verbose","<int>",2,"set verbose level:\n"
											 "0 = very quiet mode (absolutely no output)\n"
											 "1 = quiet mode (no output unless differences detected)\n"
											 "2 = default (include summary at end)\n"
											 "3 = continue after errors\n",false
											);
		setMinInt_("verbose",0);
		setMaxInt_("verbose",3);
	}
	
	ExitCodes main_(int , const char**)
	{
	
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		String in1 = getStringOption_("in1");

		String in2 = getStringOption_("in2");
		
		double acceptable_ratio = getDoubleOption_("ratio");

		double acceptable_absdiff = getDoubleOption_("absdiff");

		int verbose_level = getIntOption_("verbose");

		OpenMS::FuzzyStringComparator fsc;
		
		fsc.setAcceptableRelative(acceptable_ratio);
		fsc.setAcceptableAbsolute(acceptable_absdiff);
		fsc.setVerboseLevel(verbose_level);

		if ( fsc.compare_files(in1,in2) )
		{
			return EXECUTION_OK;
		}
		else
		{
			// TODO think about better exit codes.
			return PARSE_ERROR;
		}
	}
};

int main( int argc, const char** argv )
{
	TOPPFuzzyDiff tool;
	return tool.main(argc,argv);
}

/// @endcond



