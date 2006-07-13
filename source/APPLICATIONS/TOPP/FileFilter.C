// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/FORMAT/MzDataFile.h>

#include "TOPPBase.h"

using namespace OpenMS;
using namespace std;

/**
	@page FileFilter FileFilter
	
	@brief Extracts portions of the data from an mzData file.
	
	With this tool it is possible to exctract m/z, retention time and intensity ranges from a input mzData file
	and to write all data that lies within the given ranges to an output mzData file. It can also extract spectra
	of a certain MS level e.g. MS/MS spectra when using level '2'.
	
	@ingroup TOPP
*/


// We do not want this class to show up in the docu -> @cond
/// @cond TOPPCLASSES 

class TOPPFileFilter
	: public TOPPBase
{
	public:
		TOPPFileFilter()
			: TOPPBase("FileFilter")
		{
			
		}
	
	protected:
		void printToolUsage_()
		{
			cerr  << endl
						<< tool_name_ << " -- extracts portions of the data from an mzData file." << endl
						<< endl
						<< "Usage:" << endl
						<< " " << tool_name_ << " [options]" << endl
						<< endl
						<< "Options are:" << endl
						<< "  -in <file>        input mzData file name" << endl
						<< "  -out <file>       output mzData file name" << endl
						<< "  -mz [min]:[max]   m/z range to extract" << endl
						<< "  -rt [min]:[max]   retention time range to extract" << endl
						<< "  -int [min]:[max]  intensity range to extract" << endl
						<< "  -level i[,j]...   MS levels to extract (default: ALL)" << endl;
		}
	
		void printToolHelpOpt_()
		{
			cerr << endl
		       << tool_name_ << endl
		       << endl
		       << "INI options:" << endl
					 << "  in      input mzData file name" << endl
					 << "  out     output mzData file name" << endl
					 << "  mz      m/z range to extract" << endl
					 << "  rt      retention time range to extract" << endl
					 << "  int     intensity range to extract" << endl
					 << "  level   MS levels to extract (default: ALL)" << endl
					 << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"mz\" value=\"500:1000\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"rt\" value=\":100\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"int\" value=\"5000:\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"level\" value=\"1,2\" type=\"string\"/>" << endl;
		}
	
		void setOptionsAndFlags_()
		{
			options_["-out"] = "out";
			options_["-in"] = "in";
			options_["-mz"] = "mz";
			options_["-rt"] = "rt";
			options_["-int"] = "int";
			options_["-level"] = "level";
		}
	
		ExitCodes main_(int , char**)
		{

			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input file names and types
			String in = getParamAsString_("in");			
			writeDebug_(String("Input file: ") + in, 1);
	
			//output file names and types
			String out = getParamAsString_("out");
			writeDebug_(String("Output file: ") + out, 1);

			//ranges
			String mz, rt, it, level, tmp;
			double mz_l, mz_u, rt_l, rt_u, it_l, it_u;
			vector<UnsignedInt> levels;
			
			//initialize ranges
			mz_l = rt_l = it_l = -1 * numeric_limits<double>::max();
			mz_u = rt_u = it_u = numeric_limits<double>::max();
			
			//determine rt bounds
			rt = getParamAsString_("rt",":");
			writeDebug_(String("rt bounds: ") + rt,2);	
			
			//determine mz bounds
			mz = getParamAsString_("mz",":");
			writeDebug_(String("mz bounds: ") + mz,2);	
			
			//determine int bounds
			it = getParamAsString_("int",":");
			writeDebug_(String("int bounds: ") + it,2);	

			//determine levels
			level = getParamAsString_("level","1,2,3,4,5");
			writeDebug_(String("MS levels: ") + level,2);	
			
			//convert bounds to numbers
			try
			{
				//rt
				tmp = rt.prefix(':');
				if (tmp!="")
				{
					rt_l = tmp.toDouble();
				}
				tmp = rt.suffix(':');
				if (tmp!="")
				{
					rt_u = tmp.toDouble();
				}
				writeDebug_("rt lower/upper bound: " + String(rt_l) + " / " + String(rt_u),1);	
				
				//mz
				tmp = mz.prefix(':');
				if (tmp!="")
				{
					mz_l = tmp.toDouble();
				}
				tmp = mz.suffix(':');
				if (tmp!="")
				{
					mz_u = tmp.toDouble();
				}
				writeDebug_("mz lower/upper bound: " + String(mz_l) + " / " + String(mz_u),1);	
				
				//int
				tmp = it.prefix(':');
				if (tmp!="")
				{
					it_l = tmp.toDouble();
				}
				tmp = it.suffix(':');
				if (tmp!="")
				{
					it_u = tmp.toDouble();
				}
				writeDebug_("int lower/upper bound: " + String(it_l) + " / " + String(it_u),1);	
	
				//levels
				tmp = level;
				if (level.has(',')) //several levels given
				{
					vector<String> tmp2;
					level.split(',',tmp2);
					for (vector<String>::iterator it = tmp2.begin(); it != tmp2.end(); ++it)
					{
						levels.push_back(it->toInt());
					}
				}
				else //one level given
				{
					levels.push_back(level.toInt());
				}
				
				String tmp3("MS levels: ");
				tmp3 = tmp3 + String(*(levels.begin()));
				for (vector<UnsignedInt>::iterator it = ++levels.begin(); it != levels.end(); ++it)
				{
					tmp3 = tmp3 + ", " + String(*it);
				}
				writeDebug_(tmp3,1);	
			}
			catch(Exception::ConversionError& e)
			{
				writeLog_(String("Invalid boundary '") + tmp + "' given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;			
			}
			
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			
			MSExperiment< > exp;
			MzDataFile f;
			f.load(in,exp);						
		
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			
			//remove ms level first (might be a large amount of spectra)
			exp.erase(remove_if(exp.begin(), exp.end(), MSLevelRange<MSExperiment< >::SpectrumType>(levels, true)), exp.end());
			
			//remove rt range (discards whole spectra)
			exp.erase(remove_if(exp.begin(), exp.end(), RTRange<MSExperiment< >::SpectrumType>(rt_l, rt_u, true)), exp.end());
			
			
			for (MSExperiment< >::iterator it = exp.begin(); it!= exp.end(); ++it)
			{
				//remove int range (might be a lot more than mz)
				it->getContainer().erase(remove_if(it->begin(), it->end(), IntensityRange<MSExperiment< >::PeakType>(it_l, it_u, true)), it->end());
				
				//remove mz range
				it->getContainer().erase(remove_if(it->begin(), it->end(), MzRange<MSExperiment< >::PeakType>(mz_l, mz_u, true)), it->end());
			}
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			f.store(out,exp);
			
			return OK;
		}
};

/// @endcond


int main( int argc, char ** argv )
{
	TOPPFileFilter tool;
	return tool.main(argc,argv);
}

