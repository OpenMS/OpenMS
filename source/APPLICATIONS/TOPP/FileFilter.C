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

#include <OpenMS/APPLICATIONS/TOPPBase2.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FileFilter FileFilter
	
	@brief Extracts portions of the data from an mzData file.
	
	With this tool it is possible to exctract m/z, retention time and intensity ranges from a input mzData file
	and to write all data that lies within the given ranges to an output mzData file. It can also extract spectra
	of a certain MS level e.g. MS/MS spectra when using level '2'.
	
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileFilter
	: public TOPPBase2
{
	public:
		TOPPFileFilter()
			: TOPPBase2("FileFilter","extracts portions of the data from an mzData file")
		{
			
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerStringOption_("in","<file>","","input file in MzData format");
			registerStringOption_("out","<file>","","output file in MzData format");
			registerStringOption_("mz","[min]:[max]",":","m/z range to extract", false);
			registerStringOption_("rt","[min]:[max]",":","retention time range to extract", false);
			registerStringOption_("int","[min]:[max]",":","intensity range to extract", false);
			registerStringOption_("level","-level i[,j]...","1,2,3","MS levels to extract", false);
			registerFlag_("remove_zoom","flag that removes zoom scans");
		}
	
		ExitCodes main_(int , char**)
		{

			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			String in = getStringOption_("in");			
			String out = getStringOption_("out");
	
			//ranges
			String mz, rt, it, level, tmp;
			double mz_l, mz_u, rt_l, rt_u, it_l, it_u;
			vector<UnsignedInt> levels;		
			//initialize ranges
			mz_l = rt_l = it_l = -1 * numeric_limits<double>::max();
			mz_u = rt_u = it_u = numeric_limits<double>::max();
			
			rt = getStringOption_("rt");
			mz = getStringOption_("mz");
			it = getStringOption_("int");
			level = getStringOption_("level");
			
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
				tmp3 = tmp3 + *(levels.begin());
				for (vector<UnsignedInt>::iterator it = ++levels.begin(); it != levels.end(); ++it)
				{
					tmp3 = tmp3 + ", " + *it;
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
			
			//remove zoom scan mode (might be a lot of spectra)
			bool rem_zoom = getFlag_("remove_zoom");
			writeDebug_(String("Remove zoom: ") + String(rem_zoom),3);
			if (rem_zoom)
			{
				exp.erase(remove_if(exp.begin(), exp.end(), ScanModePredicate<MSExperiment< >::SpectrumType>(InstrumentSettings::SELECTEDIONDETECTION)), exp.end());
			}
			
			//remove rt range (discards whole spectra)
			exp.erase(remove_if(exp.begin(), exp.end(), RTRange<MSExperiment< >::SpectrumType>(rt_l, rt_u, true)), exp.end());
			
			
			for (MSExperiment< >::iterator it = exp.begin(); it!= exp.end(); ++it)
			{
				//remove int range (might be a lot more than mz)
				it->getContainer().erase(remove_if(it->begin(), it->end(), IntensityRange<MSExperiment< >::PeakType>(it_l, it_u, true)), it->end());
				
				//remove mz range
				it->getContainer().erase(remove_if(it->begin(), it->end(), MzRange<MSExperiment< >::PeakType>(mz_l, mz_u, true)), it->end());
			}
			
			//remove empty scans
			exp.erase(remove_if(exp.begin(), exp.end(), SpectrumEmptyPredicate<MSExperiment< >::SpectrumType>()), exp.end());
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			f.store(out,exp);
			
			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPFileFilter tool;
	return tool.main(argc,argv);
}

/// @endcond
