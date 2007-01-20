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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/DTAFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page DTAExtractor DTAExtractor
	
	@brief Extracts scans of an mzData file to several files in DTA format.
	
	The retention time, the m/z ratio (for MS level > 1) and the file extention are appended to the output file name.
	
	You can limit the exported spectra by m/z range, retention time range or MS level.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES 

class TOPPDTAExtractor
	: public TOPPBase
{
	public:
		TOPPDTAExtractor()
			: TOPPBase("DTAExtractor","extracts scans of an mzData file to several files in DTA format")
		{
			
		}

	protected:
		void registerOptionsAndFlags_()
		{
			registerStringOption_("in","<file>","","input file in MzData format");
			registerStringOption_("out","<file>","","base name of output files (RT, m/z and extension are appended)");
			registerStringOption_("mz","[min]:[max]",":","m/z range of precursor peaks to extract.\n"
																									 "This option is ignored for MS level 1", false);
			registerStringOption_("rt","[min]:[max]",":","retention time range of spectra to extract", false);
			registerStringOption_("level","i[,j]...","1,2,3","MS levels to extract", false);
		}
	
		ExitCodes main_(int , char**)
		{

			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			String in = getStringOption_("in");
			inputFileReadable_(in);	
			
			String out = getStringOption_("out");			

			//ranges
			String mz, rt, tmp;
			double mz_l, mz_u, rt_l, rt_u;
			vector<UnsignedInt> levels;			
			//initialize ranges
			mz_l = rt_l = -1 * numeric_limits<double>::max();
			mz_u = rt_u = numeric_limits<double>::max();
			
			rt = getStringOption_("rt");
			mz = getStringOption_("mz");
			String level = getStringOption_("level");
			
			//convert bounds to numbers
			try
			{
				//rt
				parseRange_(rt,rt_l,rt_u);
				writeDebug_("rt lower/upper bound: " + String(rt_l) + " / " + String(rt_u),1);	
				
				//mz
				parseRange_(mz,mz_l,mz_u);
				writeDebug_("mz lower/upper bound: " + String(mz_l) + " / " + String(mz_u),1);	
				
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
			f.getOptions().setRTRange(DRange<1>(rt_l,rt_u));
			f.load(in,exp);						

			DTAFile dta;
		
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			
			for (MSExperiment< >::iterator it = exp.begin(); it!= exp.end(); ++it)
			{
				//check for MS-level
				bool in_level_range = false;
				for (vector<UnsignedInt>::iterator it2 = levels.begin(); it2 != levels.end(); ++it2)
				{
					if (it->getMSLevel()==*it2)
					{
						in_level_range = true;
					}
				}
				if (!in_level_range) continue;
				
				//store spectra
				if (it->getMSLevel()>1)
				{
					double mz = it->getPrecursorPeak().getPosition()[0];
					if (mz<mz_l || mz>mz_u)
					{
						continue;
					}		
					dta.store(out+"_RT"+String(it->getRetentionTime())+"_MZ"+String(mz)+".dta", *it);
				}
				else
				{
					dta.store(out+"_RT"+String(it->getRetentionTime())+".dta", *it);
				}
			}
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, char ** argv )
{
	TOPPDTAExtractor tool;
	return tool.main(argc,argv);
}

