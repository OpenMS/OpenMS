// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzDataFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page MapNormalizer MapNormalizer
	
	@brief Normalizes peak intensities to the percentage of the maximum intensity in the HPLC-MS map.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapNormalizer
	: public TOPPBase
{
	public:
		TOPPMapNormalizer()
			: TOPPBase("MapNormalizer","Normalizes peak intensities in an MS run.")
		{
			
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file ");
			setValidFormats_("in",StringList::create("mzData"));
			registerOutputFile_("out","<file>","","output file ");
	  	setValidFormats_("out",StringList::create("mzData"));
		}
	
		ExitCodes main_(int , const char**)
		{

			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			String in = getStringOption_("in");
			String out = getStringOption_("out");
			
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			
			MSExperiment<Peak1D> exp;
			MzDataFile f;
			f.load(in,exp);						
		
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			
			//determine maximum peak
			exp.updateRanges();
			DoubleReal max = exp.getMaxInt() / 100.0;
			
			for (MSExperiment<Peak1D>::Iterator it = exp.begin(); it!= exp.end(); ++it)
			{
				if (it->getMSLevel() < 2)
				{
					for (MSExperiment<Peak1D>::SpectrumType::Iterator it2 = it->begin(); it2!= it->end(); ++it2)
					{
						it2->setIntensity( it2->getIntensity() / max);
					}
				}
			}
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			f.store(out,exp);
			
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPMapNormalizer tool;
	return tool.main(argc,argv);
}

/// @endcond
