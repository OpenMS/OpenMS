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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

#include<OpenMS/FORMAT/DTA2DFile.h>
#include<OpenMS/FORMAT/DFeatureMapFile.h>

#include<OpenMS/KERNEL/MSExperiment.h>

#include<OpenMS/DATASTRUCTURES/String.h> 

#include<OpenMS/KERNEL/DFeature.h>
#include<OpenMS/KERNEL/DFeatureMap.h>

using namespace OpenMS;
using namespace std;

/**
  
  @brief Short example application for the feature finder algorithm
    
 **/ 
int main(int argc, char *argv[])
{
	
	String filename;
	String inifile;
	
	if (argc<3)
	{
		std::cout << "Please provide a dta2d and INI file." << std::endl;
		return 1;
	} 
	else
	{
		filename = String(argv[1]);
		inifile  = String(argv[2]);
	}
	
	Param feafi_params;
	feafi_params.load(inifile);
		
	/// read data file
	DTA2DFile dta2d_file;
	MSExperiment<DPeak<1> > exp;
	dta2d_file.load(filename,exp);
	
	// Initialize feature finder 
	FeatureFinder ff;
	ff.setParam(feafi_params);
	ff.setData(exp);
		
	// run it...
	DFeatureMap<2> features;
	features = ff.run();
	
	// write features to file
	DFeatureMapFile map_file;
	map_file.store("Features.xml",features);

	return 0;
}
