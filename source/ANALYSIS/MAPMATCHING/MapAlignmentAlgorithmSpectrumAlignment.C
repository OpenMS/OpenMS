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
// $Maintainer: Vipul Patel $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>

namespace OpenMS
{

	MapAlignmentAlgorithmSpectrumAlignment::MapAlignmentAlgorithmSpectrumAlignment()
		: MapAlignmentAlgorithm(),ProgressLogger()
	{
		setName("MapAlignmentAlgorithmSpectrumAlignment");
		defaults_.setValue("gapcost",2,"gapcost",false);
		defaults_.setValue("affinegapcost", 1,"extenscion cost",false);
		defaults_.setValue("scorefunction","SteinScottImproveScore","scoring of mssprectren",false);
		setLogType(CMD);
		defaultsToParam_();
	}

	MapAlignmentAlgorithmSpectrumAlignment::~MapAlignmentAlgorithmSpectrumAlignment()
	{
	}
	
	void MapAlignmentAlgorithmSpectrumAlignment::alignPeakMaps(std::vector< MSExperiment<> >& peakmaps)
	{
 
		try
		{ 	
			updateMembers_();
			std::vector<UInt> pattern;
			std::vector<MSSpectrum<>* >versuch;			
			peakmaps[0].updateRanges(-7);
			pattern=peakmaps[0].getMSLevels();
				
				if(pattern.size()!=0)
				{	
					startProgress(0,(peakmaps.size()-1),"aligment");
					for(UInt i=0; i< peakmaps[0].size();++i)
					{
						if(peakmaps[0][i].getMSLevel()==1)
						{ 
							versuch.push_back(&(peakmaps[0][i]));
						}
					}
					for(UInt i = 1 ; i < peakmaps.size();++i )
									{									
									preparealign(versuch,peakmaps[i]);
									setProgress(i);
std::cout<< std::endl;
									
									}
					endProgress();
				}
				else
					{		
						//throw Exception::FileEmpty(__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				
				
		}
		catch (Exception::OutOfRange& e) 
		{
			throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);	
		}
			

	}
	void MapAlignmentAlgorithmSpectrumAlignment::updateMembers_()
	 {
		
		gap_	=(int)param_.getValue("gapcost");
		e_		=(int)param_.getValue("affinegapcost");
		c1 = Factory<PeakSpectrumCompareFunctor>::create((String)param_.getValue("scorefunction"));
		
	 }
	


} 

