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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/MaxReducer.h>
using namespace std;
namespace OpenMS
{
	
  MaxReducer::MaxReducer()
  	: DataReducer()
  {
		setName(MaxReducer::getProductName());

		defaults_.setValue("Peaksperstep",10);
		
		defaultsToParam_();		
  }
  
  MaxReducer::~MaxReducer()
  {
  }
  
 	void MaxReducer::applyReduction(const ExperimentType& in, ExperimentType& out)
 	{
		UnsignedInt peaks_per_bin = (UnsignedInt)(param_.getValue("Peaksperstep"));
 		
		out.resize(in.size());
		UnsignedInt out_spec = 0;

		//variables
	// 	UnsignedInt peaks_per_bin;
		UnsignedInt counter;
		SpectrumType::ConstIterator begin;
		SpectrumType::ConstIterator end;
		SpectrumType::ConstIterator max;	

		for(ExperimentType::ConstIterator spec_it = in.begin(); spec_it !=in.end(); ++spec_it)
		{
			out[out_spec].setRetentionTime(spec_it->getRetentionTime());
			out[out_spec].setMSLevel(spec_it->getMSLevel()); 
			
			//init
// // 			peaks_per_bin = std::max( (int)(spec_it->size() * ratio), 1);
// 			peaks_per_bin = ratio;
// 			cout<<peaks_per_bin<<endl;

			begin  = spec_it->begin();
			end  = begin;
			while (end != spec_it->end())
			{
				counter=0;
				while (end != spec_it->end() && counter < peaks_per_bin)
				{
					++end;
					++counter;
				}
				
				max = begin;
				while (begin != end)
				{
					if (begin->getIntensity() > max->getIntensity())
					{
						max = begin;
					}
					++begin;
				}
				out[out_spec].push_back(*max);
			}

			if(!out[out_spec].empty())
			{	
				++out_spec;
			}
		}
		out.resize(out_spec);
		out.updateRanges(); 
	}
	
}// namespace openms
