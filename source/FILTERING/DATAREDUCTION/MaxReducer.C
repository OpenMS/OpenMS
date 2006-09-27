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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/DATAREDUCTION/MaxReducer.h>
using namespace std;
namespace OpenMS
{
	
  MaxReducer::MaxReducer(): DataReducer()
  {	
		name_= MaxReducer::getName();
  }
  MaxReducer::~MaxReducer()
  {
  }
 	void MaxReducer::applyReduction(const ExperimentType& in, ExperimentType& out)
 	{
		string name = "max_reduced_";
		name+=in.getName();
		out.setName(name);	
		double ratio = param_.getValue("Ratio");
 		ratio = ratio * 0.01;
		std::vector<PeakT> mz_values;
		SpectrumType base;
		 for(ExperimentType::ConstIterator spec_it = in.RTBegin(in.getMinRT()); 
				spec_it !=in.RTEnd(in.getMaxRT()); 
				++spec_it)
		 {
			 UnsignedInt  mz_counter = 1;
			 double ratio1 = (double)spec_it->size() * ratio;
			 if(ratio1 <=1)
			 {
				 ratio1 = 1;
			 }
			 for(SpectrumType::ConstIterator it  = spec_it->begin(); 
					 it!=spec_it->end(); 
					 ++it)
			 {
				 mz_values.push_back(*it);
				 if(mz_counter%(int)ratio1 ==0 ^ mz_counter == spec_it->size()+1)
				 {
					 if(mz_values.size() == 0)
					 {
						 mz_values.push_back(*it);	
					 }
					 BaseSpectrum::Iterator it1 = findMaxIntensity(mz_values);
					 base.push_back(*it1);		
					 base.setRetentionTime(spec_it->getRetentionTime(),0,0);
					 base.setMSLevel(spec_it->getMSLevel()); 
					 mz_values.erase(mz_values.begin(),mz_values.end());
				 }
				 mz_counter++;
				 
			 }
			 
			if(mz_values.size()!=0)
			{
				mz_values.erase(mz_values.begin(),mz_values.end());
			}	
			if(!base.empty())
			{	
				out.push_back(base);
				base.erase(base.begin(),base.end());
			}
		 }
		 out.updateRanges();
		 
	}
	MaxReducer::BaseSpectrum::Iterator  MaxReducer::findMaxIntensity(std::vector<PeakT >& peaks)
	{
		BaseSpectrum::Iterator it = peaks.begin();
		std::vector<PeakT >::iterator i;
		for( i = peaks.begin();i!=peaks.end();i++)
		{
			if(i->getIntensity()>it->getIntensity())
			{
				it = i;
			}
		}
		return it;
	}
	
}// namespace openms
