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

#include <OpenMS/FILTERING/DATAREDUCTION/SumReducer.h>
using namespace std;
namespace OpenMS
{

  SumReducer::SumReducer(): DataReducer()
  {
		name_ = SumReducer::getName();
  }
  SumReducer::~SumReducer()
  {
  }
	void SumReducer::applyReduction(const ExperimentType& in, ExperimentType& out)
 	{
		string name = "sum_reduced_";
		name+=in.getName();
		out.setName(name);
		int ratio = param_.getValue("Ratio");
		double reduction = (double)ratio * 0.01;
		std::vector<PeakT> mz_values;
		SpectrumType base;
		
		for(ExperimentType::ConstIterator spec_it = in.RTBegin(in.getMinRT()); 
				spec_it !=in.RTEnd(in.getMaxRT()); 
				++spec_it)
		{
			int mz_counter =1;
			
			double  distance = (spec_it->end()-1)->getPosition()[0]- spec_it->begin()->getPosition()[0];
			double reduction1 =distance * reduction;
			if(reduction1<=1)
			{
				reduction1 = 1;
			}
			for(SpectrumType::ConstIterator it = spec_it->begin(); 
					it!=spec_it->end(); 
					++it)
			{
				mz_values.push_back(*it);
				
				if(it->getPosition()[0]>=(spec_it->begin()->getPosition()[0]+ (double)mz_counter * reduction1) ^ 
					 it->getPosition()[0] == (spec_it->end()-1)->getPosition()[0])
				{
					if(mz_values.size() == 0)
					{
						mz_values.push_back(*it);	
					}
					base.push_back(findSumIntensity(mz_values));		
					base.setRetentionTime(spec_it->getRetentionTime(),0,0);
					base.setMSLevel(spec_it->getMSLevel()); 						
					mz_values.erase(mz_values.begin(),mz_values.end());
					mz_counter++;
					
				}
			}
			if(!base.empty())
			{
				out.push_back(base);
				base.erase(base.begin(),base.end());
			}
		}
		out.updateRanges();
	}

	SumReducer::PeakT  SumReducer::findSumIntensity(std::vector<PeakT >& peaks)
	{
	 	double sumintensity = 0;
 		BaseSpectrum::Iterator it = peaks.begin();
		std::vector<PeakT >::iterator i;
 		for( i = peaks.begin();i!=peaks.end();i++)
 		{
 			sumintensity = sumintensity + i->getIntensity();
			if(i->getIntensity()>=it->getIntensity())
			{	
				it = i;
			}
		}
		it->setIntensity(sumintensity);
		return *it;
}

	
}// namespace openms

