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

#include <OpenMS/FILTERING/DATAREDUCTION/SumReducer.h>

using namespace std;

namespace OpenMS
{

  SumReducer::SumReducer(): DataReducer()
  {
		setName(SumReducer::getProductName());

		defaults_.setValue("range_per_step",2,"m/z window in which the peaks are summed up.",false);
		
		defaultsToParam_();		
  }
  
  SumReducer::~SumReducer()
  {
  }
  
	void SumReducer::applyReduction(const ExperimentType& in, ExperimentType& out)
 	{
		double reduction = (double)param_.getValue("range_per_step") ;
	// 	cout<<reduction<<endl;
		double distance;
		double sum;
		SpectrumType::ConstIterator begin;
		SpectrumType::ConstIterator end;
		SpectrumType::ConstIterator max;
		
		//experiment size
		out.resize(in.size());
		UInt out_spec = 0;
		
		for(ExperimentType::ConstIterator spec_it = in.begin(); spec_it !=in.end(); ++spec_it)
		{
			out[out_spec].setRT(spec_it->getRT());
			out[out_spec].setMSLevel(spec_it->getMSLevel()); 	
			
			distance  = std::max(((spec_it->end()-1)->getMZ()- spec_it->begin()->getMZ()) * reduction,0.01);
			begin  = spec_it->begin();
			end  = begin;
			sum = 0.0;
			
			while (end != spec_it->end())
			{
				while (end != spec_it->end() && end->getMZ() <= begin->getMZ()+distance )
				{
					++end;
				}
				if (begin==end)
				{
					continue;
				}
				max = begin;
				sum = 0;
				while (begin != end)
				{
					sum += begin->getIntensity();
					if (begin->getIntensity() > max->getIntensity())
					{
						max = begin;
					}
					++begin;
				}
				out[out_spec].push_back(*max);
			// 	cout<<*max<<endl;
// 				cout<<sum<<endl;

	// 			cout<<spec_it->getRT()<<endl;
				out[out_spec].back().setIntensity(sum);
			}
			if(!out[out_spec].empty())
			{	
				++out_spec;
			}
		}
		out.resize(out_spec);
		out.updateRanges(); 
// 		cout<<param_.getValue("range_per_step")<<"\t"<<out.getSize()<<endl;
	}
	
}// namespace openms

