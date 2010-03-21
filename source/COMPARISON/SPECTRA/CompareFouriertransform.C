// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Erhan Kenar $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------
//									
#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>

#include <cmath>



namespace OpenMS
{
	CompareFouriertransform::CompareFouriertransform()
  : PeakSpectrumCompareFunctor()
  {
		setName(CompareFouriertransform::getProductName());
		defaults_.setValue("epsilon", 0.2, "defines the absolut error of the mass spectrometer");
		defaultsToParam_();
  }
	
	CompareFouriertransform::CompareFouriertransform(const CompareFouriertransform& source)
  : PeakSpectrumCompareFunctor(source)
  {
  }

	CompareFouriertransform::~CompareFouriertransform()
  {
  }
	
	CompareFouriertransform& CompareFouriertransform::operator = (const CompareFouriertransform& source)
  {
		if (this != &source)
		{
    	PeakSpectrumCompareFunctor::operator = (source);
		}
    return *this;
  }
	double CompareFouriertransform::operator () (const PeakSpectrum& ) const
	{
		return 0;
	}
	double CompareFouriertransform::operator () (const PeakSpectrum& spec1 , const PeakSpectrum& spec2 ) const
	{
		const MSSpectrum<>::FloatDataArrays& temp1 = spec1.getFloatDataArrays();
		if(temp1.size()== 0)
		{
 
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Input need to be a fouriertransformation, try first transform ()");
			//transform(spec1); 
		}
	
		UInt i=	searchTransformation_(spec1);
		
		const MSSpectrum<>::FloatDataArrays& temp2 = spec2.getFloatDataArrays();
		if(temp2.size()== 0)
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Input need to be a fouriertransformation, try first transform ()");
		}
		UInt j=	searchTransformation_(spec2);
		
		if(temp1[i].size() != temp2[j].size())
		{
		//	std::cout<< temp1[i].size() << temp2[j].size() << std::endl; 
			return 0.0;
		}
		else
		{	
			
			Real sum=0;
			for (Size k=0; k< temp1[i].size();++k)
			{
		 
				sum=sum + temp1[i][k]-temp2[j][k];
			}
		//	std::cout << sum << " summe " << std::endl;
			if(sum !=0)
			{
				return 0;
			}
			else
			{
				return 1;
			}
		}
	}
	UInt CompareFouriertransform::searchTransformation_(const PeakSpectrum&  spec) const
  {
		const MSSpectrum<>::FloatDataArrays& temp = spec.getFloatDataArrays();
		UInt i=0;
		while(i< temp.size())
		{
			if(temp[i].getName()=="Fouriertransformation")
			{
				break;
			}
			else
			{
				++i;
			}
		}
		//find entry
		if(i == temp.size() || temp[i].getName()!= "Fouriertransformation")
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Input need to be a fouriertransformation, try first transform ()");			//transform_(spec);
			//return i+1;
		
		}
		else
		{
			return i;
		}
  }
  void CompareFouriertransform::transform(PeakSpectrum & spec)
  {
  
  	double* data = new double [2*spec.size()];
  	//normalize first the intensity!!!
  	DoubleReal int_sum=0;
  	for (Size p = 0 ;p<spec.size(); ++p)
  	{
  		int_sum+=spec[p].getIntensity();
  	}
  	//copy the peaks two times
  	Size i = 0;
  	for (Size p=0; p<spec.size(); ++p)
  	{
  		data[i] = spec[p].getIntensity()/int_sum;
  		++i;
  	}
  	for (SignedSize p=spec.size()-1; p>=0; --p)
  	{
  		data[i] = spec[p].getIntensity()/int_sum;
  		++i;
  	}
  	
  	gsl_fft_real_wavetable * real;
  	gsl_fft_real_workspace * work;
  	work = gsl_fft_real_workspace_alloc (spec.size());
  	real = gsl_fft_real_wavetable_alloc (spec.size());
  	gsl_fft_real_transform (data,1,spec.size(),real, work);
  	gsl_fft_real_wavetable_free (real);
  	gsl_fft_real_workspace_free (work);
  	
  	MSSpectrum<>::FloatDataArrays& temp = spec.getFloatDataArrays();
  	i= temp.size();
  	temp.resize(i+1);
  	temp[i].setName("Fouriertransformation");
  	UInt j=0;
  	while(j < spec.size())
  	{
  		temp[i].push_back(data[j]);
  		if(j==0) ++j; //we only intress in the real part of FFT
  		else j=j+2;
  	}
  	delete[] data;
  }

}

