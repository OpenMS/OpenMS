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

#include <OpenMS/METADATA/ProcessingMethod.h>


using namespace std;

namespace OpenMS
{

	ProcessingMethod::ProcessingMethod():
		MetaInfoInterface(),
		deisotoping_(false),
		charge_deconvolution_(false),
		method_(SpectrumSettings::UNKNOWN),
		intensity_cutoff_(0.0)
	{
		
	}
	
	ProcessingMethod::ProcessingMethod(const ProcessingMethod& source):
		MetaInfoInterface(source),
	  deisotoping_(source.deisotoping_),
	  charge_deconvolution_(source.charge_deconvolution_),
	  method_(source.method_),
	  intensity_cutoff_(source.intensity_cutoff_)
	{
	  
	}
	
	ProcessingMethod::~ProcessingMethod()
	{
	  
	}
	
	ProcessingMethod& ProcessingMethod::operator = (const ProcessingMethod& source)
	{
	  if (&source == this) return *this;
	  
	  deisotoping_ = source.deisotoping_;
	  charge_deconvolution_ = source.charge_deconvolution_;
	  method_ = source.method_;
	  intensity_cutoff_ = source.intensity_cutoff_;
	  MetaInfoInterface::operator=(source);
	  
	  return *this;
	}
	
	bool ProcessingMethod::operator== (const ProcessingMethod& rhs) const
	{
		return 
			deisotoping_ == rhs.deisotoping_ &&
	    charge_deconvolution_ == rhs.charge_deconvolution_ &&
	    method_ == rhs.method_ &&
	    intensity_cutoff_ == rhs.intensity_cutoff_ &&
			MetaInfoInterface::operator==(rhs)
			;
	}
	
	bool ProcessingMethod::operator!= (const ProcessingMethod& rhs) const
	{
		return !(operator==(rhs));
	}
	
	bool ProcessingMethod::getDeisotoping() const 
	{
	  return deisotoping_; 
	}
	
	void ProcessingMethod::setDeisotoping(bool deisotoping)
	{
	  deisotoping_ = deisotoping; 
	}
	
	bool ProcessingMethod::getChargeDeconvolution() const 
	{
	  return charge_deconvolution_; 
	}
	
	void ProcessingMethod::setChargeDeconvolution(bool charge_deconvolution)
	{
	  charge_deconvolution_ = charge_deconvolution; 
	}
	
	SpectrumSettings::SpectrumType ProcessingMethod::getSpectrumType() const 
	{
	  return method_; 
	}
	
	void ProcessingMethod::setSpectrumType(SpectrumSettings::SpectrumType method)
	{
	  method_ = method; 
	}

	float ProcessingMethod::getIntensityCutoff() const
	{
		return intensity_cutoff_;
	}
	
	void ProcessingMethod::setIntensityCutoff(float cutoff)
	{
		intensity_cutoff_ = cutoff;
	}
}

