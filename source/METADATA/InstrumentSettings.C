// -*- Mode: C++; tab-width: 2; -*-s
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

#include <OpenMS/METADATA/InstrumentSettings.h>

using namespace std;

namespace OpenMS
{

	const std::string InstrumentSettings::NamesOfScanMode[] = {"Unknown","SELECTEDIONDETECTION","MASSSCAN"};

	InstrumentSettings::InstrumentSettings():
		MetaInfoInterface(),
		scan_mode_(SCANMODENULL),
		polarity_(IonSource::POLNULL),
		mz_range_start_(0.0),
		mz_range_stop_(0.0)
	{
		
	}
	
	InstrumentSettings::InstrumentSettings(const InstrumentSettings& source):
		MetaInfoInterface(source),
	  scan_mode_(source.scan_mode_),
	  polarity_(source.polarity_),
	  mz_range_start_(source.mz_range_start_),
	  mz_range_stop_(source.mz_range_stop_)
	{
	  
	}
	
	InstrumentSettings::~InstrumentSettings()
	{
	  
	}
	
	InstrumentSettings& InstrumentSettings::operator = (const InstrumentSettings& source)
	{
		if (&source == this) return *this;
	  
	  scan_mode_ = source.scan_mode_;
	  polarity_ = source.polarity_;
	  mz_range_start_ = source.mz_range_start_;
	  mz_range_stop_ = source.mz_range_stop_;
	  MetaInfoInterface::operator=(source);
	  
	  return *this;
	}

  bool InstrumentSettings::operator== (const InstrumentSettings& rhs) const
  {
  	return 
	 	  scan_mode_ == rhs.scan_mode_ &&
		  polarity_ == rhs.polarity_ &&
		  mz_range_start_ == rhs.mz_range_start_ &&
		  mz_range_stop_ == rhs.mz_range_stop_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }
  
  bool InstrumentSettings::operator!= (const InstrumentSettings& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	InstrumentSettings::ScanMode InstrumentSettings::getScanMode() const 
	{
	  return scan_mode_; 
	}
	
	void InstrumentSettings::setScanMode(InstrumentSettings::ScanMode scan_mode)
	{
	  scan_mode_ = scan_mode; 
	}
	
	IonSource::Polarity InstrumentSettings::getPolarity() const 
	{
	  return polarity_; 
	}
	
	void InstrumentSettings::setPolarity(IonSource::Polarity polarity)
	{
	  polarity_ = polarity; 
	}
	
	float InstrumentSettings::getMzRangeStart() const 
	{
	  return mz_range_start_; 
	}
	
	void InstrumentSettings::setMzRangeStart(float mz_range_start)
	{
	  mz_range_start_ = mz_range_start; 
	}
	
	float InstrumentSettings::getMzRangeStop() const 
	{
	  return mz_range_stop_; 
	}
	
	void InstrumentSettings::setMzRangeStop(float mz_range_stop)
	{
	  mz_range_stop_ = mz_range_stop; 
	}
	
}

