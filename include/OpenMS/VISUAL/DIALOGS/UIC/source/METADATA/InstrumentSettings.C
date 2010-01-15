// -*- mode: C++; tab-width: 2; -*-s
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/InstrumentSettings.h>

using namespace std;

namespace OpenMS
{
	const std::string InstrumentSettings::NamesOfScanMode[] = {"Unknown","MassSpectrum","SelectedIonMonitoring","SelectedReactionMonitoring","ConsecutiveReactionMonitoring","ConstantNeutralGain","ConstantNeutralLoss","Precursor","EnhancedMultiplyCharged","TimeDelayedFragmentation","ElectromagneticRadiation","Emission","Absorbtion"};

	InstrumentSettings::InstrumentSettings():
		MetaInfoInterface(),
		scan_mode_(UNKNOWN),
		zoom_scan_(false),
		polarity_(IonSource::POLNULL),
		scan_windows_()
	{
	}
	
	InstrumentSettings::InstrumentSettings(const InstrumentSettings& source):
		MetaInfoInterface(source),
	  scan_mode_(source.scan_mode_),
		zoom_scan_(source.zoom_scan_),
	  polarity_(source.polarity_),
		scan_windows_(source.scan_windows_)
	{
	}
	
	InstrumentSettings::~InstrumentSettings()
	{
	}
	
	InstrumentSettings& InstrumentSettings::operator = (const InstrumentSettings& source)
	{
		if (&source == this) return *this;
	  
	  scan_mode_ = source.scan_mode_;
		zoom_scan_  = source.zoom_scan_;
	  polarity_ = source.polarity_;
	  scan_windows_ = source.scan_windows_;
	  MetaInfoInterface::operator=(source);
	  
	  return *this;
	}

  bool InstrumentSettings::operator== (const InstrumentSettings& rhs) const
  {
  	return 
	 	  scan_mode_ == rhs.scan_mode_ &&
			zoom_scan_  == rhs.zoom_scan_ &&
		  polarity_ == rhs.polarity_ &&
		  scan_windows_ == rhs.scan_windows_ &&
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
	
	const std::vector< ScanWindow >&  InstrumentSettings::getScanWindows() const
	{
	  return scan_windows_;
	}
	
	std::vector< ScanWindow >&  InstrumentSettings::getScanWindows()
	{
	  return scan_windows_;
	}
	
	void InstrumentSettings::setScanWindows(std::vector< ScanWindow >  scan_windows)
	{
	  scan_windows_ =  scan_windows;
	}

	bool InstrumentSettings::getZoomScan() const
	{
		return zoom_scan_;
	}
	
	void InstrumentSettings::setZoomScan(bool zoom_scan)
	{
		zoom_scan_ = zoom_scan;
	}

}

