// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/DataProcessing.h>


using namespace std;

namespace OpenMS
{
	const std::string DataProcessing::NamesOfProcessingAction[] = {"Data processing action",
																																 "Charge deconvolution",
																																 "Deisotoping",
																																 "Smoothing",
																																 "Charge calculation",
																																 "Precursor recalculation",
																																 "Baseline reduction",
																																 "Peak picking",
																																 "Retention time alignment",
																																 "Calibration of m/z positions",
																																 "Intensity normalization",
																																 "Data filtering",
																																 "Quantitation",
																																 "Feature grouping",
																																 "Identification mapping",
																																 "File format conversion",
																																 "Conversion to mzData format",
																																 "Conversion to mzML format",
																																 "Conversion to mzXML format",
																																 "Conversion to DTA format"
																																};

	DataProcessing::DataProcessing():
		MetaInfoInterface(),
		software_(),
		processing_actions_(),
		completion_time_()
	{
		
	}
	
	DataProcessing::DataProcessing(const DataProcessing& rhs):
		MetaInfoInterface(rhs),
		software_(rhs.software_),
		processing_actions_(rhs.processing_actions_),
		completion_time_(rhs.completion_time_)
	{
	  
	}
	
	DataProcessing::~DataProcessing()
	{
	  
	}
	
	DataProcessing& DataProcessing::operator = (const DataProcessing& rhs)
	{
	  if (&rhs == this) return *this;
	  
	  MetaInfoInterface::operator=(rhs);
		software_ = rhs.software_;
		processing_actions_ = rhs.processing_actions_;
		completion_time_ = rhs.completion_time_;
	  
	  return *this;
	}
	
	bool DataProcessing::operator== (const DataProcessing& rhs) const
	{
			return 
				software_ == rhs.software_ &&
				processing_actions_ == rhs.processing_actions_ &&
				completion_time_ == rhs.completion_time_ &&
				MetaInfoInterface::operator==(rhs)
				;
	}
	
	bool DataProcessing::operator!= (const DataProcessing& rhs) const
	{
		return !(operator==(rhs));
	}


	const Software& DataProcessing::getSoftware() const 
	{
	  return software_; 
	}
	
	Software& DataProcessing::getSoftware()
	{
	  return software_; 
	}
	
	void DataProcessing::setSoftware(const Software& software)
	{
	  software_ = software; 
	}	

	const DateTime& DataProcessing::getCompletionTime() const
	{
	  return completion_time_;
	}

	void DataProcessing::setCompletionTime(const DateTime& completion_time)
	{
	  completion_time_ = completion_time;
	}

  const set<DataProcessing::ProcessingAction>& DataProcessing::getProcessingActions() const
	{
		return processing_actions_;
	}
	
  set<DataProcessing::ProcessingAction>& DataProcessing::getProcessingActions()
	{
		return processing_actions_;
	}
	
  void DataProcessing::setProcessingActions(const set<DataProcessing::ProcessingAction>& processing_actions)
	{
		processing_actions_ = processing_actions;
	}

}

