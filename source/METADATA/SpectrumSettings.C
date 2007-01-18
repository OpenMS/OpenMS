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

#include <OpenMS/METADATA/SpectrumSettings.h>


using namespace std;

namespace OpenMS
{

	const std::string SpectrumSettings::NamesOfSpectrumType[] = {"Unknown","Peak data","Raw data"};

	SpectrumSettings::SpectrumSettings():
		type_(UNKNOWN), 
		identification_()
	{
	  
	}
	
	SpectrumSettings::SpectrumSettings(const SpectrumSettings& source):
		type_(source.type_),
	  comment_(source.comment_),
	  instrument_settings_(source.instrument_settings_),
	  source_file_(source.source_file_),
	  acquisition_info_(source.acquisition_info_),
	  meta_info_descriptions_(source.meta_info_descriptions_),
	  precursor_(source.precursor_),
	  identification_(source.identification_)
	{
	  
	}
	
	SpectrumSettings::~SpectrumSettings()
	{
	  
	}
	
	SpectrumSettings& SpectrumSettings::operator = (const SpectrumSettings& source)
	{
	  if (&source == this) return *this;
	  
	  type_ = source.type_;
    comment_ = source.comment_;
    instrument_settings_ = source.instrument_settings_;
    acquisition_info_ = source.acquisition_info_;
    source_file_ = source.source_file_;
    meta_info_descriptions_ = source.meta_info_descriptions_;
    precursor_ = source.precursor_;
    identification_ = source.identification_;
	  
	  return *this;
	}

  bool SpectrumSettings::operator== (const SpectrumSettings& rhs) const
  {
  	return
	   	type_ == rhs.type_ &&
	    comment_ == rhs.comment_ &&
	    instrument_settings_ == rhs.instrument_settings_ &&
	    acquisition_info_ == rhs.acquisition_info_ &&
		  source_file_ == rhs.source_file_ &&
	    meta_info_descriptions_ == rhs.meta_info_descriptions_ &&
	    precursor_ == rhs.precursor_ &&
	    identification_ == rhs.identification_
  		;
  }
  
  bool SpectrumSettings::operator!= (const SpectrumSettings& rhs) const
  {
  	return !(operator==(rhs));
 	}

	SpectrumSettings::SpectrumType SpectrumSettings::getType() const
	{
		return type_;	
	}

	void SpectrumSettings::setType(SpectrumSettings::SpectrumType type)
	{
		type_ = type;
	}
	
	const String& SpectrumSettings::getComment() const 
	{
	  return comment_; 
	}
	
	void SpectrumSettings::setComment(const String& comment)
	{
	  comment_ = comment; 
	}
	
	const InstrumentSettings& SpectrumSettings::getInstrumentSettings() const 
	{
	  return instrument_settings_; 
	}
	
	InstrumentSettings&  SpectrumSettings::getInstrumentSettings()
	{
	  return instrument_settings_; 
	}
	
	void SpectrumSettings::setInstrumentSettings(const InstrumentSettings& instrument_settings)
	{
	  instrument_settings_ = instrument_settings; 
	}
	
	const AcquisitionInfo& SpectrumSettings::getAcquisitionInfo() const 
	{
	  return acquisition_info_; 
	}
	
	AcquisitionInfo&  SpectrumSettings::getAcquisitionInfo()
	{
	  return acquisition_info_; 
	}
	
	void SpectrumSettings::setAcquisitionInfo(const AcquisitionInfo& acquisition_info)
	{
	  acquisition_info_ = acquisition_info; 
	}
	
	const SourceFile& SpectrumSettings::getSourceFile() const 
	{
	  return source_file_; 
	}
	
	SourceFile&  SpectrumSettings::getSourceFile()
	{
	  return source_file_; 
	}
	
	void SpectrumSettings::setSourceFile(const SourceFile& source_file)
	{
	  source_file_ = source_file; 
	}
	
	const std::map<String,MetaInfoDescription>& SpectrumSettings::getMetaInfoDescriptions() const 
	{
	  return meta_info_descriptions_; 
	}
	
	std::map<String,MetaInfoDescription>&  SpectrumSettings::getMetaInfoDescriptions()
	{
	  return meta_info_descriptions_; 
	}
	
	void SpectrumSettings::setMetaInfoDescriptions(const std::map<String,MetaInfoDescription>& meta_info_descriptions)
	{
	  meta_info_descriptions_ = meta_info_descriptions; 
	}
	
	const Precursor& SpectrumSettings::getPrecursor() const 
	{
	  return precursor_; 
	}
	
	Precursor&  SpectrumSettings::getPrecursor()
	{
	  return precursor_; 
	}
	
	void SpectrumSettings::setPrecursor(const Precursor& precursor)
	{
	  precursor_ = precursor; 
	}
	
	std::ostream& operator << (std::ostream& os, const SpectrumSettings& /*spec*/)
	{
		os << "-- SPECTRUMSETTINGS BEGIN --"<<std::endl;
		os << "-- SPECTRUMSETTINGS END --"<<std::endl;		
		return os;
	}
	
  const std::vector<Identification>& SpectrumSettings::getIdentifications() const 
  {
  	return identification_;
  }

  std::vector<Identification>& SpectrumSettings::getIdentifications() 
  {
  	return identification_;
  }

	void SpectrumSettings::setIdentifications(const std::vector<Identification>& identification)
	{
		identification_ = identification;
	}	
}




