// -*- mode: C++; tab-width: 2; -*-
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ChromatogramSettings.h>


using namespace std;

namespace OpenMS
{

	ChromatogramSettings::ChromatogramSettings()
		: MetaInfoInterface(),
			native_id_(),
			comment_(),
			instrument_settings_(),
			source_file_(),
			acquisition_info_(),
			precursor_(),
			product_(),
			data_processing_(),
			type_()
	{
	}

	ChromatogramSettings::ChromatogramSettings(const ChromatogramSettings& source)
		: MetaInfoInterface(source),
			native_id_(source.native_id_),
			comment_(source.comment_),
			instrument_settings_(source.instrument_settings_),
			source_file_(source.source_file_),
			acquisition_info_(source.acquisition_info_),
			precursor_(source.precursor_),
			product_(source.product_),
			data_processing_(source.data_processing_),
			type_(source.type_)
	{
	}
	
	ChromatogramSettings::~ChromatogramSettings()
	{	  
	}
	
	ChromatogramSettings& ChromatogramSettings::operator = (const ChromatogramSettings& source)
	{
	  if (&source == this) return *this;
	  
	  MetaInfoInterface::operator=(source);
	  native_id_ = source.native_id_;
    comment_ = source.comment_;
    instrument_settings_ = source.instrument_settings_;
    acquisition_info_ = source.acquisition_info_;
    source_file_ = source.source_file_;
    precursor_ = source.precursor_;
	  product_ = source.product_;
	  data_processing_ = source.data_processing_;
		type_ = source.type_;
		
	  return *this;
	}

  bool ChromatogramSettings::operator== (const ChromatogramSettings& rhs) const
  {
  	return
	  	MetaInfoInterface::operator==(rhs) &&
	 		native_id_ == rhs.native_id_ &&
	    comment_ == rhs.comment_ &&
	    instrument_settings_ == rhs.instrument_settings_ &&
	    acquisition_info_ == rhs.acquisition_info_ &&
		  source_file_ == rhs.source_file_ &&
	    precursor_ == rhs.precursor_ &&
	    product_ == rhs.product_ &&
			data_processing_ == rhs.data_processing_ &&
			type_ == rhs.type_
  		;
  }
  
  bool ChromatogramSettings::operator!= (const ChromatogramSettings& rhs) const
  {
  	return !(operator==(rhs));
 	}

	const String& ChromatogramSettings::getComment() const 
	{
	  return comment_; 
	}
	
	void ChromatogramSettings::setComment(const String& comment)
	{
	  comment_ = comment; 
	}
	
	const InstrumentSettings& ChromatogramSettings::getInstrumentSettings() const 
	{
	  return instrument_settings_; 
	}
	
	InstrumentSettings&  ChromatogramSettings::getInstrumentSettings()
	{
	  return instrument_settings_; 
	}
	
	void ChromatogramSettings::setInstrumentSettings(const InstrumentSettings& instrument_settings)
	{
	  instrument_settings_ = instrument_settings; 
	}
	
	const AcquisitionInfo& ChromatogramSettings::getAcquisitionInfo() const 
	{
	  return acquisition_info_; 
	}
	
	AcquisitionInfo&  ChromatogramSettings::getAcquisitionInfo()
	{
	  return acquisition_info_; 
	}
	
	void ChromatogramSettings::setAcquisitionInfo(const AcquisitionInfo& acquisition_info)
	{
	  acquisition_info_ = acquisition_info; 
	}
	
	const SourceFile& ChromatogramSettings::getSourceFile() const 
	{
	  return source_file_; 
	}
	
	SourceFile&  ChromatogramSettings::getSourceFile()
	{
	  return source_file_; 
	}
	
	void ChromatogramSettings::setSourceFile(const SourceFile& source_file)
	{
	  source_file_ = source_file; 
	}
	
	const Precursor& ChromatogramSettings::getPrecursor() const 
	{
	  return precursor_; 
	}
	
	Precursor&  ChromatogramSettings::getPrecursor()
	{
	  return precursor_; 
	}
	
	void ChromatogramSettings::setPrecursor(const Precursor& precursor)
	{
	  precursor_ = precursor; 
	}

	const Product& ChromatogramSettings::getProduct() const 
	{
	  return product_; 
	}
	
	Product&  ChromatogramSettings::getProduct()
	{
	  return product_; 
	}
	
	void ChromatogramSettings::setProduct(const Product& product)
	{
	  product_ = product; 
	}

	std::ostream& operator << (std::ostream& os, const ChromatogramSettings& /*spec*/)
	{
		os << "-- CHROMATOGRAMSETTINGS BEGIN --"<<std::endl;
		os << "-- CHROMATOGRAMSETTINGS END --"<<std::endl;		
		return os;
	}
	
  const String& ChromatogramSettings::getNativeID() const
	{
		return native_id_;
	}
	
  void ChromatogramSettings::setNativeID(const String& native_id)
	{
		native_id_ = native_id;
	}

	const vector<DataProcessing>& ChromatogramSettings::getDataProcessing() const 
	{
	  return data_processing_; 
	}
	
	vector<DataProcessing>&  ChromatogramSettings::getDataProcessing()
	{
	  return data_processing_; 
	}
	
	void ChromatogramSettings::setDataProcessing(const vector<DataProcessing>& processing_method)
	{
	  data_processing_ = processing_method; 
	}

	ChromatogramSettings::ChromatogramType ChromatogramSettings::getChromatogramType() const
	{
		return type_;
	}

	void ChromatogramSettings::setChromatogramType(ChromatogramType type)
	{
		type_ = type;
	}
}




