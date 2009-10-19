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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumSettings.h>


using namespace std;

namespace OpenMS
{

	const std::string SpectrumSettings::NamesOfSpectrumType[] = {"Unknown","Peak data","Raw data"};

	SpectrumSettings::SpectrumSettings()
		: MetaInfoInterface(),
			type_(UNKNOWN),
			native_id_(),
			comment_(),
			instrument_settings_(),
			source_file_(),
			acquisition_info_(),
			precursors_(),
			products_(),
			identification_(),
			data_processing_()
	{
	}

	SpectrumSettings::SpectrumSettings(const SpectrumSettings& source)
		: MetaInfoInterface(source),
			type_(source.type_),
			native_id_(source.native_id_),
			comment_(source.comment_),
			instrument_settings_(source.instrument_settings_),
			source_file_(source.source_file_),
			acquisition_info_(source.acquisition_info_),
			precursors_(source.precursors_),
			products_(source.products_),
			identification_(source.identification_),
			data_processing_(source.data_processing_)
	{
	}
	
	SpectrumSettings::~SpectrumSettings()
	{	  
	}
	
	SpectrumSettings& SpectrumSettings::operator = (const SpectrumSettings& source)
	{
	  if (&source == this) return *this;
	  
	  MetaInfoInterface::operator=(source);
	  type_ = source.type_;
	  native_id_ = source.native_id_;
    comment_ = source.comment_;
    instrument_settings_ = source.instrument_settings_;
    acquisition_info_ = source.acquisition_info_;
    source_file_ = source.source_file_;
    precursors_ = source.precursors_;
    identification_ = source.identification_;
	  products_ = source.products_;
	  data_processing_ = source.data_processing_;
		
	  return *this;
	}

  bool SpectrumSettings::operator== (const SpectrumSettings& rhs) const
  {
  	return
	  	MetaInfoInterface::operator==(rhs) &&
	   	type_ == rhs.type_ &&
	 		native_id_ == rhs.native_id_ &&
	    comment_ == rhs.comment_ &&
	    instrument_settings_ == rhs.instrument_settings_ &&
	    acquisition_info_ == rhs.acquisition_info_ &&
		  source_file_ == rhs.source_file_ &&
	    precursors_ == rhs.precursors_ &&
	    identification_ == rhs.identification_ &&
	    products_ == rhs.products_ &&
			data_processing_ == rhs.data_processing_
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
	
	const vector<Precursor>& SpectrumSettings::getPrecursors() const 
	{
	  return precursors_; 
	}
	
	vector<Precursor>&  SpectrumSettings::getPrecursors()
	{
	  return precursors_; 
	}
	
	void SpectrumSettings::setPrecursors(const vector<Precursor>& precursors)
	{
	  precursors_ = precursors; 
	}

	const vector<Product>& SpectrumSettings::getProducts() const 
	{
	  return products_; 
	}
	
	vector<Product>&  SpectrumSettings::getProducts()
	{
	  return products_; 
	}
	
	void SpectrumSettings::setProducts(const vector<Product>& products)
	{
	  products_ = products; 
	}

	std::ostream& operator << (std::ostream& os, const SpectrumSettings& /*spec*/)
	{
		os << "-- SPECTRUMSETTINGS BEGIN --"<<std::endl;
		os << "-- SPECTRUMSETTINGS END --"<<std::endl;		
		return os;
	}
	
  const std::vector<PeptideIdentification>& SpectrumSettings::getPeptideIdentifications() const 
  {
  	return identification_;
  }

  std::vector<PeptideIdentification>& SpectrumSettings::getPeptideIdentifications() 
  {
  	return identification_;
  }

	void SpectrumSettings::setPeptideIdentifications(const std::vector<PeptideIdentification>& identification)
	{
		identification_ = identification;
	}

  const String& SpectrumSettings::getNativeID() const
	{
		return native_id_;
	}
	
  void SpectrumSettings::setNativeID(const String& native_id)
	{
		native_id_ = native_id;
	}

	const vector<DataProcessing>& SpectrumSettings::getDataProcessing() const 
	{
	  return data_processing_; 
	}
	
	vector<DataProcessing>&  SpectrumSettings::getDataProcessing()
	{
	  return data_processing_; 
	}
	
	void SpectrumSettings::setDataProcessing(const vector<DataProcessing>& processing_method)
	{
	  data_processing_ = processing_method; 
	}
	
}




