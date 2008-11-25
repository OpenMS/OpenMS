// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/ExperimentalSettings.h>

using namespace std;

namespace OpenMS
{

	const std::string ExperimentalSettings::NamesOfNativeIDType[] = {"Unknown","Thermo","Waters","WIFF","Bruker/Agilent","Bruker BAF","Bruker FID","multiple peak lists","single peak list","scan number","spectrum identifier"};


	ExperimentalSettings::ExperimentalSettings():
		MetaInfoInterface(),
		DocumentIdentifier(),
		native_id_type_(UNKNOWN),
		sample_(),
		source_files_(),
		contacts_(),
		instrument_(),
		data_processing_(),
		hplc_(),
		datetime_(),
		comment_(),
		protein_identifications_()
	{
	}
	
	ExperimentalSettings::ExperimentalSettings(const ExperimentalSettings& source):
		MetaInfoInterface(source),
		DocumentIdentifier(source),
		native_id_type_(source.native_id_type_),
	  sample_(source.sample_),
	  source_files_(source.source_files_),
	  contacts_(source.contacts_),
	  instrument_(source.instrument_),
	  data_processing_(source.data_processing_),
	  hplc_(source.hplc_),
	  datetime_(source.datetime_),
	  comment_(source.comment_),
		protein_identifications_(source.protein_identifications_)
	{
	}
	
	ExperimentalSettings::~ExperimentalSettings()
	{
	}
	
	ExperimentalSettings& ExperimentalSettings::operator = (const ExperimentalSettings& source)
	{
	  if (&source == this) return *this;
	  
	  native_id_type_ = source.native_id_type_;
    sample_ = source.sample_;
    source_files_ = source.source_files_;
    contacts_ = source.contacts_;
    instrument_ = source.instrument_;
    data_processing_ = source.data_processing_;
    hplc_ = source.hplc_;
    datetime_ = source.datetime_;
    comment_ = source.comment_;
    protein_identifications_ = source.protein_identifications_;
    MetaInfoInterface::operator=(source);
		DocumentIdentifier::operator=(source);
	  
	  return *this;
	}

  bool ExperimentalSettings::operator== (const ExperimentalSettings& rhs) const
  {
  	return
  		native_id_type_ == rhs.native_id_type_ &&
	    sample_ == rhs.sample_ &&
	    source_files_ == rhs.source_files_ &&
	    contacts_ == rhs.contacts_ &&
	    instrument_ == rhs.instrument_ &&
	    data_processing_ == rhs.data_processing_ &&
	    hplc_ == rhs.hplc_ &&
	    datetime_ == rhs.datetime_ &&
    	protein_identifications_ == rhs.protein_identifications_ &&
    	comment_ == rhs.comment_ &&
  		MetaInfoInterface::operator==(rhs) &&
			DocumentIdentifier::operator==(rhs)
  		;
  }

  bool ExperimentalSettings::operator!= (const ExperimentalSettings& rhs) const
  {
  	return !(operator==(rhs));
 	}
 
	const Sample& ExperimentalSettings::getSample() const 
	{
	  return sample_; 
	}
	
	Sample&  ExperimentalSettings::getSample()
	{
	  return sample_; 
	}
	
	void ExperimentalSettings::setSample(const Sample& sample)
	{
	  sample_ = sample; 
	}
	
	const vector<SourceFile>& ExperimentalSettings::getSourceFiles() const 
	{
	  return source_files_; 
	}
	
	vector<SourceFile>& ExperimentalSettings::getSourceFiles()
	{
	  return source_files_; 
	}
	
	void ExperimentalSettings::setSourceFiles(const vector<SourceFile>& source_file)
	{
	  source_files_ = source_file; 
	}
	
	const vector<ContactPerson>& ExperimentalSettings::getContacts() const 
	{
	  return contacts_; 
	}
	
	vector<ContactPerson>&  ExperimentalSettings::getContacts()
	{
	  return contacts_; 
	}
	
	void ExperimentalSettings::setContacts(const std::vector<ContactPerson>& contacts)
	{
	  contacts_ = contacts; 
	}
	
	const Instrument& ExperimentalSettings::getInstrument() const 
	{
	  return instrument_; 
	}
	
	Instrument&  ExperimentalSettings::getInstrument()
	{
	  return instrument_; 
	}
	
	void ExperimentalSettings::setInstrument(const Instrument& instrument)
	{
	  instrument_ = instrument; 
	}
	
	const vector<DataProcessing>& ExperimentalSettings::getDataProcessing() const 
	{
	  return data_processing_; 
	}
	
	vector<DataProcessing>&  ExperimentalSettings::getDataProcessing()
	{
	  return data_processing_; 
	}
	
	void ExperimentalSettings::setDataProcessing(const vector<DataProcessing>& processing_method)
	{
	  data_processing_ = processing_method; 
	}

	const DateTime& ExperimentalSettings::getDateTime() const
	{
		return datetime_;
	}
	
	void ExperimentalSettings::setDateTime(const DateTime& date)
	{
		datetime_ = date;
	}


	const HPLC& ExperimentalSettings::getHPLC() const
	{
		return hplc_;
	}
	
	HPLC& ExperimentalSettings::getHPLC()
	{
		return hplc_;
	}
	
	void ExperimentalSettings::setHPLC(const HPLC& hplc)
	{
		hplc_ = hplc;
	}

	std::ostream& operator << (std::ostream& os, const ExperimentalSettings& /*exp*/)
	{
		os << "-- EXPERIMENTALSETTINGS BEGIN --"<<std::endl;
		os << "-- EXPERIMENTALSETTINGS END --"<<std::endl;
		return os;
	}
	
 	const vector<ProteinIdentification>& ExperimentalSettings::getProteinIdentifications() const
 	{
  	return protein_identifications_;	   		
 	}	
 		    	
  vector<ProteinIdentification>& ExperimentalSettings::getProteinIdentifications()
  {
  	return protein_identifications_;	
  }
  
  void ExperimentalSettings::setProteinIdentifications(const vector<ProteinIdentification>& protein_identifications)
  {
  	protein_identifications_ = protein_identifications;
  }
  
  void ExperimentalSettings::addProteinIdentification(ProteinIdentification& protein_identification)
  {
  	protein_identifications_.push_back(protein_identification);
  }

	const String& ExperimentalSettings::getComment() const 
	{
	  return comment_; 
	}
	
	void ExperimentalSettings::setComment(const String& comment)
	{
	  comment_ = comment; 
	}

  ExperimentalSettings::NativeIDType ExperimentalSettings::getNativeIDType() const
	{
		return native_id_type_;
	}
	
  void ExperimentalSettings::setNativeIDType(ExperimentalSettings::NativeIDType type)
	{
		native_id_type_ = type;
	}

}

