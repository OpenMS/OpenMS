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

#include <OpenMS/METADATA/ExperimentalSettings.h>

using namespace std;

namespace OpenMS
{

	const std::string ExperimentalSettings::NamesOfExperimentType[] = {"Unknown","MS","MS_MS","HPLC_MS","HPLC_MS_MS"};

	ExperimentalSettings::ExperimentalSettings():
		MetaInfoInterface(),
		sample_(),
		source_file_(),
		contacts_(),
		instrument_(),
		software_(),
		processing_method_(),
		hplc_(),
		type_(UNKNOWN),
		date_(),
		comment_(),
		protein_identifications_()
	{
	  
	}
	
	ExperimentalSettings::ExperimentalSettings(const ExperimentalSettings& source):
		MetaInfoInterface(source),
	  sample_(source.sample_),
	  source_file_(source.source_file_),
	  contacts_(source.contacts_),
	  instrument_(source.instrument_),
	  software_(source.software_),
	  processing_method_(source.processing_method_),
	  hplc_(source.hplc_),
	  type_(source.type_),
	  date_(source.date_),
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
	  
    sample_ = source.sample_;
    source_file_ = source.source_file_;
    contacts_ = source.contacts_;
    instrument_ = source.instrument_;
    software_ = source.software_;
    processing_method_ = source.processing_method_;
    hplc_ = source.hplc_;
    type_ = source.type_;
    date_ = source.date_;
    comment_ = source.comment_;
    protein_identifications_ = source.protein_identifications_;
    MetaInfoInterface::operator=(source);
	  
	  return *this;
	}

  bool ExperimentalSettings::operator== (const ExperimentalSettings& rhs) const
  {
  	return  
	    sample_ == rhs.sample_ &&
	    source_file_ == rhs.source_file_ &&
	    contacts_ == rhs.contacts_ &&
	    instrument_ == rhs.instrument_ &&
	    software_ == rhs.software_ &&
	    processing_method_ == rhs.processing_method_ &&
	    hplc_ == rhs.hplc_ &&
	    type_ == rhs.type_ &&
	    date_ == rhs.date_ &&
    	protein_identifications_ == rhs.protein_identifications_ &&
    	comment_ == rhs.comment_ &&
  		MetaInfoInterface::operator==(rhs)
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
	
	const SourceFile& ExperimentalSettings::getSourceFile() const 
	{
	  return source_file_; 
	}
	
	SourceFile&  ExperimentalSettings::getSourceFile()
	{
	  return source_file_; 
	}
	
	void ExperimentalSettings::setSourceFile(const SourceFile& source_file)
	{
	  source_file_ = source_file; 
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
	
	const Software& ExperimentalSettings::getSoftware() const 
	{
	  return software_; 
	}
	
	Software&  ExperimentalSettings::getSoftware()
	{
	  return software_; 
	}
	
	void ExperimentalSettings::setSoftware(const Software& software)
	{
	  software_ = software; 
	}
	
	const ProcessingMethod& ExperimentalSettings::getProcessingMethod() const 
	{
	  return processing_method_; 
	}
	
	ProcessingMethod&  ExperimentalSettings::getProcessingMethod()
	{
	  return processing_method_; 
	}
	
	void ExperimentalSettings::setProcessingMethod(const ProcessingMethod& processing_method)
	{
	  processing_method_ = processing_method; 
	}

	ExperimentalSettings::ExperimentType ExperimentalSettings::getType() const
	{
		return type_;
	}
	
	void ExperimentalSettings::setType(ExperimentalSettings::ExperimentType type)
	{
		type_ = type;
	}

	const Date& ExperimentalSettings::getDate() const
	{
		return date_;
	}
	
	void ExperimentalSettings::setDate(const Date& date)
	{
		date_ = date;
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
}

