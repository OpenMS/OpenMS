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

#ifndef OPENMS_FORMAT_HANDLERS_MZXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZXMLHANDLER_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <stack>

namespace OpenMS
{
	class MetaInfoInterface;
	
	namespace Internal
	{
  /**
  	@brief XML handlers for MzXMLFile

		MapType has to be a MSExperiment or have the same interface.
		Do not use this class. It is only needed in MzXMLFile.
  */
	template <typename MapType>
  class MzXMLHandler
		: public SchemaHandler
  {
    public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzXMLHandler(MapType& exp, const String& filename, ProgressLogger& logger)
			: SchemaHandler(TAG_NUM,MAP_NUM,filename), // number of tags, number of maps
				exp_(&exp),	
				cexp_(0),
				spec_(0),
				analyzer_(0),
				decoder_(),
				peak_count_(0),
				char_rest_(),
				last_scan_num_(0),
				spec_write_counter_(1),
				logger_(logger)
  		{
				fillMaps_(Schemes::MzXML[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}

      /// Constructor for a read-only handler
      MzXMLHandler(const MapType& exp, const String& filename, const ProgressLogger& logger)
			: SchemaHandler(TAG_NUM,MAP_NUM,filename), // number of tags, number of maps
				exp_(0), 
				cexp_(&exp),
				spec_(),
				decoder_(),
				peak_count_(0),
				char_rest_(),
				last_scan_num_(0),
				spec_write_counter_(1),
				logger_(logger)
  		{
				fillMaps_(Schemes::MzXML[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}

  		/// Destructor
      virtual ~MzXMLHandler(){}
      //@}
			
			// Docu in base class
      virtual void endElement( const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, const unsigned int length);

  		///Write the contents to a stream
			void writeTo(std::ostream& os);
			
			void setOptions(const PeakFileOptions& opt) { options_ = opt; }

    protected:
			/// map pointer for reading
			MapType* exp_;
			/// map pointer for writing
			const MapType* cexp_;
	
			/** @brief indices for tags used by mzData
	
				Used to access is_parser_in_tag_.
				If you add tags, also add them to XMLSchemes.h.
				Add no elements to the enum after TAG_NUM.
			*/
			enum Tags { TAGNULL=0, MSRUN, INDEX, OFFSET, SHA1, PARENTFILE, INSTRUMENT, DATAPROCESSING,
									SEPARATION, SPOTTING, SCAN, SCANORIGIN, PRECURSORMZ, MALDI, PEAKS, NAMEVALUE,
									COMMENT, SOFTWARE, INDEXOFFSET, OPERATOR, MANUFACTURER, MODEL, IONISATION,
									ANALYZER, DETECTOR, RESOLUTION, MZXML, PROCESSING, SEPARATIONTECH, TAG_NUM};
	
	
			/** @brief indices for attributes used by MzXML
	
				Used to access enum2str_() with ATTMAP.
				If you add terms, also add them to XMLSchemes.h
			*/
			enum Attributes {ATTNULL, POLARITY, SCANTYPE, CENTROIDED, DEISOTOPED,
											 DECONVOLUTED, RETTIME,IONENERGY, COLLENERGY, PRESSURE,
											 STARTMZ, ENDMZ, LOWMZ, HIGHMZ, BASEPEAKMZ, BASEPEAKINT,
											 TOTIONCURRENT, PEAKSCOUNT, NUM, MSLEVEL, SCANCOUNT,
											 FILENAME, FILETYPE, SOFTWAREVERSION, NAME, TYPE,
											 COMPLETION_TIME, PRECURSOR_INTENSITY, PRECURSOR_CHARGE,
											 FIRST_NAME, LAST_NAME, EMAIL, PHONE, URI, VALUE, CATEGORY,
											 PRECISION, BYTEORDER, PAIRORDER, SCHEMA, SPOT_INTEGRATION,
											 INTENSITY_CUTOFF, STARTTIME, ENDTIME, FILESHA1, PARENTFILEID,
											 PRECURSOR_SCANNUM, WINDOW_WIDENESS, PLATEID, SPOTID,
											 LASER_SHOOT_COUNT, LASER_FREQUENCY, LASER_INTESITY};
	
			/** @brief indices for enum2str-maps used by mzXML
	
				Used to access enum2str_().
				If you add maps, also add them to XMLSchemes.h.
				Add no elements to the enum after MAP_NUM.
				Each map corresponds to a string in XMLSchemes.h.
			*/
			enum MapTypes {	POLARITYMAP=0, IONTYPEMAP, TYPEMAP, ANALYZERTYPEMAP, SCANMODEMAP,
											ATTMAP, TAGMAP, RESMETHODMAP, PEAKPROCMAP, PRECISIONMAP,MAP_NUM};
	
			/// Possible precisions for Base64 data encoding
			enum Precision { UNKNOWN_PRECISION, REAL, DOUBLE};
	
			typedef typename MapType::SpectrumType SpectrumType;
			typedef typename SpectrumType::PeakType PeakType;
			typedef typename SpectrumType::Iterator  PeakIterator;
			typedef typename SpectrumType::PrecursorPeakType PrecursorPeakType;
	
			PeakFileOptions options_;
			
			/**@name temporary datastructures to hold parsed data */
	    //@{
			SpectrumType* spec_;
			MassAnalyzer* analyzer_;
			MetaInfoDescription* meta_;
			String meta_id_;
			Base64 decoder_;
			UInt peak_count_;
			Precision precision_;
			String char_rest_;
			UInt last_scan_num_;
			//@}
	
			/// spectrum counter (spectra without peaks are not written)
			UInt spec_write_counter_;
			
			/// Progress logging class
			const ProgressLogger& logger_;
			
			/// Add name, value and description to a given MetaInfo object
			void setAddInfo_(MetaInfoInterface& info, const String& name, const String& value, const String& description)
			{
				info.metaRegistry().registerName(name, description);
				info.setMetaValue(name,value);
			}
	
			/// write metaInfo to xml (usually in nameValue-tag)
			inline void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta,
																	int indent=4, String tag="nameValue")
			{
				std::vector<String> keys;  // Vector to hold keys to meta info
				meta.getKeys(keys);
	
				for (std::vector<String>::const_iterator it = keys.begin(); it!=keys.end(); ++it)
					if ( (*it)[0] != '#')  // internally used meta info start with '#'
				{
					String name = *it;
					os << String(indent,'\t') << "<" << tag << " name=\"";
					if (tag=="processingOperation" && name.find('#')!=std::string::npos)
					{
						std::vector<String> parts;
						name.split('#',parts);
						os << parts[0] << "\" type=\"" << parts[1];
					}
					else
					{
						os << name;
					}
					os << "\" value=\""
						 << meta.getMetaValue(*it) << "\"/>\n";
				}
			}
		
		private:
			MzXMLHandler(); // not impelmented -> private
  };







	//--------------------------------------------------------------------------------



	template <typename MapType>
  void MzXMLHandler<MapType>::characters(const XMLCh* const chars, unsigned int /*length*/)
  {
    if (skip_tag_.top()) return;
    
		char* transcoded_chars = xercesc::XMLString::transcode(chars);
  		
		if(is_parser_in_tag_[PEAKS])
		{
			//chars may be split to several chunks => concatenate them
			char_rest_ += transcoded_chars;
		}
		else if (	is_parser_in_tag_[OFFSET] ||
							is_parser_in_tag_[INDEXOFFSET] ||
							is_parser_in_tag_[SHA1])
		{
			
		}
		else if (	is_parser_in_tag_[PRECURSORMZ])
		{
			if (spec_ != 0)
			{
				spec_->getPrecursorPeak().getPosition()[0] = asFloat_(transcoded_chars);
			}
		}
		else if (	is_parser_in_tag_[COMMENT])
		{
			if (is_parser_in_tag_[INSTRUMENT])
			{
				 setAddInfo_(exp_->getInstrument(),"#Comment" , transcoded_chars, "Instrument.Comment");
			}
			else if (is_parser_in_tag_[DATAPROCESSING])
			{
				setAddInfo_(exp_->getProcessingMethod(),"#Comment", transcoded_chars,"DataProcessing.Comment");
			}
			else if (is_parser_in_tag_[SCAN])
			{
				if (spec_ != 0)
				{
					spec_->setComment( transcoded_chars );
				}
			}
			else if (String(transcoded_chars).trim()!="")
			{
				std::cerr << "Unhandled characters: \"" << transcoded_chars << "\"" << std::endl;
			}
		}
		else if (String(transcoded_chars).trim()!="")
		{
				std::cerr << "Unhandled characters: \"" << transcoded_chars << "\"" << std::endl;
		}
  	//std::cout << " -- !Chars -- " << std::endl;
		
		xercesc::XMLString::release(&transcoded_chars);
  }

	/**
		@brief MzXML parsing routine. Source: http://sashimi.sourceforge.net/schema_revision/mzXML_2.1/Doc/mzXML_2.1_tutorial.pdf
	 */
	template <typename MapType>
  void MzXMLHandler<MapType>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
  	//std::cout << " -- Start -- "<< xercesc::XMLString::transcode(qname) << " -- " << std::endl;
  	//std::cout << " skip size: " << 	skip_tag_.size();
  	//if (skip_tag_.size()) std::cout << " skip current: " << skip_tag_.top();
  	//std::cout << std::endl;
		
		int tag = enterTag(qname, attributes);
		if (skip_tag_.top()) return;
		
		String tmp_str;
		switch(tag)
		{
			case MSRUN:
				tmp_str = getAttributeAsString_(SCANCOUNT, false, qname);
				if (tmp_str!="") 
				{
					exp_->reserve( asUInt_(tmp_str) );
				}
				logger_.startProgress(0,asUInt_(tmp_str),"loading mzXML file");
				// fall through to next tag because of different MzXML versions... -> no break
			case MZXML:
				// look for schema information
				if (atts_->getIndex(xercesc::XMLString::transcode(enum2str_(ATTMAP,SCHEMA).c_str()))!=-1)
				{
					tmp_str = getAttributeAsString_(SCHEMA, false, qname);
					//std::cout << "SCHEMA: " << tmp_str << std::endl;
					if (tmp_str!="")
					{
						for (UInt index=0; index<Schemes::MzXML_num; ++index)
						{
							if (tmp_str.hasSubstring(Schemes::MzXML[index][0]))
							{
								schema_ = index;
								// refill maps with older schema
								for (UInt i=0; i<str2enum_array_.size(); i++)	str2enum_array_[i].clear();
								for (UInt i=0; i<enum2str_array_.size(); i++)	enum2str_array_[i].clear();
								fillMaps_(Schemes::MzXML[schema_]);
								break;
							}
						}
					}
				}
				break;
			case PARENTFILE:
				tmp_str = getAttributeAsString_(FILENAME, true, qname);
				if (tmp_str != "") 
				{
					exp_->getSourceFile().setNameOfFile( tmp_str.c_str() );
				}
				
				tmp_str = getAttributeAsString_(FILETYPE, true, qname);
				if (tmp_str != "") 
				{
					exp_->getSourceFile().setFileType( tmp_str );
				}
				
				tmp_str = getAttributeAsString_(FILESHA1, true, qname);
				if (tmp_str != "") 
				{
					exp_->getSourceFile().setSha1(tmp_str);
				}
				break;
			case INSTRUMENT:
				{
					if (attributes.getLength()==0) break;  // attributes only in mzXML 1.0
					
					exp_->getInstrument().setModel( xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode(enum2str_(TAGMAP,MODEL).c_str()))));
					exp_->getInstrument().setVendor( xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode(enum2str_(TAGMAP,MANUFACTURER).c_str()))));
					
					MassAnalyzer analyzer;
					String str = enum2str_(TAGMAP,ANALYZER);
					analyzer.setType((MassAnalyzer::AnalyzerType)str2enum_(ANALYZERTYPEMAP,xercesc::XMLString::transcode(atts_->getValue(xercesc::XMLString::transcode(str.c_str()))),str.c_str()));
					exp_->getInstrument().getMassAnalyzers().push_back(analyzer);
					str = enum2str_(TAGMAP,IONISATION);
					exp_->getInstrument().getIonSource().setIonizationMethod((IonSource::IonizationMethod)str2enum_(IONTYPEMAP,xercesc::XMLString::transcode(atts_->getValue(xercesc::XMLString::transcode(str.c_str()))),str.c_str()));
				}
				break;
			case SOFTWARE:
				if (is_parser_in_tag_[DATAPROCESSING])
				{
					tmp_str = getAttributeAsString_(SOFTWAREVERSION, true, qname);
					if (tmp_str != "") 
					{
						exp_->getSoftware().setVersion(tmp_str);
					}
					
					tmp_str = getAttributeAsString_(NAME, true, qname);
					if (tmp_str != "") 
					{
						exp_->getSoftware().setName(tmp_str);
					}
					
					tmp_str = getAttributeAsString_(TYPE, true, qname);
					if (tmp_str != "") 
					{
						exp_->getSoftware().setComment(tmp_str);
					}
					
					tmp_str = getAttributeAsString_(COMPLETION_TIME, false, qname);
					if (tmp_str != "") 
					{
						exp_->getSoftware().setCompletionTime( asDateTime_(tmp_str) );
					}
				}
				else if (is_parser_in_tag_[INSTRUMENT])
				{
					// not part of METADATA -> putting it into MetaInfo
					std::string swn = "#InstSoftware", swn_d = "Instrument software name",
						swv = "#InstSoftwareVersion", swv_d = "Instrument software version",
						swt = "#InstSoftwareType", swt_d = "Instrument software type",
						cmpl = "#InstSoftwareTime", cmpl_d = "Instrument software completion time";
					MetaInfoRegistry& registry =	MetaInfo().registry();
					registry.registerName(swn,swn_d);
					registry.registerName(swv,swv_d);
					registry.registerName(swt,swt_d);
					registry.registerName(cmpl,cmpl_d);
					
					tmp_str = getAttributeAsString_(NAME, true, qname);
					if (tmp_str != "") 
					{
						exp_->getInstrument().setMetaValue(swn, tmp_str);
					}
					
					tmp_str = getAttributeAsString_(SOFTWAREVERSION, true, qname);
					if (tmp_str != "") 
					{
						exp_->getInstrument().setMetaValue(swv, tmp_str);
					}
					
					tmp_str = getAttributeAsString_(TYPE, true, qname);
					if (tmp_str != "") 
					{
						exp_->getInstrument().setMetaValue(swt, tmp_str);
					}
					
					tmp_str = getAttributeAsString_(COMPLETION_TIME, false, qname);
					if  (tmp_str != "") 
					{
						DateTime time = asDateTime_(tmp_str);
						time.get(tmp_str);
						exp_->getInstrument().setMetaValue(cmpl,tmp_str);
					}
				}
				break;
			case PEAKS:
				{
					checkAttribute_(PRECISION,enum2str_(PRECISIONMAP,REAL),
																		enum2str_(PRECISIONMAP,DOUBLE));
					const String str = enum2str_(ATTMAP,PRECISION);
					precision_ = (Precision) str2enum_(PRECISIONMAP,xercesc::XMLString::transcode(atts_->getValue(xercesc::XMLString::transcode(str.c_str()))),str.c_str());
					checkAttribute_(BYTEORDER,"network");
					checkAttribute_(PAIRORDER,"m/z-int");
				}
				break;
			case PRECURSORMZ:
				{
					if (!spec_) break;
					
					PrecursorPeakType& peak = spec_->getPrecursorPeak();
					
					tmp_str = getAttributeAsString_(PRECURSOR_INTENSITY, false, qname);
					if (tmp_str != "") 
					{
						peak.setIntensity( asFloat_(tmp_str) );
					}
					
					tmp_str = getAttributeAsString_(PRECURSOR_CHARGE, false, qname);
					if (tmp_str != "")
					{
						peak.setCharge(asInt_(tmp_str));
					}
					
					tmp_str = getAttributeAsString_(PRECURSOR_SCANNUM, false, qname);
					if (tmp_str != "")
					{
						// ignore
					}
					
					tmp_str = getAttributeAsString_(WINDOW_WIDENESS, false, qname);
					if (tmp_str != "")
					{
						spec_->getPrecursor().setWindowSize(asDouble_(tmp_str));
					}
				}
				break;
			case MALDI:
				{
					tmp_str = getAttributeAsString_(PLATEID, true, qname);
					if (tmp_str != "") 
					{
						// ignored
					}
					
					tmp_str = getAttributeAsString_(SPOTID, true, qname);
					if (tmp_str != "") 
					{
						// ignored
					}
					
					tmp_str = getAttributeAsString_(LASER_SHOOT_COUNT, false, qname);
					if (tmp_str != "") 
					{
						// ignored
					}
					
					tmp_str = getAttributeAsString_(LASER_FREQUENCY, false, qname);
					if (tmp_str != "") 
					{
						// ignored
					}
					
					tmp_str = getAttributeAsString_(LASER_INTESITY, false, qname);
					if (tmp_str != "") 
					{
						// ignored
					}
				}
				break;
			case SCAN:
				{
					if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					tmp_str = getAttributeAsString_(NUM, true, qname);
					
					//std::cout << "Last: " << last_scan_num_ << "  - Scan num: " << tmp_str << std::endl; 
					
					if (asUInt_(tmp_str) != last_scan_num_ + 1)
					{
						// num tag starts from 1 and must be consecutive
						error("non-consecutive numbers in 'scan' tags");
					}
					++last_scan_num_;
					
					SpectrumType spec;
					
					tmp_str = getAttributeAsString_(MSLEVEL, true, qname);
					if (tmp_str != "")
					{
						spec.setMSLevel(asInt_(tmp_str));
					}
					else
					{
						break;
					}
					
					tmp_str = getAttributeAsString_(PEAKSCOUNT, true, qname);
					if (tmp_str != "")
					{
						peak_count_ = asInt_(tmp_str);
					}
					
					//optinal attributes
					for (UInt i=0; i<attributes.getLength(); i++)
					{
						int att = str2enum_(ATTMAP,xercesc::XMLString::transcode(attributes.getQName(i)),"scan attribute");
						String value = xercesc::XMLString::transcode(attributes.getValue(i));
						InstrumentSettings& sett = spec.getInstrumentSettings();
						switch (att)
							{
							case POLARITY:
								sett.setPolarity( (IonSource::Polarity) str2enum_(POLARITYMAP,value,"polarity") );
								break;
							case SCANTYPE:
								sett.setScanMode( (InstrumentSettings::ScanMode) str2enum_(SCANMODEMAP,value,"scan mode") );
								break;
							case RETTIME:
								value.trim();
								spec.setRT( asFloat_(value.substr(2,value.size()-3)));
								//std::cout << spec_->getRT() << std::endl;
								break;
							case STARTMZ:
								sett.setMzRangeStart( asDouble_(value));
								break;
							case ENDMZ:
								sett.setMzRangeStop( asDouble_(value));
								break;
							case DEISOTOPED:
							  exp_->getProcessingMethod().setDeisotoping(asBool_(value));
							  break;
							case DECONVOLUTED:
								exp_->getProcessingMethod().setChargeDeconvolution(asBool_(value));
								break;
							case CENTROIDED:
								if (asBool_(value)) spec.setType(SpectrumSettings::PEAKS);
								break;
							case COLLENERGY:
								spec.getPrecursor().setActivationEnergy(asFloat_(value));
								break;
							}
					}
					
					// check if the scan is in the desired range
					//std::cout << "Range: " << options_.getRTRange() << " RT: " << spec.getRT() << std::endl;
					if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec.getRT()))
					 || options_.hasMSLevels() && !options_.containsMSLevel(spec.getMSLevel()))
					{
						// skip this tag
						spec_ = 0;
						skipTag_();
					} 
					else 
					{
						logger_.setProgress(exp_->size());
						exp_->push_back(spec);
						spec_ = &(exp_->back());
					}
				}
				break;
			case SCANORIGIN:
				tmp_str = getAttributeAsString_(PARENTFILEID, true, qname);
				if (tmp_str != "") 
				{
					// ignore
				}
				
				tmp_str = getAttributeAsString_(NUM, true, qname);
				if (tmp_str != "") 
				{
					// ignore
				}
				break;
			case OPERATOR:
				{
					String first = getAttributeAsString_(FIRST_NAME, true, qname);  
					String last = getAttributeAsString_(LAST_NAME, true, qname);    
					
					exp_->getContacts().insert(exp_->getContacts().end(), ContactPerson());
					exp_->getContacts().back().setFirstName(first);
					exp_->getContacts().back().setLastName(last);
					
					tmp_str = getAttributeAsString_(EMAIL, false, qname);
					if (tmp_str != "") 
					{
						exp_->getContacts().back().setEmail(tmp_str);
					}
					
					String contact_info;
					tmp_str = getAttributeAsString_(PHONE, false, qname);
					if (tmp_str != "") 
					{
						contact_info = "PHONE: " + tmp_str;
					}
					
					tmp_str = getAttributeAsString_(URI, false, qname);
					if (tmp_str != "") 
					{
						contact_info += String(contact_info == "" ? "" : " ") + "URI: " + tmp_str;
					}
					
					// if either one of phone or uri was specified
					if (contact_info != "")
					{
						exp_->getContacts().back().setContactInfo(contact_info);
					}
				}
				break;
			case MANUFACTURER:
				{
					tmp_str = getAttributeAsString_(CATEGORY, true, qname);
					if (tmp_str == enum2str_(TAGMAP,MANUFACTURER))
					{
						tmp_str = getAttributeAsString_(VALUE, true, qname);
						if (tmp_str != "")
						{
							exp_->getInstrument().setVendor(tmp_str);
						}
					}
					else
					{
						error("unknown category in 'manufacturer' tag");
					}
				}
				break;
			case MODEL:
				{
					tmp_str = getAttributeAsString_(CATEGORY, true, qname);
					if (tmp_str == enum2str_(TAGMAP,MODEL))
					{
						tmp_str = getAttributeAsString_(VALUE, true, qname);
						if (tmp_str != "")
						{
							exp_->getInstrument().setModel(tmp_str);
						}
					}
					else
					{
						error("unknown category in 'model' tag");
					}
				}
				break;
			case IONISATION:
				{
					tmp_str = getAttributeAsString_(CATEGORY, true, qname);
					if (tmp_str == enum2str_(TAGMAP,IONISATION))
					{
						tmp_str = getAttributeAsString_(VALUE, true, qname);
						if (tmp_str != "")
						{
							exp_->getInstrument().getIonSource().setIonizationMethod((IonSource::IonizationMethod) str2enum_(IONTYPEMAP, tmp_str, "ionization type") );
						}
					}
					else
					{
						error("unknown category in 'ionization' tag");
					}
				}
				break;
			case ANALYZER:
				{
					tmp_str = getAttributeAsString_(CATEGORY, true, qname);
					if (tmp_str == enum2str_(TAGMAP,ANALYZER))
					{
						tmp_str = getAttributeAsString_(VALUE, false, qname);
						if (tmp_str != "")
						{
							exp_->getInstrument().getMassAnalyzers().insert(exp_->getInstrument().getMassAnalyzers().end(), MassAnalyzer());
							analyzer_ = &(exp_->getInstrument().getMassAnalyzers().back());
							analyzer_->setType( (MassAnalyzer::AnalyzerType) str2enum_(ANALYZERTYPEMAP, tmp_str, "analyzer type")
							);
						}
					}
					else
					{
						error("unknown category in 'analyzer' tag");
					}
				}
				break;
			case DETECTOR:
				{
					tmp_str = getAttributeAsString_(CATEGORY, true, qname);
					if (tmp_str == enum2str_(TAGMAP,DETECTOR))
					{
						tmp_str = getAttributeAsString_(VALUE, true, qname);
						if (tmp_str != "")
						{
							IonDetector& ion_d = exp_->getInstrument().getIonDetector();
							ion_d.setType( (IonDetector::Type) str2enum_(TYPEMAP, tmp_str, "detector type") );
						}
					}
					else
					{
						error("unknown category in 'detector' tag");
					}
				}
				break;
			case RESOLUTION:
				{
					tmp_str = getAttributeAsString_(CATEGORY, true, qname);
					if (tmp_str == enum2str_(TAGMAP,RESOLUTION))
					{
						tmp_str = getAttributeAsString_(VALUE, true, qname);
						if (tmp_str != "")
						{
							if (analyzer_ == 0) break;
							analyzer_->setResolutionMethod(
								(MassAnalyzer::ResolutionMethod)
								str2enum_(RESMETHODMAP, tmp_str, "resolution method"));
						}
					}
					else
					{
						error("unknown category in 'resolution' tag");
					}
				}
				break;
			case DATAPROCESSING:
					//optinal attributes
					for (UInt i=0; i<attributes.getLength(); i++)
					{
						int att = str2enum_(ATTMAP,xercesc::XMLString::transcode(attributes.getQName(i)),"dataprocessing attribute");
						String value = xercesc::XMLString::transcode(attributes.getValue(i));
						switch (att)
							{
								case DEISOTOPED:
								  exp_->getProcessingMethod().setDeisotoping(asBool_(value));
									break;
								case DECONVOLUTED:
									exp_->getProcessingMethod().setChargeDeconvolution(asBool_(value));
									break;
								case CENTROIDED:
									exp_->getProcessingMethod().setSpectrumType((SpectrumSettings::SpectrumType)
										str2enum_(PEAKPROCMAP,value,"peak processing"));
									break;
								case SPOT_INTEGRATION:
									// TODO
									break;
								case INTENSITY_CUTOFF:
									exp_->getProcessingMethod().setIntensityCutoff(asDouble_(value));
									break;
							}
					}
					break;
			case NAMEVALUE:
				{
					String name = getAttributeAsString_(NAME, true, qname);
					String value = getAttributeAsString_(VALUE, true, qname);
					
					if (name == "")
					{
						break;
					}
					else if (value == "")
					{
						break;
					}
					
					if (is_parser_in_tag_[INSTRUMENT])
					{
						setAddInfo_(exp_->getInstrument(), name, value, "Instrument.Comment");
					}
					else if (is_parser_in_tag_[SCAN] )
					{
						setAddInfo_(	*spec_, name, value, "Instrument.Comment");
					}
					else
					{
						std::cout << " Warning Unhandled tag: \"" << enum2str_(TAGMAP,NAMEVALUE) << "\"" << std::endl;
					}
				}
				break;
			case PROCESSING:
				{
					String name = getAttributeAsString_(NAME, true, qname);
					String type = getAttributeAsString_(TYPE, true, qname);
					String value = getAttributeAsString_(VALUE, true, qname);
					
					if (value != "" && name != "")
					{
						setAddInfo_(exp_->getProcessingMethod(), name + "#" + type, value, "Processing.Comment");
					}
				}
				break;
		}
		
		//std::cout << " -- !Start -- " << std::endl;
	}


	template <typename MapType>
	void MzXMLHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
  	//std::cout << " -- End -- " << xercesc::XMLString::transcode(qname) << " -- " << std::endl;
  	
  	bool skip = skip_tag_.top();
		int tag = leaveTag(qname);
		if (skip) return;
		
		if (tag==MZXML)
		{
			logger_.endProgress();
		}
		
		if (tag==INSTRUMENT && analyzer_)
		{
			analyzer_ = 0;
		}
		
		if (tag==PEAKS && spec_ != 0)
		{
			//std::cout << "reading scan" << std::endl;
			if (char_rest_=="") // no peaks
			{
				return;
			}
			if (precision_==DOUBLE)		//precision 64
			{
				std::vector<DoubleReal> data;
				decoder_.decode(char_rest_, Base64::BIGENDIAN, data);
				char_rest_ = "";
				PeakType peak;
				//push_back the peaks into the container
				for (UInt n = 0 ; n < ( 2 * peak_count_) ; n += 2)
				{
					// check if peak in in the specified range
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(data[n])))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(data[n+1]))))
					{
						peak.setPosition(data[n]);
						peak.setIntensity(data[n+1]);
						spec_->push_back(peak);
					}
 				}
			}
			else	//precision 32
			{
				std::vector<Real> data;
				decoder_.decode(char_rest_, Base64::BIGENDIAN, data);
				char_rest_ = "";
				PeakType peak;
				//push_back the peaks into the container
				for (UInt n = 0 ; n < (2 * peak_count_) ; n += 2)
				{
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(data[n])))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(data[n+1]))))
					{
						peak.setPosition(data[n]);
						peak.setIntensity(data[n+1]);
						spec_->push_back(peak);
					}
				}
			}
		}
		//std::cout << " -- End -- " << std::endl;
  }

	template <typename MapType>
	void MzXMLHandler<MapType>::writeTo(std::ostream& os)
	{
		//determine how many spectra there are (count only those with peaks)
		UInt count_tmp_  = 0;
		for (UInt s=0; s<cexp_->size(); s++)
		{
			const SpectrumType& spec = (*cexp_)[s];
			if (spec.size()!=0) ++count_tmp_;
		}
		logger_.startProgress(0,cexp_->size(),"storing mzXML file");
		os  << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
			 << "<mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_2.0\" "
			 << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
			 << "xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_2.0 "
			 << "http://sashimi.sourceforge.net/schema_revision/mzXML_2.0/mzXML_idx_2.0.xsd\">\n"
			 << "\t<msRun scanCount=\"" << count_tmp_ << "\">\n"
			 << "\t\t<parentFile fileName=\"" << cexp_->getSourceFile().getNameOfFile()
			 //file type is an enum in mzXML => search for 'raw' string
			 << "\" fileType=\"";
			 String tmp_string = cexp_->getSourceFile().getFileType();
			 tmp_string.toLower();
			 if (tmp_string.hasSubstring("raw"))
			 {
			 	os << "RAWData";
			 }
			 else
			 {
			 	os << "processedData";
			 }
			 //Sha1 checksum must have 40 characters => create a fake if it is unknown
			 os << "\" fileSha1=\"";
			 tmp_string = cexp_->getSourceFile().getSha1();
			 if (cexp_->getSourceFile().getSha1().size()!=40)
			 {
			 	 os << "0000000000000000000000000000000000000000";
			 }
			 else
			 {
			   os << cexp_->getSourceFile().getSha1();
			 }
			 os  << "\"/>\n";

		if (cexp_->getInstrument() != Instrument())
		{
			const Instrument& inst = cexp_->getInstrument();
			os << "\t\t<msInstrument>\n"
				 << "\t\t\t<msManufacturer category=\"msManufacturer\" value=\""
				 <<	inst.getVendor() << "\"/>\n"
				 << "\t\t\t<msModel category=\"msModel\" value=\""
				 << inst.getModel() << "\"/>\n"
				 << "\t\t\t<msIonisation category=\"msIonisation\" value=\""
				 << enum2str_(IONTYPEMAP,inst.getIonSource().getIonizationMethod())
				 << "\"/>\n";

			const std::vector<MassAnalyzer>& analyzers = inst.getMassAnalyzers();
			if ( analyzers.size()>0 )
			{
				os << "\t\t\t<msMassAnalyzer category=\"msMassAnalyzer\" value=\""
					 << enum2str_(ANALYZERTYPEMAP,analyzers[0].getType())  << "\"/>\n";
			}
			else
			{
				std::cout << " Warning: mzXML supports only one analyzer! Skipping the other " << (analyzers.size()-1) << "mass analyzers." << std::endl;
			}
			os << "\t\t\t<msDetector category=\"msDetector\" value=\""
				 << enum2str_(TYPEMAP,inst.getIonDetector().getType()) << "\"/>\n";
			try
			{
				String type = inst.getMetaValue("#InstSoftwareType").toString();
				//invalid type is resetted to 'processing' as it fits all actions
				if (type!="acquisition" && type!="conversion" && type!="processing")
				{
					type = "processing";
				}
				String name = inst.getMetaValue("#InstSoftware").toString();
				String version = inst.getMetaValue("#InstSoftwareVersion").toString();
				String str = inst.getMetaValue("#InstSoftwareTime").toString();
				String time(str);
				time.substitute(' ', 'T');
				os << "\t\t\t<software type=\"" << type
					 << "\" name=\"" << name
					 << "\" version=\"" << version << "\"";
				if (time != "")
				{
					os << " completionTime=\"" << time << "\"";
				}
				os << "/>\n";
			}
			catch(Exception::InvalidValue exception)
			{

			}
			
			if ( analyzers.size()>0 )
			{
				if (analyzers[0].getResolutionMethod())
					os << "\t\t\t<msResolution category=\"msResolution\" value=\""
				 		 << enum2str_(RESMETHODMAP,analyzers[0].getResolutionMethod()) << "\"/>\n";
			}
			else
			{
				std::cout << "Warning: mzXML supports only one analyzer! Skipping the other " << (analyzers.size()-1) << "mass analyzers." << std::endl;
			}
			
			if ( cexp_->getContacts().size()>0 )
			{
				const ContactPerson& cont = cexp_->getContacts()[0];
				
				os << "\t\t\t<operator first=\"" << cont.getFirstName() << "\" last=\"" << cont.getLastName();
				
				String info = cont.getContactInfo();
				std::string::size_type phone = info.find("PHONE:");
				std::string::size_type uri = info.find("URI:");
				if (phone != std::string::npos)
				{
					UInt end = uri != std::string::npos ? uri : info.size();
					os << "\" phone=\"" << info.substr(phone + 6, end - phone + 6);
				}
				
				if (cont.getEmail() != "")
				{
					os << "\" email=\"" << cont.getEmail();
				}
				
				if (uri != std::string::npos)
				{
					UInt uri = info.find("URI:");
					os << "\" URI=\"" << info.substr(uri+4).trim();
				}
				
				os << "\"/>\n";
			}
			writeUserParam_(os,inst,3);
			try
			{
				DataValue com = inst.getMetaValue("#Comment");
				if (!com.isEmpty()) os << "\t\t\t<comment>" << com << "</comment>\n";
			}
			catch(Exception::InvalidValue exception)
			{

			}
			os << "\t\t</msInstrument>\n";
		}

		const Software& software = cexp_->getSoftware();
		os << "\t\t<dataProcessing deisotoped=\""
			 << cexp_->getProcessingMethod().getDeisotoping()
			 << "\" chargeDeconvoluted=\""
			 << cexp_->getProcessingMethod().getChargeDeconvolution()
			 << "\" centroided=\""
			 << enum2str_(PEAKPROCMAP,cexp_->getProcessingMethod().getSpectrumType())
			 << "\" intensityCutoff=\""
			 << cexp_->getProcessingMethod().getIntensityCutoff()
			 << "\">\n"
			 << "\t\t\t<software type=\"processing\" name=\"" << software.getName()
			 << "\" version=\"" << software.getVersion();

		if (software.getCompletionTime() != DateTime())
		{
			String tmp;
			software.getCompletionTime().get(tmp);
			String qtmp(tmp);
			qtmp.substitute(' ', 'T');
			os << "\" completionTime=\"" << qtmp;
		}
		os << "\"/>\n";
		writeUserParam_(os,cexp_->getProcessingMethod(),3,"processingOperation");

		try
		{
			DataValue com = cexp_->getProcessingMethod().getMetaValue("#Comment");
			if (!com.isEmpty()) os << "\t\t\t<comment>" << com << "</comment>\n";
		}
		catch(Exception::InvalidValue exception)
		{

		}

		os << "\t\t</dataProcessing>\n";
		
		std::stack<UInt> open_scans;
		
		// write scans
		for (UInt s=0; s<cexp_->size(); s++)
		{
			logger_.setProgress(s);
			const SpectrumType& spec = (*cexp_)[s];
						
			int ms_level = spec.getMSLevel();
			open_scans.push(ms_level);
			
			os << String(ms_level+1,'\t')
				 << "<scan num=\"" << spec_write_counter_++ << "\" msLevel=\""
				 << ms_level << "\" peaksCount=\""
				 << spec.size() << "\" polarity=\""
				 << enum2str_(POLARITYMAP,spec.getInstrumentSettings().getPolarity());
			
			if (spec.getInstrumentSettings().getScanMode())
			{
				os << "\" scanType=\""
					 << enum2str_(SCANMODEMAP,spec.getInstrumentSettings().getScanMode());
			}
			os << "\" retentionTime=\"PT"
				 << spec.getRT() << "S\"";
			if (spec.getInstrumentSettings().getMzRangeStart()!=0)
				os << " startMz=\"" << spec.getInstrumentSettings().getMzRangeStart() << "\"";
			if (spec.getInstrumentSettings().getMzRangeStop()!=0)
				os << " endMz=\"" << spec.getInstrumentSettings().getMzRangeStop() << "\"";
			os << ">\n";

			const PrecursorPeakType& peak = spec.getPrecursorPeak();
			if (peak!= PrecursorPeakType())
			{
				os << String(ms_level+2,'\t') << "<precursorMz precursorIntensity=\""
					 << peak.getIntensity();
				if (peak.getCharge()!=0)
					os << "\" precursorCharge=\"" << peak.getCharge();
				os << "\">"
				 	 << peak.getPosition()[0] << "</precursorMz>\n";
			}

			if (spec.size() > 0)
			{
				os << String(ms_level+2,'\t') << "<peaks precision=\"32\""
					 << " byteOrder=\"network\" pairOrder=\"m/z-int\">";
				
				//std::cout << "Writing scan " << s << std::endl;
				std::vector<Real> tmp;
				for (UInt i=0; i<spec.size(); i++)
				{
					tmp.push_back(spec.getContainer()[i].getMZ());
					tmp.push_back(spec.getContainer()[i].getIntensity());
				}
				
				std::string encoded;
				decoder_.encode(tmp, Base64::BIGENDIAN, encoded);
				os << encoded << "</peaks>\n";
			}
			else
			{
				os << String(ms_level+2,'\t') << "<peaks precision=\"32\""
					 << " byteOrder=\"network\" pairOrder=\"m/z-int\" xsi:nil=\"1\"/>\n";
			}
			
			writeUserParam_(os,spec,ms_level+2);
			if (spec.getComment() != "")
			{
				os << String(ms_level+2,'\t') << "<comment>" << spec.getComment() << "</comment>\n";
			}
			
			//check MS level of next scan and close scans (scans can be nested)
			int next_ms_level = 0;
			if (s < cexp_->size()-1)
			{
				next_ms_level = ((*cexp_)[s+1]).getMSLevel();
			}
			//std::cout << "scan: " << s << " this: " << ms_level << " next: " << next_ms_level << std::endl;
			if (next_ms_level <= ms_level)
			{
				for (Int i = 0; i<= ms_level-next_ms_level && !open_scans.empty(); ++i)
				{
					os << String(ms_level-i+1,'\t') << "</scan>\n";
					open_scans.pop();
				}
			}
		}

		os << "\t</msRun>\n"
			 << "\t<indexOffset>0</indexOffset>\n"
			 << "</mzXML>\n";
		
		logger_.endProgress();
	}

	} // namespace Internal

} // namespace OpenMS

#endif
