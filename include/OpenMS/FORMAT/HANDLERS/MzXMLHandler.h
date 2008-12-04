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

#ifndef OPENMS_FORMAT_HANDLERS_MZXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZXMLHANDLER_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
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
			: public XMLHandler
	  {
	    public:
	      /**@name Constructors and destructor */
	      //@{
	      /// Constructor for a write-only handler
	      MzXMLHandler(MapType& exp, const String& filename, const String& version, ProgressLogger& logger)
				: XMLHandler(filename,version),
					exp_(&exp),	
					cexp_(0),
					decoder_(),
					peak_count_(0),
					char_rest_(),
					skip_spectrum_(false),
					spec_write_counter_(1),
					logger_(logger)
	  		{
	  			cv_terms_.resize(6);
	  			//Polarity
					String("any;+;-").split(';',cv_terms_[0]);
					//Scan type
					String(";zoom;Full;SIM;SRM;CRM;CNG;CNL;PRODUCT;PRECURSOR;ER").split(';',cv_terms_[1]);
					//Ionization method
					String(";ESI;EI;CI;FAB;TSP;MALDI;FD;FI;PD;SI;TI;API;ISI;CID;CAD;HN;APCI;APPI;ICP").split(';',cv_terms_[2]);
					//Mass analyzer
					String(";Quadrupole;Quadrupole Ion Trap;;;TOF;Magnetic Sector;FT-ICR;").split(';',cv_terms_[3]);
					//Detector
					String(";EMT;Daly;;Faraday Cup;;;;Channeltron").split(';',cv_terms_[4]);
					//Resolution method
					String(";FWHM;TenPercentValley;Baseline").split(';',cv_terms_[5]);
				}
	
	      /// Constructor for a read-only handler
	      MzXMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger)
				: XMLHandler(filename,version),
					exp_(0), 
					cexp_(&exp),
					decoder_(),
					peak_count_(0),
					char_rest_(),
					skip_spectrum_(false),
					spec_write_counter_(1),
					logger_(logger)
	  		{
	  			cv_terms_.resize(6);
	  			//Polarity
					String("any;+;-").split(';',cv_terms_[0]);
					//Scan type
					String(";zoom;Full;SIM;SRM;CRM;CNG;CNL;PRODUCT;PRECURSOR;ER").split(';',cv_terms_[1]);
					//Ionization method
					String(";ESI;EI;CI;FAB;TSP;MALDI;FD;FI;PD;SI;TI;API;ISI;CID;CAD;HN;APCI;APPI;ICP").split(';',cv_terms_[2]);
					//Mass analyzer
					String(";Quadrupole;Quadrupole Ion Trap;;;TOF;Magnetic Sector;FT-ICR;").split(';',cv_terms_[3]);
					//Detector
					String(";EMT;Daly;;Faraday Cup;;;;Channeltron").split(';',cv_terms_[4]);
					//Resolution method
					String(";FWHM;TenPercentValley;Baseline").split(';',cv_terms_[5]);
				}

	  		/// Destructor
	      virtual ~MzXMLHandler()
	      {
	      }
	      //@}
				
				// Docu in base class
	      virtual void endElement( const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname);
				
				// Docu in base class
	      virtual void startElement(const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname, const xercesc::Attributes& attributes);
				
				// Docu in base class
	      virtual void characters(const XMLCh* const chars, const unsigned int length);
	
	  		///Write the contents to a stream
				void writeTo(std::ostream& os);
				
				///Sets the options
				void setOptions(const PeakFileOptions& options)
				{
					options_ = options; 
				}
	
	    protected:
				
        /// Peak type
        typedef typename MapType::PeakType PeakType;
        /// Spectrum type
        typedef MSSpectrum<PeakType, std::allocator<PeakType> > SpectrumType;        
        
				typedef typename SpectrumType::Iterator  PeakIterator;
				typedef typename SpectrumType::PrecursorPeakType PrecursorPeakType;
				
				/// map pointer for reading
				MapType* exp_;
				/// map pointer for writing
				const MapType* cexp_;
				
				/// Options for loading and storing
				PeakFileOptions options_;
				
				/**@name temporary datastructures to hold parsed data */
		    //@{
				Base64 decoder_;
				UInt peak_count_;
				String precision_;
				String char_rest_;
				//@}
				
				/// Flag that indicates whether this spectrum should be skipped (due to options)
				bool skip_spectrum_;
				
				/// spectrum counter (spectra without peaks are not written)
				UInt spec_write_counter_;
				
				/// Progress logging class
				const ProgressLogger& logger_;
		
				/// write metaInfo to xml (usually in nameValue-tag)
				inline void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, int indent=4, String tag="nameValue")
				{
					std::vector<String> keys;  // Vector to hold keys to meta info
					meta.getKeys(keys);
		
					for (std::vector<String>::const_iterator it = keys.begin(); it!=keys.end(); ++it)
					{
						if ( (*it)[0] != '#')  // internally used meta info start with '#'
						{
							os << String(indent,'\t') << "<" << tag << " name=\"" << *it << "\" value=\"" << meta.getMetaValue(*it) << "\"/>\n";
						}
					}
				}
			
			private:
				/// Not implemented
				MzXMLHandler();
	  };
	
		//--------------------------------------------------------------------------------
	
		template <typename MapType>
	  void MzXMLHandler<MapType>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	  {
	  	static const XMLCh* s_value = xercesc::XMLString::transcode("value");
	  	static const XMLCh* s_count = xercesc::XMLString::transcode("scanCount");
	  	static const XMLCh* s_type = xercesc::XMLString::transcode("type");
	  	static const XMLCh* s_name = xercesc::XMLString::transcode("name");
	  	static const XMLCh* s_version = xercesc::XMLString::transcode("version");
	  	static const XMLCh* s_filename = xercesc::XMLString::transcode("fileName");
	  	static const XMLCh* s_filetype = xercesc::XMLString::transcode("fileType");
	  	static const XMLCh* s_filesha1 = xercesc::XMLString::transcode("fileSha1");
	  	static const XMLCh* s_completiontime = xercesc::XMLString::transcode("completionTime");
	  	static const XMLCh* s_precision = xercesc::XMLString::transcode("precision");
	  	static const XMLCh* s_byteorder = xercesc::XMLString::transcode("byteOrder");
	  	static const XMLCh* s_pairorder = xercesc::XMLString::transcode("pairOrder");
	  	static const XMLCh* s_precursorintensity = xercesc::XMLString::transcode("precursorIntensity");
	  	static const XMLCh* s_precursorcharge = xercesc::XMLString::transcode("precursorCharge");
	  	static const XMLCh* s_windowwideness = xercesc::XMLString::transcode("windowWideness");
	  	static const XMLCh* s_mslevel = xercesc::XMLString::transcode("msLevel");
	  	static const XMLCh* s_peakscount = xercesc::XMLString::transcode("peaksCount");
	  	static const XMLCh* s_polarity = xercesc::XMLString::transcode("polarity");
	  	static const XMLCh* s_scantype = xercesc::XMLString::transcode("scanType");
	  	static const XMLCh* s_retentiontime = xercesc::XMLString::transcode("retentionTime");
	  	static const XMLCh* s_collisionenergy = xercesc::XMLString::transcode("collisionEnergy");
	  	static const XMLCh* s_startmz = xercesc::XMLString::transcode("startMz");
	  	static const XMLCh* s_endmz = xercesc::XMLString::transcode("endMz");
	  	static const XMLCh* s_first = xercesc::XMLString::transcode("first");
	  	static const XMLCh* s_last = xercesc::XMLString::transcode("last");
	  	static const XMLCh* s_phone = xercesc::XMLString::transcode("phone");
	  	static const XMLCh* s_email = xercesc::XMLString::transcode("email");
	  	static const XMLCh* s_uri = xercesc::XMLString::transcode("URI");
	  	static const XMLCh* s_num = xercesc::XMLString::transcode("num");
	  	static const XMLCh* s_intensitycutoff = xercesc::XMLString::transcode("intensityCutoff");
	  	static const XMLCh* s_centroided = xercesc::XMLString::transcode("centroided");
	  	static const XMLCh* s_deisotoped = xercesc::XMLString::transcode("deisotoped");
	  	static const XMLCh* s_chargedeconvoluted = xercesc::XMLString::transcode("chargeDeconvoluted");
	  	
	  	static UInt scan_count = 0;
	  	
	  	String tag = sm_.convert(qname);
	  	open_tags_.push_back(tag);
	  	//std::cout << " -- Start -- "<< tag << " -- " << std::endl;
	  	
	  	//Skip all tags until the the next scan
	  	if (skip_spectrum_ && tag!="scan") return;
	  	
			if (tag=="msRun")
			{
				Int count = 0;
				optionalAttributeAsInt_(count, attributes, s_count);
				exp_->reserve(count);
				logger_.startProgress(0,count,"loading mzXML file");
				scan_count = 0;
				//start and end time are xs:duration. This makes no senes => ignore them
			}
			else if (tag=="parentFile")
			{
				SourceFile sf;
				sf.setNameOfFile(attributeAsString_(attributes, s_filename));
				sf.setFileType(attributeAsString_(attributes, s_filetype));
				sf.setChecksum(attributeAsString_(attributes, s_filesha1), SourceFile::SHA1);		
				exp_->getSourceFiles().push_back(sf);			
			}
			else if (tag=="software")
			{
				String& parent_tag = *(open_tags_.end()-2);
				if (parent_tag=="dataProcessing")
				{
					exp_->getDataProcessing().back().getSoftware().setVersion(attributeAsString_(attributes, s_version));
					exp_->getDataProcessing().back().getSoftware().setName(attributeAsString_(attributes, s_name));
					exp_->getDataProcessing().back().setMetaValue("#type",String(attributeAsString_(attributes, s_type)));
					
					String time;
					optionalAttributeAsString_(time,attributes,s_completiontime);
					exp_->getDataProcessing().back().setCompletionTime( asDateTime_(time) );
				}
				else if (parent_tag=="msInstrument")
				{
					exp_->getInstrument().getSoftware().setVersion(attributeAsString_(attributes, s_version));
					exp_->getInstrument().getSoftware().setName(attributeAsString_(attributes, s_name));
				}
			}
			else if (tag=="peaks")
			{
				precision_ = attributeAsString_(attributes, s_precision);
				if (precision_!="32" && precision_!="64")
				{
					error(LOAD, String("Invalid precision '") + precision_ + "' in element 'peaks'");
				}
				String byte_order;
				optionalAttributeAsString_(byte_order, attributes, s_byteorder);
				if (byte_order!="network")
				{
					error(LOAD, String("Invalid or missing byte order '") + byte_order + "' in element 'peaks'. Must be 'network'!");
				}
				String pair_order;
				optionalAttributeAsString_(pair_order, attributes, s_pairorder);
				if (pair_order!="m/z-int")
				{
					error(LOAD, String("Invalid or missing pair order '") + pair_order + "' in element 'peaks'. Must be 'm/z-int'!");
				}
			}
			else if (tag=="precursorMz")
			{
				try
				{
					exp_->back().getPrecursorPeak().setIntensity( attributeAsDouble_(attributes, s_precursorintensity) );
				}
				catch (Exception::ParseError& e)
				{
					error(LOAD, "Mandatory attribute precursorMz not found! Setting precursor intensity to 0 - trying to continue");
					exp_->back().getPrecursorPeak().setIntensity(0.0);
				}
				
				Int charge = 0;
				optionalAttributeAsInt_(charge, attributes, s_precursorcharge);
				exp_->back().getPrecursorPeak().setCharge(charge);
				
				DoubleReal window = 0;
				optionalAttributeAsDouble_(window, attributes, s_windowwideness);
				exp_->back().getPrecursor().setWindowSize(window);
			}
			else if (tag=="scan")
			{
				skip_spectrum_ = false;
				
				if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				// check if the scan is in the desired MS / RT range
				UInt ms_level = attributeAsInt_(attributes, s_mslevel);
				//parse retention time and convert it from xs:duration to seconds
				DoubleReal retention_time = 0.0;
				String time_string = "";
				if (optionalAttributeAsString_(time_string, attributes, s_retentiontime))
				{
					time_string = time_string.suffix('T');
					//std::cout << "Initial trim: " << time_string << std::endl;
					if (time_string.has('H'))
					{
						retention_time += 3600*asDouble_(time_string.prefix('H'));
						time_string = time_string.suffix('H');
						//std::cout << "After H: " << time_string << std::endl;
					}
					if (time_string.has('M'))
					{
						retention_time += 60*asDouble_(time_string.prefix('M'));
						time_string = time_string.suffix('M');
						//std::cout << "After M: " << time_string << std::endl;
					}
					if (time_string.has('S'))
					{
						retention_time += asDouble_(time_string.prefix('S'));
						time_string = time_string.suffix('S');
						//std::cout << "After S: " << time_string << std::endl;
					}
				}

				logger_.setProgress(scan_count);
				
				if ( (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(retention_time)))
				 || (options_.hasMSLevels() && !options_.containsMSLevel(ms_level)) )
				{
					// skip this tag
					skip_spectrum_ = true;
					++scan_count;			
					return;
				}

				//Add a new spectrum and set MS level and RT
				exp_->resize(exp_->size()+1);
				exp_->back().setMSLevel(ms_level);
				exp_->back().setRT(retention_time);
				exp_->back().setNativeID(String("scan=") + attributeAsString_(attributes, s_num));
				//peak count == twice the scan size
				peak_count_ = attributeAsInt_(attributes, s_peakscount);
				exp_->back().reserve(peak_count_);
				
				//centroided, chargeDeconvoluted, deisotoped are ignored

				//other optional attributes
				InstrumentSettings::ScanWindow window;
				optionalAttributeAsDouble_(window.begin, attributes, s_startmz);
				optionalAttributeAsDouble_(window.end, attributes, s_endmz);
				if (window.begin!=0.0 || window.end!=0.0 )
				{
					exp_->back().getInstrumentSettings().getScanWindows().push_back(window);
				}
				
				DoubleReal tmp = 0.0;
				optionalAttributeAsDouble_(tmp, attributes, s_collisionenergy);
				exp_->back().getPrecursor().setActivationEnergy(tmp);
				
				String polarity = "any";
				optionalAttributeAsString_(polarity, attributes, s_polarity);
				exp_->back().getInstrumentSettings().setPolarity( (IonSource::Polarity) cvStringToEnum_(0,polarity,"polarity") );
				
				String type = "";
				optionalAttributeAsString_(type, attributes, s_scantype);
				if (type=="EMS") type=""; //Enhanced MS (ABI - Sashimi converter)
				if (type=="EPI") type=""; //Enhanced Product Ion (ABI - Sashimi converter)
				exp_->back().getInstrumentSettings().setScanMode( (InstrumentSettings::ScanMode) cvStringToEnum_(1,type,"scanType") );
				
				++scan_count;
			}
			else if (tag=="operator")
			{
				exp_->getContacts().resize(1);
				exp_->getContacts().back().setFirstName(attributeAsString_(attributes, s_first));
				exp_->getContacts().back().setLastName(attributeAsString_(attributes, s_last));
				
				String tmp = "";
				optionalAttributeAsString_(tmp, attributes,s_email);
				exp_->getContacts().back().setEmail(tmp);

				tmp = "";
				optionalAttributeAsString_(tmp, attributes,s_phone);
				if (tmp != "") 
				{
					exp_->getContacts().back().setMetaValue("#phone", tmp);
				}
				
				tmp = "";
				optionalAttributeAsString_(tmp, attributes,s_uri);
				exp_->getContacts().back().setURL(tmp);
			}
			else if (tag=="msManufacturer")
			{
				exp_->getInstrument().setVendor(attributeAsString_(attributes, s_value));
			}
			else if (tag=="msModel")
			{
				exp_->getInstrument().setModel(attributeAsString_(attributes, s_value));
			}
			else if (tag=="msIonisation")
			{
				exp_->getInstrument().getIonSources().resize(1);
				exp_->getInstrument().getIonSources()[0].setIonizationMethod((IonSource::IonizationMethod) cvStringToEnum_(2, attributeAsString_(attributes, s_value), "msIonization") );
			}
			else if (tag=="msMassAnalyzer")
			{
				exp_->getInstrument().getMassAnalyzers().resize(1);
				exp_->getInstrument().getMassAnalyzers()[0].setType( (MassAnalyzer::AnalyzerType) cvStringToEnum_(3, attributeAsString_(attributes, s_value), "msMassAnalyzer") );
			}
			else if (tag=="msDetector")
			{
				exp_->getInstrument().getIonDetectors().resize(1);
				exp_->getInstrument().getIonDetectors()[0].setType( (IonDetector::Type) cvStringToEnum_(4, attributeAsString_(attributes, s_value), "msDetector") );
			}
			else if (tag=="msResolution")
			{
				exp_->getInstrument().getMassAnalyzers()[0].setResolutionMethod( (MassAnalyzer::ResolutionMethod) cvStringToEnum_(5, attributeAsString_(attributes, s_value), "msResolution"));
			}
			else if (tag=="dataProcessing")
			{
				exp_->getDataProcessing().push_back(DataProcessing());
				String boolean = "";
				optionalAttributeAsString_(boolean, attributes, s_deisotoped);
				if (boolean == "true" || boolean == "1")
				{
					exp_->getDataProcessing().back().getProcessingActions().insert(DataProcessing::DEISOTOPING);
				}
				
				boolean = "";
				optionalAttributeAsString_(boolean, attributes, s_chargedeconvoluted);
				if (boolean == "true" || boolean == "1")
				{
					exp_->getDataProcessing().back().getProcessingActions().insert(DataProcessing::CHARGE_DECONVOLUTION);
				}
				
				DoubleReal cutoff = 0.0;
				optionalAttributeAsDouble_(cutoff, attributes, s_intensitycutoff);
				if (cutoff!=0.0)
				{
					exp_->getDataProcessing().back().setMetaValue("#intensity_cutoff",cutoff);
				}
					
				boolean = "";
				optionalAttributeAsString_(boolean, attributes, s_centroided);
				if (boolean == "true" || boolean == "1")
				{
					exp_->getDataProcessing().back().getProcessingActions().insert(DataProcessing::PEAK_PICKING);
				}
			}
			else if (tag=="nameValue")
			{
				String name = "";
				optionalAttributeAsString_(name, attributes, s_name);
				if (name == "") return;

				String value = "";
				optionalAttributeAsString_(value, attributes, s_value);
				
				String& parent_tag = *(open_tags_.end()-2);
								
				if (parent_tag == "msInstrument")
				{
					exp_->getInstrument().setMetaValue(name, value);
				}
				else if (parent_tag == "scan")
				{
					exp_->back().setMetaValue(name, value);
				}
				else
				{
					std::cout << " Warning: Unexpected tag 'nameValue' in tag '" << parent_tag << "'" << std::endl;
				}
			}
			else if (tag=="processingOperation")
			{
				String name = "";
				optionalAttributeAsString_(name, attributes, s_name);
				if (name == "") return;

				String value = "";
				optionalAttributeAsString_(value, attributes, s_value);
				
				exp_->getDataProcessing().back().setMetaValue(name, value);
			}
			
			//std::cout << " -- !Start -- " << std::endl;
		}
	
	
		template <typename MapType>
		void MzXMLHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	  {
	  	//std::cout << " -- End -- " << sm_.convert(qname) << " -- " << std::endl;
	  	
	  	static const XMLCh* s_mzxml = xercesc::XMLString::transcode("mzXML");
	  	static const XMLCh* s_peaks = xercesc::XMLString::transcode("peaks");
	  	
			open_tags_.pop_back();
			
			//abort if this scan should be skipped
			if (skip_spectrum_) return;
			
			if (equal_(qname,s_mzxml))
			{
				logger_.endProgress();
			}
			else if (equal_(qname,s_peaks))
			{
				//std::cout << "reading scan" << std::endl;
				if (char_rest_=="") // no peaks
				{
					return;
				}
				
				//remove whitespaces from binary data
				//this should not be necessary, but linebreaks inside the base64 data are unfortunately no exception
				char_rest_.removeWhitespaces();
				
				if (precision_=="64")
				{
					std::vector<DoubleReal> data;
					decoder_.decode(char_rest_, Base64::BYTEORDER_BIGENDIAN, data);
					char_rest_ = "";
					PeakType peak;
					//push_back the peaks into the container
					for (UInt n = 0 ; n < ( 2 * peak_count_) ; n += 2)
					{
						// check if peak in in the specified m/z  and intensity range
						if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(data[n])))
						 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(data[n+1]))))
						{
							peak.setPosition(data[n]);
							peak.setIntensity(data[n+1]);
							exp_->back().push_back(peak);
						}
	 				}
				}
				else	//precision 32
				{
					std::vector<Real> data;
					decoder_.decode(char_rest_, Base64::BYTEORDER_BIGENDIAN, data);
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
							exp_->back().push_back(peak);
						}
					}
				}
			}
			//std::cout << " -- End -- " << std::endl;
			sm_.clear();
	  }
	
		template <typename MapType>
	  void MzXMLHandler<MapType>::characters(const XMLCh* const chars, unsigned int /*length*/)
	  {
	  	//Abort if this spectrum should be skipped
	    if (skip_spectrum_) return;
			
			char* transcoded_chars = sm_.convert(chars);
			
			if(open_tags_.back()=="peaks")
			{
				//chars may be split to several chunks => concatenate them
				char_rest_ += transcoded_chars;
			}
			else if (	open_tags_.back()=="offset" || open_tags_.back()=="indexOffset" || open_tags_.back()=="sha1")
			{
				
			}
			else if (	open_tags_.back()=="precursorMz")
			{
				exp_->back().getPrecursorPeak().getPosition()[0] = asDouble_(transcoded_chars);
			}
			else if (	open_tags_.back()=="comment")
			{
				String parent_tag = *(open_tags_.end()-2);
				//std::cout << "- Comment of parent " << parent_tag << std::endl;
					
				if (parent_tag=="msInstrument")
				{
					exp_->getInstrument().setMetaValue("#Comment" , String(transcoded_chars));
				}
				else if (parent_tag=="dataProcessing")
				{
					//this is currently ignored
				}
				else if (parent_tag=="scan")
				{
					exp_->back().setComment( transcoded_chars );
				}
				else if (String(transcoded_chars).trim()!="")
				{
					warning(LOAD, String("Unhandled comment '") + transcoded_chars + "' in element '" + open_tags_.back() + "'");
				}
			}
			else if (String(transcoded_chars).trim()!="")
			{
				warning(LOAD, String("Unhandled character content '") + transcoded_chars + "' in element '" + open_tags_.back() + "'");
			}
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
			if (count_tmp_==0) ++count_tmp_;
			logger_.startProgress(0,cexp_->size(),"storing mzXML file");
			os  << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
				 << "<mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_2.1\" "
				 << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
				 << "xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_2.1 "
				 << "http://sashimi.sourceforge.net/schema_revision/mzXML_2.1/mzXML_idx_2.1.xsd\">\n"
				 << "\t<msRun scanCount=\"" << count_tmp_ << "\">\n";
			
			//----------------------------------------------------------------------------------------
			// parent files
			//----------------------------------------------------------------------------------------
			if (cexp_->getSourceFiles().size()==0)
			{
					os << "\t\t<parentFile fileName=\"\" fileType=\"processedData\" fileSha1=\"0000000000000000000000000000000000000000\"/>\n";
			}
			else 
			{
				for (UInt i=0; i< cexp_->getSourceFiles().size(); ++i)
				{
					const SourceFile& sf = cexp_->getSourceFiles()[i];
					os << "\t\t<parentFile fileName=\"" << sf.getNameOfFile() << "\" fileType=\"";
					//file type is an enum in mzXML => search for 'raw' string
					String tmp_string = sf.getFileType();
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
					tmp_string = sf.getChecksum();
					if (sf.getChecksum().size()!=40 || sf.getChecksumType()!=SourceFile::SHA1)
					{
						os << "0000000000000000000000000000000000000000";
					}
					else
					{
						os << sf.getChecksum();
					}
					os  << "\"/>\n";
				}
			}

			//----------------------------------------------------------------------------------------
			//instrument
			//----------------------------------------------------------------------------------------
			if (cexp_->getInstrument() != Instrument() || cexp_->getContacts().size()!=0)
			{
				const Instrument& inst = cexp_->getInstrument();
				os << "\t\t<msInstrument>\n"
					 << "\t\t\t<msManufacturer category=\"msManufacturer\" value=\"" <<	inst.getVendor() << "\"/>\n" << "\t\t\t<msModel category=\"msModel\" value=\"" << inst.getModel() << "\"/>\n";
				if (inst.getIonSources().size()==0 || !inst.getIonSources()[0].getIonizationMethod())
				{
					os << "\t\t\t<msIonisation category=\"msIonisation\" value=\"\"/>\n";
				}
				else
				{
					os << "\t\t\t<msIonisation category=\"msIonisation\" value=\"" << cv_terms_[2][inst.getIonSources()[0].getIonizationMethod()] << "\"/>\n";
				}
				const std::vector<MassAnalyzer>& analyzers = inst.getMassAnalyzers();
				if (analyzers.size()==0 || !analyzers[0].getResolutionMethod())
				{
					os << "\t\t\t<msMassAnalyzer category=\"msMassAnalyzer\" value=\"\"/>\n";
				}
				else
				{
					os << "\t\t\t<msMassAnalyzer category=\"msMassAnalyzer\" value=\"" << cv_terms_[3][analyzers[0].getType()]  << "\"/>\n";
				}
				if (inst.getIonDetectors().size()==0 || !inst.getIonDetectors()[0].getType())
				{
					os << "\t\t\t<msDetector category=\"msDetector\" value=\"\"/>\n";
				}
				else
				{
					os << "\t\t\t<msDetector category=\"msDetector\" value=\"" << cv_terms_[4][inst.getIonDetectors()[0].getType()] << "\"/>\n";					
				}
				os << "\t\t\t<software type=\"acquisition\" name=\"" << inst.getSoftware().getName() << "\" version=\"" << inst.getSoftware().getVersion() << "\"/>\n";
				if (analyzers.size()==0 || !analyzers[0].getResolutionMethod())
				{
					os << "\t\t\t<msResolution category=\"msResolution\" value=\"\"/>\n";
				}
				else
				{
					os << "\t\t\t<msResolution category=\"msResolution\" value=\"" << cv_terms_[5][analyzers[0].getResolutionMethod()] << "\"/>\n";
				}
				
				if ( cexp_->getContacts().size()>0 )
				{
					const ContactPerson& cont = cexp_->getContacts()[0];
					
					os << "\t\t\t<operator first=\"" << cont.getFirstName() << "\" last=\"" << cont.getLastName() <<  "\"" ;
					
					if (cont.getEmail() != "")
					{
						os << " email=\"" << cont.getEmail() << "\"";
					}
					
					if (cont.getURL() != "")
					{
						os << " URI=\"" << cont.getURL() << "\"";
					}
					
					if (cont.metaValueExists("#phone"))
					{
						os << " phone=\"" << (String)(cont.getMetaValue("#phone")) << "\"";
					}
					
					os << "/>\n";
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
			
			//----------------------------------------------------------------------------------------
			//data processing
			//----------------------------------------------------------------------------------------
			if (cexp_->getDataProcessing().size()==0)
			{
				os << "\t\t<dataProcessing>\n"
					 << "\t\t\t<software type=\"processing\" name=\"\" version=\"\"/>\n"
					 << "\t\t</dataProcessing>\n";
			}
			else
			{
				for (UInt i=0; i<cexp_->getDataProcessing().size(); ++i)
				{
					const DataProcessing& data_processing = cexp_->getDataProcessing()[i];
					os << "\t\t<dataProcessing deisotoped=\""
						 << data_processing.getProcessingActions().count(DataProcessing::DEISOTOPING)
						 << "\" chargeDeconvoluted=\""
						 << data_processing.getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION)
						 << "\" centroided=\""
						 << data_processing.getProcessingActions().count(DataProcessing::PEAK_PICKING)
						 << "\"";
					if (data_processing.metaValueExists("#intensity_cutoff"))
					{
						os  << " intensityCutoff=\"" << data_processing.getMetaValue("#intensity_cutoff").toString() << "\"";
					}
					os << ">\n"
						 << "\t\t\t<software type=\"";
					if (data_processing.metaValueExists("#type"))
					{
						os << data_processing.getMetaValue("#type").toString();
					}
					else
					{
						os << "processing";
					}
					
					os << "\" name=\"" << data_processing.getSoftware().getName()
						 << "\" version=\"" << data_processing.getSoftware().getVersion();
			
					if (data_processing.getCompletionTime() != DateTime())
					{
						os << "\" completionTime=\"" << data_processing.getCompletionTime().get().substitute(' ', 'T');
					}
					os << "\"/>\n";
					writeUserParam_(os,data_processing,3,"processingOperation");
			
					os << "\t\t</dataProcessing>\n";	
				}
			}
			
			//check if the nativeID of all spectra are numbers or numbers prefixed with 'scan='
			//If not we need to renumber all spectra.
			bool all_numbers = true;
			bool all_empty = true;
			bool all_prefixed_numbers = true;
			for (UInt s=0; s<cexp_->size(); s++)
			{
				String native_id = (*cexp_)[s].getNativeID();
				if (!native_id.hasPrefix("scan="))
				{
					all_prefixed_numbers = false;
				}
				else
				{
					native_id = native_id.substr(5);
				}
				try
				{
					native_id.toInt();
				}
				catch (Exception::ConversionError&)
				{
					all_numbers = false;
					all_prefixed_numbers = false;
					if (native_id!="")
					{
						all_empty = false;
					}
				}
			}
			//If we need to renumber and the nativeIDs were not empty, warn the user
			if (!all_numbers && !all_empty)
			{
				warning(STORE, "Not all spectrum native IDs are numbers or correctly prefixed with 'scan='. The spectra are renumbered and the native IDs are lost!");
			}
			
			// write scans
			std::stack<UInt> open_scans;
			for (UInt s=0; s<cexp_->size(); s++)
			{
				logger_.setProgress(s);
				const SpectrumType& spec = (*cexp_)[s];
							
				UInt ms_level = spec.getMSLevel();
				open_scans.push(ms_level);

				UInt spectrum_id = s+1;
				if (all_prefixed_numbers)
				{
					spectrum_id = spec.getNativeID().substr(5).toInt();
				}
				else if (all_numbers)
				{
					spectrum_id = spec.getNativeID().toInt();
				}

				os << String(ms_level+1,'\t')
					 << "<scan num=\"" << spectrum_id << "\" msLevel=\""
					 << ms_level << "\" peaksCount=\""
					 << spec.size() << "\" polarity=\"";
				if (spec.getInstrumentSettings().getPolarity()==IonSource::POSITIVE)
				{
					os << "+";
				}
				else if (spec.getInstrumentSettings().getPolarity()==IonSource::NEGATIVE)
				{
					os << "-";
				}
				else
				{
					os << "any";
				}
				
				if (spec.getInstrumentSettings().getScanMode()!=0)
				{
					os << "\" scanType=\"" << cv_terms_[1][spec.getInstrumentSettings().getScanMode()];
				}
				os << "\" retentionTime=\"";
				if (spec.getRT()<0) os << "-";
				os << "PT"<< std::fabs(spec.getRT()) << "S\"";
				if (spec.getInstrumentSettings().getScanWindows().size() > 0)
				{
					os << " startMz=\"" << spec.getInstrumentSettings().getScanWindows()[0].begin << "\" endMz=\"" << spec.getInstrumentSettings().getScanWindows()[0].end << "\"";
				}
				if (spec.getInstrumentSettings().getScanWindows().size() > 1)
				{
					warning(STORE, "The MzXML format can store only one scan window for each scan. Only the first one is stored!");
				}
				os << ">\n";
	
				const PrecursorPeakType& peak = spec.getPrecursorPeak();
				if (peak!= PrecursorPeakType())
				{
					os << String(ms_level+2,'\t') << "<precursorMz precursorIntensity=\""
						 << peak.getIntensity();
					if (peak.getCharge()!=0) os << "\" precursorCharge=\"" << peak.getCharge();
					os << "\">" << peak.getPosition()[0] << "</precursorMz>\n";
				}
	
				if (spec.size() > 0)
				{
					os << String(ms_level+2,'\t') << "<peaks precision=\"32\"" << " byteOrder=\"network\" pairOrder=\"m/z-int\">";
					
					//std::cout << "Writing scan " << s << std::endl;
					std::vector<Real> tmp;
					for (UInt i=0; i<spec.size(); i++)
					{
						tmp.push_back(spec[i].getMZ());
						tmp.push_back(spec[i].getIntensity());
					}
					
					std::string encoded;
					decoder_.encode(tmp, Base64::BYTEORDER_BIGENDIAN, encoded);
					os << encoded << "</peaks>\n";
				}
				else
				{
					os << String(ms_level+2,'\t') << "<peaks precision=\"32\"" << " byteOrder=\"network\" pairOrder=\"m/z-int\" xsi:nil=\"1\"/>\n";
				}
				
				writeUserParam_(os,spec,ms_level+2);
				if (spec.getComment() != "")
				{
					os << String(ms_level+2,'\t') << "<comment>" << spec.getComment() << "</comment>\n";
				}
				
				//check MS level of next scan and close scans (scans can be nested)
				UInt next_ms_level = 0;
				if (s < cexp_->size()-1)
				{
					next_ms_level = ((*cexp_)[s+1]).getMSLevel();
				}
				//std::cout << "scan: " << s << " this: " << ms_level << " next: " << next_ms_level << std::endl;
				if (next_ms_level <= ms_level)
				{
					for (UInt i = 0; i<= ms_level-next_ms_level && !open_scans.empty(); ++i)
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
			spec_write_counter_ = 1;
		}

	} // namespace Internal

} // namespace OpenMS

#endif
