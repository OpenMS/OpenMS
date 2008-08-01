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

#ifndef OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataExpSettHandler.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <sstream>

namespace OpenMS
{
	namespace Internal
	{

		/**
			@brief XML handler for MzDataFile
			
			MapType has to be a MSExperiment or have the same interface.
			Do not use this class. It is only needed in MzDataFile.
		*/
		template <typename MapType>
		class MzDataHandler
			: public XMLHandler
		{
		 public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzDataHandler(MapType& exp, const String& filename, const String& version, ProgressLogger& logger)
				: XMLHandler(filename, version),
					exp_(&exp),
					cexp_(0),
					peak_count_(0),
					meta_id_descs_(),
					exp_sett_(),
					decoder_(),
					spec_write_counter_(1),
					in_description_(false),
					skip_spectrum_(false),
					logger_(logger)
	  	{
	  		cv_terms_.resize(4);
				// EnergyUnits
				String(";eV;Percent").split(';',cv_terms_[0]);
				// ScanMode
				String(";SelectedIonDetection;MassScan").split(';',cv_terms_[1]);
				// Polarity
				String(";Positive;Negative").split(';',cv_terms_[2]);
				// ActivationMethod
				String(";CID;PSD;PD;SID").split(';',cv_terms_[3]);
			}

      /// Constructor for a read-only handler
      MzDataHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger)
				: XMLHandler(filename, version),
					exp_(0),
					cexp_(&exp),
					peak_count_(0),
					meta_id_descs_(),
					exp_sett_(),
					decoder_(),
					spec_write_counter_(1),
					in_description_(false),
					skip_spectrum_(false),
					logger_(logger)
  		{
	  		cv_terms_.resize(4);
				// EnergyUnits
				String(";eV;Percent").split(';',cv_terms_[0]);
				// ScanMode
				String(";SelectedIonDetection;MassScan").split(';',cv_terms_[1]);
				// Polarity
				String(";Positive;Negative").split(';',cv_terms_[2]);
				// ActivationMethod
				String(";CID;PSD;PD;SID").split(';',cv_terms_[3]);
			}

      /// Destructor
      virtual ~MzDataHandler()
      {
      }
      //@}


			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, const unsigned int length);

  		/// Writes the contents to a stream
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
			
			/// map pointer for reading
			MapType* exp_;
			/// map pointer for writing
			const MapType* cexp_;
			
			///Options that can be set for loading/storing
			PeakFileOptions options_;
		
			/**@name temporary datastructures to hold parsed data */
			//@{
			/// The number of peaks in the current spectrum
			UInt peak_count_;
			/// The current spectrum
			SpectrumType spec_;
			/// An array of pairs MetaInfodescriptions and their ids
			std::vector< std::pair<String, MetaInfoDescription> > meta_id_descs_;
			/// encoded data which is read and has to be decoded
			std::vector<String> data_to_decode_;
			/// floating point numbers which have to be encoded and written
			std::vector<Real> data_to_encode_;
			std::vector<std::vector<Real> > decoded_list_;
			std::vector<std::vector<DoubleReal> > decoded_double_list_;
			std::vector<String> precisions_;
			std::vector<String> endians_;
			/// stream to collect experimental settings
			std::stringstream exp_sett_;
			//@}

			/// Decoder/Encoder for Base64-data in MzData
			Base64 decoder_;

			/// spectrum counter (needed because spectra without peaks are not written)
			UInt spec_write_counter_;
			
			/// Flag that indicates that of the parser is in a description
			bool in_description_;
			
			/// Flag that indicates wether this spectrum should be skipped (due to options)
			bool skip_spectrum_;
			
			/// Progress logger
			const ProgressLogger& logger_;

			/// fills the current spectrum with peaks and meta data
			void fillData_();

			/** 
				@brief read attributes of MzData's cvParamType
	
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:1000001" name="@p name" value="@p value"/&gt;
				@p name and sometimes @p value are defined in the MzData ontology.
			*/
			void cvParam_(const String& name, const String& value);

			/**
				@brief write binary data to stream (first one)
			
				The @p name and @p id are only used if the @p tag is @em supDataArrayBinary or @em supDataArray.
			*/
			inline void writeBinary_(std::ostream& os, UInt size, const String& tag, const String& name="", int id=-1)
			{
				os 	<< "\t\t\t<" << tag;
				if (tag=="supDataArrayBinary" || tag=="supDataArray")
				{
					os << " id=\"" << id << "\"";
				}
				os << ">\n";
				if (tag=="supDataArrayBinary" || tag=="supDataArray")
				{
					os << "\t\t\t\t<arrayName>" << name << "</arrayName>\n";
				}

				std::string str;				
				decoder_.encode(data_to_encode_, Base64::BYTEORDER_LITTLEENDIAN, str);
				data_to_encode_.clear();
				os << "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\""
					 << size << "\">"
					 << str
					 << "</data>\n\t\t\t</" << tag << ">\n";
			}
			
		};



		//--------------------------------------------------------------------------------

		template <typename MapType>
		void MzDataHandler<MapType>::characters(const XMLCh* const chars, unsigned int /*length*/)
		{
			// skip current spectrum
			if (skip_spectrum_) return;
			
			char* transcoded_chars = sm_.convert(chars);
			
			// collect Experimental Settings
			if (in_description_)
			{
				exp_sett_ << transcoded_chars;
				return;
			}
			
			//current tag
			const String& current_tag = open_tags_.back();
			
			//determine the parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
			
			if (current_tag == "comments" && parent_tag=="spectrumDesc")
			{
				spec_.setComment( transcoded_chars );
			}
			else if (current_tag == "data")
			{
				//chars may be split to several chunks => concatenate them
				data_to_decode_.back() += transcoded_chars;
			}
			else if (current_tag == "arrayName" && parent_tag=="supDataArrayBinary")
			{
				spec_.getMetaDataArrays().back().setName(transcoded_chars);
			}
			else if (current_tag == "nameOfFile" && parent_tag == "supSourceFile")
			{
				meta_id_descs_.back().second.getSourceFile().setNameOfFile( transcoded_chars );
			}
			else if (current_tag == "pathToFile" && parent_tag == "supSourceFile")
			{
				meta_id_descs_.back().second.getSourceFile().setPathToFile( transcoded_chars );
			}
			else if (current_tag == "fileType" && parent_tag == "supSourceFile")
			{
				meta_id_descs_.back().second.getSourceFile().setFileType( transcoded_chars );
			}
			else
			{
				String trimmed_transcoded_chars = transcoded_chars;
				trimmed_transcoded_chars.trim();
				if (trimmed_transcoded_chars!="")
				{
					warning(String("Unhandled character content in tag '") + current_tag + "': " + trimmed_transcoded_chars);
				}
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
		{
			static const XMLCh* s_name = xercesc::XMLString::transcode("name");
			static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			static const XMLCh* s_id = xercesc::XMLString::transcode("id");
			static const XMLCh* s_count = xercesc::XMLString::transcode("count");
			static const XMLCh* s_spectrumtype = xercesc::XMLString::transcode("spectrumType");
			static const XMLCh* s_methodofcombination = xercesc::XMLString::transcode("methodOfCombination");
			static const XMLCh* s_acqnumber = xercesc::XMLString::transcode("acqNumber");
			static const XMLCh* s_mslevel = xercesc::XMLString::transcode("msLevel");
			static const XMLCh* s_mzrangestart = xercesc::XMLString::transcode("mzRangeStart");
			static const XMLCh* s_mzrangestop = xercesc::XMLString::transcode("mzRangeStop");
			static const XMLCh* s_supdataarrayref = xercesc::XMLString::transcode("supDataArrayRef");
			static const XMLCh* s_precision = xercesc::XMLString::transcode("precision");
			static const XMLCh* s_endian = xercesc::XMLString::transcode("endian");
			static const XMLCh* s_length = xercesc::XMLString::transcode("length");
			static const XMLCh* s_comment = xercesc::XMLString::transcode("comment");
			
			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			//std::cout << "Start: '" << tag << "'" << std::endl;
			
			//determine the parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
							
			//do nothing as until a new spectrum is reached
			if (tag!="spectrum" && skip_spectrum_) return;
			
			// collect Experimental Settings
			if (in_description_)
			{
				exp_sett_ << '<' << tag;
				UInt n=attributes.getLength();
				for (UInt i=0; i<n; ++i)
				{
					exp_sett_ << ' ' << sm_.convert(attributes.getQName(i)) << "=\""	<< sm_.convert(attributes.getValue(i)) << '\"';
				}
				exp_sett_ << '>';
				return;
			}
			
			if (tag=="description")
			{
				exp_sett_ << "<description>";
				in_description_ = true; 
			}
			else if (tag=="cvParam")
			{
				String accession = attributeAsString_(attributes, s_accession);
				String value = "";
				optionalAttributeAsString_(value, attributes, s_value);
				cvParam_(accession, value);
			}
			else if (tag=="supDataDesc")
			{
				String comment;
				optionalAttributeAsString_(comment, attributes, s_comment);
				meta_id_descs_.back().second.setComment(comment);
			}
			else if (tag=="userParam")
			{
				String name = attributeAsString_(attributes, s_name);
				String value = "";
				optionalAttributeAsString_(value, attributes, s_value);
				
				if(parent_tag=="spectrumInstrument")
				{
					spec_.getInstrumentSettings().setMetaValue(name, value);
				}
				else if(parent_tag=="acquisition")
				{
					spec_.getAcquisitionInfo().back().setMetaValue(name, value);
				}
				else if (parent_tag=="ionSelection")
				{
					spec_.getPrecursorPeak().setMetaValue(name, value);
				}
				else if (parent_tag=="activation")
				{
					spec_.getPrecursor().setMetaValue(name, value);
				}
				else if (parent_tag=="supDataDesc")
				{
					meta_id_descs_.back().second.setMetaValue(name, value);
				}
				else
				{
					warning("Invalid userParam: name=\"" + name + ", value=\"" + value + "\"");
				}
			}
			else if (tag=="supDataArrayBinary")
			{
				
				//create MetaDataArray
				typename MapType::SpectrumType::MetaDataArray mda;
				//Assign the right MetaInfoDescription ("supDesc" tag)
				String id = attributeAsString_(attributes, s_id);
				for (UInt i=0;i<meta_id_descs_.size(); ++i)
				{
					if (meta_id_descs_[i].first==id)
					{
						mda.MetaInfoDescription::operator=(meta_id_descs_[i].second);
						break;
					}
				}
				//append MetaDataArray
				spec_.getMetaDataArrays().push_back(mda);
			}
			else if (tag=="spectrum")
			{
				spec_ = SpectrumType();
			}
			else if (tag=="spectrumList")
			{
				if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		  	//std::cout << Date::now() << " Reserving space for spectra" << std::endl;
		  	UInt count = attributeAsInt_(attributes, s_count);
		  	exp_->reserve(count);
		  	logger_.startProgress(0,count,"loading mzData file");
		  	//std::cout << Date::now() << " done" << std::endl;
			}
			else if (tag=="acqSpecification")
			{
				String tmp_type = attributeAsString_(attributes, s_spectrumtype);
				if  (tmp_type == "discrete")
				{
					spec_.setType(SpectrumSettings::PEAKS);
				}
				else if (tmp_type == "continuous")
				{
					spec_.setType(SpectrumSettings::RAWDATA);
				}
				else
				{
					spec_.setType(SpectrumSettings::UNKNOWN);
					warning(String("Invalid MzData/SpectrumList/Spectrum/SpectrumDescription/SpectrumSettings/acqSpecification/SpectrumType '") + tmp_type + "'.");
				}
				
				spec_.getAcquisitionInfo().setMethodOfCombination(attributeAsString_(attributes, s_methodofcombination));
			}
			else if (tag=="acquisition")
			{
				spec_.getAcquisitionInfo().insert(spec_.getAcquisitionInfo().end(), Acquisition());
				spec_.getAcquisitionInfo().back().setNumber(attributeAsInt_(attributes, s_acqnumber));
			}
			else if (tag=="spectrumInstrument" || tag=="acqInstrument")
			{
				spec_.setMSLevel(attributeAsInt_(attributes, s_mslevel));
				DoubleReal start=0.0, stop=0.0;
				optionalAttributeAsDouble_(start, attributes, s_mzrangestart);
				optionalAttributeAsDouble_(stop, attributes, s_mzrangestop);
				spec_.getInstrumentSettings().setMzRangeStart(start);
				spec_.getInstrumentSettings().setMzRangeStop(stop);
				
				if (options_.hasMSLevels() && !options_.containsMSLevel(spec_.getMSLevel()))
				{
					skip_spectrum_ = true;
				}
			}
			else if (tag=="supDesc")
			{
				meta_id_descs_.push_back(std::make_pair(attributeAsString_(attributes, s_supdataarrayref),MetaInfoDescription()));
			}
			else if (tag=="data")
			{
				// store precision for later
				precisions_.push_back(attributeAsString_(attributes, s_precision));
				endians_.push_back(attributeAsString_(attributes, s_endian));

				//reserve enough space in spectrum
				if (parent_tag=="mzArrayBinary")
				{
					peak_count_ = attributeAsInt_(attributes, s_length);
					spec_.getContainer().reserve(peak_count_);
				}			
			}
			else if (tag=="mzArrayBinary")
			{
				data_to_decode_.resize(data_to_decode_.size()+1);
			}
			else if (tag=="intenArrayBinary")
			{
					data_to_decode_.resize(data_to_decode_.size()+1);
			}
			else if (tag=="arrayName")
			{
					// Note: name is set in closing tag as it is CDATA
					data_to_decode_.resize(data_to_decode_.size()+1);
			}
			//std::cout << "end startelement" << std::endl;
		}


		template <typename MapType>
		void MzDataHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static UInt scan_count = 0;
			
			static const XMLCh* s_description = xercesc::XMLString::transcode("description");
			static const XMLCh* s_spectrum = xercesc::XMLString::transcode("spectrum");
			static const XMLCh* s_mzdata = xercesc::XMLString::transcode("mzData");
			
			open_tags_.pop_back();			
			//std::cout << "End: '" << sm_.convert(qname) << "'" << std::endl;
			
			if (in_description_)	// collect Experimental Settings
			{
				exp_sett_ << "</" << sm_.convert(qname) << ">\n";
				if (!equal_(qname,s_description)) return;
			}
			
			if(equal_(qname,s_description))
			{
				// initialize parser
				xercesc::XMLPlatformUtils::Initialize();
				xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
				parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
				parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);
				
				MzDataExpSettHandler handler(exp_->getExperimentalSettings(),file_);
				handler.resetErrors();
				parser->setContentHandler(&handler);
				parser->setErrorHandler(&handler);
				
				String tmp(exp_sett_.str().c_str());
				xercesc::MemBufInputSource source((const XMLByte*)(tmp.c_str()), tmp.size(), "dummy");
				parser->parse(source);
				delete(parser);
				
				in_description_ = false;
			}
			else if(equal_(qname,s_spectrum))
			{
				if (!skip_spectrum_)
				{
					fillData_();
					exp_->push_back(spec_);
				}
				skip_spectrum_ = false;
				logger_.setProgress(++scan_count);
				decoded_list_.clear();
				decoded_double_list_.clear();
				data_to_decode_.clear();
				precisions_.clear();
				endians_.clear();
				meta_id_descs_.clear();
			}
			else if(equal_(qname,s_mzdata))
			{
				logger_.endProgress();
				scan_count = 0;
			}

			sm_.clear();
		}

		template <typename MapType>
		void MzDataHandler<MapType>::cvParam_(const String& accession, const String& value)
		{
			String error = "";

			//determine the parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
			
			if(parent_tag=="spectrumInstrument")
			{
				if (accession=="PSI:1000036") //Scan Mode
				{
					spec_.getInstrumentSettings().setScanMode((InstrumentSettings::ScanMode)cvStringToEnum_(1, value,"scan mode"));
				}
				else if (accession=="PSI:1000038") //Time in minutes
				{
					spec_.setRT(asFloat_(value)*60); //Minutes to seconds
					if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRT())))
					{
						skip_spectrum_=true;
					}
				}
				else if (accession=="PSI:1000039") //Time in seconds
				{
					spec_.setRT(asFloat_(value));
					if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRT())))
					{
						skip_spectrum_=true;
					}
				}
				else if (accession=="PSI:1000037") //Polarity
				{
					spec_.getInstrumentSettings().setPolarity((IonSource::Polarity)cvStringToEnum_(2, value,"polarity") );				
				}
				else 
				{
			  	error = "SpectrumDescription.SpectrumSettings.SpectrumInstrument";
				}
			}
			else if(parent_tag=="ionSelection")
			{
				if (accession=="PSI:1000040") //m/z
				{
					spec_.getPrecursorPeak().getPosition()[0] = asFloat_(value);			
				}
				else if (accession=="PSI:1000041") //Charge
				{
					spec_.getPrecursorPeak().setCharge(asInt_(value));		
				}
				else if (accession=="PSI:1000042") //Intensity
				{
					spec_.getPrecursorPeak().setIntensity(asFloat_(value));		
				}
				else if (accession=="PSI:1000043") //Intensity unit (not really handled)
				{
					spec_.getPrecursorPeak().setMetaValue("#IntensityUnits", value);
				}
				else
				{
					error = "PrecursorList.Precursor.IonSelection.UserParam";
				}
			}
			else if (parent_tag=="activation") 
			{
				if (accession=="PSI:1000044") //Method
				{
					spec_.getPrecursor().setActivationMethod((Precursor::ActivationMethod)cvStringToEnum_(3, value,"activation method"));
				}
				else if (accession=="PSI:1000045") //Energy
				{
					spec_.getPrecursor().setActivationEnergy( asFloat_(value) );
				}
				else if (accession=="PSI:1000046") //Energy unit
				{
					spec_.getPrecursor().setActivationEnergyUnit((Precursor::EnergyUnits)cvStringToEnum_(0, value, "energy unit"));
				}
				else 
				{
					error = "PrecursorList.Precursor.Activation.UserParam";
				}
			}
			else if (parent_tag=="supDataDesc") 
			{
				//no terms defined in ontology
				error = "supDataDesc.UserParam";
			}
			else if (parent_tag=="acquisition")
			{
				 //no terms defined in ontology
				 error = "spectrumDesc.spectrumSettings.acquisitionSpecification.acquisition.UserParam";
			}
			else
			{
				warning(String("Unexpected cvParam: accession=\"") + accession + ", value=\"" + value + "\" in tag " + parent_tag);
			}

			if (error != "")
			{
				warning(String("Invalid cvParam: accession=\"") + accession + ", value=\"" + value + "\" in " + error);
			}
			//std::cout << "End of MzDataHander::cvParam_" << std::endl;
		}

		template <typename MapType>
		void MzDataHandler<MapType>::fillData_()
		{
			std::vector<Real> decoded;
			std::vector<DoubleReal> decoded_double;
			
			// data_to_decode is an encoded spectrum, represented as
			// vector of base64-encoded strings:
			// Each string represents one property (e.g. mzData) and decodes
			// to a vector of property values - one value for every peak in the spectrum.
			for (UInt i=0; i<data_to_decode_.size(); i++)
			{
				//remove whitespaces from binary data
				//this should not be necessary, but linebreaks inside the base64 data are unfortunately no exception
				data_to_decode_[i].removeWhitespaces();
				
				if (precisions_[i]=="64")	// precision 64 Bit
				{
					if (endians_[i]=="big")
					{
						//std::cout << "nr. " << i << ": decoding as high-precision big endian" << std::endl;
						decoder_.decode(data_to_decode_[i], Base64::BYTEORDER_BIGENDIAN, decoded_double);
					}
					else
					{
						//std::cout << "nr. " << i << ": decoding as high-precision little endian" << std::endl;
						decoder_.decode(data_to_decode_[i], Base64::BYTEORDER_LITTLEENDIAN, decoded_double);
					}
					// push_back the decoded double data - and an empty one into
					// the dingle-precision vector, so that we don't mess up the index
					//std::cout << "list size: " << decoded_double.size() << std::endl;
					decoded_double_list_.push_back(decoded_double);
					decoded_list_.push_back(std::vector<float>());
				}
				else
				{											// precision 32 Bit
					if (endians_[i]=="big")
					{
						//std::cout << "nr. " << i << ": decoding as low-precision big endian" << std::endl;
						decoder_.decode(data_to_decode_[i], Base64::BYTEORDER_BIGENDIAN, decoded);
					}
					else
					{
						//std::cout << "nr. " << i << ": decoding as low-precision little endian" << std::endl;
						decoder_.decode(data_to_decode_[i], Base64::BYTEORDER_LITTLEENDIAN, decoded);
					}
					//std::cout << "list size: " << decoded.size() << std::endl;
					decoded_list_.push_back(decoded);
					decoded_double_list_.push_back(std::vector<double>());
				}
			}

			// this works only if MapType::PeakType is at least DPeak
			{
				//store what precision is used for intensity and m/z
				bool mz_precision_64 = true;
				if (precisions_[0]=="32")
				{
					mz_precision_64 = false;
				}
				bool int_precision_64 = true;
				if (precisions_[1]=="32")
				{
					int_precision_64 = false;
				}
				
				//reserve space for meta data arrays (peak count)
				for (UInt i=0;i<spec_.getMetaDataArrays().size();++i)
				{
					spec_.getMetaDataArrays()[i].reserve(peak_count_);
				}
				
				//push_back the peaks into the container				
				for (UInt n = 0 ; n < peak_count_ ; n++)
				{
					DoubleReal mz = mz_precision_64 ? decoded_double_list_[0][n] : decoded_list_[0][n];
					DoubleReal intensity = int_precision_64 ? decoded_double_list_[1][n] : decoded_list_[1][n];
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(mz)))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(intensity))))
					{
						PeakType tmp;
						tmp.setIntensity(intensity);
						tmp.setPosition(mz);
						spec_.push_back(tmp);
						//load data from meta data arrays
						for (UInt i=0;i<spec_.getMetaDataArrays().size();++i)
						{
							spec_.getMetaDataArrays()[i].push_back(precisions_[2+i]=="64" ? decoded_double_list_[2+i][n] : decoded_list_[2+i][n]);
						}
					}
				}
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::writeTo(std::ostream& os)
		{
			logger_.startProgress(0,cexp_->size(),"storing mzData file");
			
			os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
				 << "<mzData version=\"1.05\" accessionNumber=\"\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"http://psidev.sourceforge.net/ms/xml/mzdata/mzdata.xsd\">\n";

			// delegate control to ExperimentalSettings handler
			Internal::MzDataExpSettHandler handler( cexp_->getExperimentalSettings(),"");
			handler.writeTo(os);

			if (cexp_->size()!=0)
			{
				os << "\t<spectrumList count=\"" << cexp_->size() << "\">\n";
				int spectrum_ref = -1;
				for (UInt s=0; s<cexp_->size(); s++)
				{
					logger_.setProgress(s);
					//std::cout << "writing scan" << std::endl;
					
					const SpectrumType& spec = (*cexp_)[s];
	
					os << "\t\t<spectrum id=\"" << spec_write_counter_++ << "\">\n"
						 << "\t\t\t<spectrumDesc>\n"
						 << "\t\t\t\t<spectrumSettings>\n";
	
					if (!spec.getAcquisitionInfo().empty())
					{
						os << "\t\t\t\t\t<acqSpecification spectrumType=\"";
						if (spec.getType()==SpectrumSettings::PEAKS)
						{
							os << "discrete";
						}
						else if (spec.getType()==SpectrumSettings::RAWDATA)
						{
							os << "continuous";
						}
	
						os << "\" methodOfCombination=\""
							 << spec.getAcquisitionInfo().getMethodOfCombination() << "\" count=\""
							 << spec.getAcquisitionInfo().size() << "\">\n";
						for (UInt i=0; i<spec.getAcquisitionInfo().size(); ++i)
						{
							const Acquisition& ac = spec.getAcquisitionInfo()[i];
							os << "\t\t\t\t\t\t<acquisition acqNumber=\"" << ac.getNumber() << "\">\n";
							writeUserParam_(os, ac, 7);
							os << "\t\t\t\t\t\t</acquisition>\n";
						}
						os << "\t\t\t\t\t</acqSpecification>\n";
					}
	
					const InstrumentSettings& iset = spec.getInstrumentSettings();
					os << "\t\t\t\t\t<spectrumInstrument msLevel=\"" << spec.getMSLevel()
						 << "\"";
	
					if (spec.getMSLevel()==1) spectrum_ref = spec_write_counter_-1;
					if (iset.getMzRangeStart() != 0 && iset.getMzRangeStop() != 0)
					{
						os << " mzRangeStart=\""
							 << iset.getMzRangeStart() << "\" mzRangeStop=\""
							 << iset.getMzRangeStop() << "\"";
					}
					os << ">\n";
	
					writeCVS_(os, spec.getInstrumentSettings().getScanMode(), 1, "1000036", "ScanMode",6);
					writeCVS_(os, spec.getInstrumentSettings().getPolarity(), 2, "1000037", "Polarity",6);
					//Retiontion time already in TimeInSeconds
					writeCVS_(os, spec.getRT(), "1000039", "TimeInSeconds",6);
					writeUserParam_(os, spec.getInstrumentSettings(), 6);
					os 	<< "\t\t\t\t\t</spectrumInstrument>\n\t\t\t\t</spectrumSettings>\n";
	
					typedef typename SpectrumType::PrecursorPeakType PrecursorPeak;
					if (spec.getPrecursorPeak() != PrecursorPeak()
							|| spec.getPrecursor() != Precursor())
					{
						os	<< "\t\t\t\t<precursorList count=\"1\">\n"
								<< "\t\t\t\t\t<precursor msLevel=\"2\" spectrumRef=\""
								<< spectrum_ref << "\">\n";
						os << "\t\t\t\t\t\t<ionSelection>\n";
						if (spec.getPrecursorPeak() != PrecursorPeak())
						{
							const PrecursorPeak& peak = spec.getPrecursorPeak();
							writeCVS_(os, peak.getPosition()[0], "1000040", "MassToChargeRatio",7);
							writeCVS_(os, peak.getCharge(), "1000041", "ChargeState",7);
							writeCVS_(os, peak.getIntensity(), "1000042", "Intensity",7);
							if (peak.metaValueExists("#IntensityUnits"))
							{
								writeCVS_(os, String(peak.getMetaValue("#IntensityUnits")), "1000043", "IntensityUnits",7);
							}
							writeUserParam_(os, peak, 7);
						}
						os << "\t\t\t\t\t\t</ionSelection>\n";
						os << "\t\t\t\t\t\t<activation>\n";
						if (spec.getPrecursor() != Precursor())
						{
							const Precursor& prec = spec.getPrecursor();
							writeCVS_(os, prec.getActivationMethod(), 3, "1000044", "Method",7);
							writeCVS_(os, prec.getActivationEnergy(), "1000045", "CollisionEnergy",7);
							writeCVS_(os, prec.getActivationEnergyUnit(), 0,"1000046", "EnergyUnits",7);
							writeUserParam_(os, prec,7);
						}
						os << "\t\t\t\t\t\t</activation>\n";
						os << "\t\t\t\t\t</precursor>\n"
							 << "\t\t\t\t</precursorList>\n";
					}
					os << "\t\t\t</spectrumDesc>\n";
	
					// write the supplementary data?
					if (options_.getWriteSupplementalData())
					{
						//write meta data array descriptions
						for (UInt i=0; i<spec.getMetaDataArrays().size(); ++i)
						{
							const MetaInfoDescription& desc = spec.getMetaDataArrays()[i];
							os << "\t\t\t<supDesc supDataArrayRef=\"" << (i+1) << "\">\n";
							if (!desc.isMetaEmpty())
							{
								os << "\t\t\t\t<supDataDesc";
								if (desc.getComment()!="")
								{
									os << " comment=\"" << desc.getComment() << "\"";
								}
								os << ">\n";
								writeUserParam_(os, desc, 5);
								os << "\t\t\t\t</supDataDesc>\n";
							}
							if (desc.getSourceFile()!=SourceFile())
							{
								os << "\t\t\t\t<supSourceFile>\n"
						 				<< "\t\t\t\t\t<nameOfFile>" << desc.getSourceFile().getNameOfFile()
										<< "</nameOfFile>\n"
						 				<< "\t\t\t\t\t<pathToFile>" << desc.getSourceFile().getPathToFile()
										<< "</pathToFile>\n";
								if (desc.getSourceFile().getFileType()!="")	os << "\t\t\t\t\t<fileType>"
									<< desc.getSourceFile().getFileType()	<< "</fileType>\n";
								os << "\t\t\t\t</supSourceFile>\n";
							}
							os << "\t\t\t</supDesc>\n";
						}
					}
					
					//write m/z and intensity arrays
					data_to_encode_.clear();
					for (UInt i=0; i<spec.size(); i++)
					{
						data_to_encode_.push_back(spec.getContainer()[i].getPosition()[0]);
					}
					
					writeBinary_(os,spec.size(),"mzArrayBinary");
	
					// intensity
					data_to_encode_.clear();
					for (UInt i=0; i<spec.size(); i++)
					{
						data_to_encode_.push_back(spec.getContainer()[i].getIntensity());
					}
					
					writeBinary_(os,spec.size(),"intenArrayBinary");
	
					// write the supplementary data?
					if (options_.getWriteSupplementalData())
					{
						//write supplemental data arrays
						for (UInt i=0; i<spec.getMetaDataArrays().size(); ++i)
						{
							const typename MapType::SpectrumType::MetaDataArray& mda = spec.getMetaDataArrays()[i];
							//check if spectrum and meta data array have the same length
							if (mda.size()!=spec.size())
							{
								error(String("Length of meta data array (index:'")+i+"' name:'"+mda.getName()+"') differs from spectrum length. meta data array: " + mda.size() + " / spectrum: " + spec.size() +" .");
							}
							//encode meta data array
							data_to_encode_.clear();
							for (UInt j=0; j<mda.size(); j++)
							{
								data_to_encode_.push_back (mda[j]);
							}
							//write meta data array
							writeBinary_(os,mda.size(),"supDataArrayBinary",mda.getName(),i+1);
						}
					}
	
					os <<"\t\t</spectrum>\n";
				}
			}
			else
			{
				os << "\t<spectrumList count=\"1\">\n";
				os <<"\t\t<spectrum id=\"1\">\n";
				os <<"\t\t\t<spectrumDesc>\n";
				os <<"\t\t\t\t<spectrumSettings>\n";
				os <<"\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n";
				os <<"\t\t\t\t</spectrumSettings>\n";
				os <<"\t\t\t</spectrumDesc>\n";
				os <<"\t\t\t<mzArrayBinary>\n";
				os <<"\t\t\t\t<data length=\"0\" endian=\"little\" precision=\"32\"></data>\n";
				os <<"\t\t\t</mzArrayBinary>\n";
				os <<"\t\t\t<intenArrayBinary>\n";
				os <<"\t\t\t\t<data length=\"0\" endian=\"little\" precision=\"32\"></data>\n";
				os <<"\t\t\t</intenArrayBinary>\n";
				os <<"\t\t</spectrum>\n";
			}
			os << "\t</spectrumList>\n</mzData>\n";
			
			logger_.endProgress();
		}

	} // namespace Internal

} // namespace OpenMS

#endif
