// -*- Mode: C++; tab-width: 2; -*-
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

#ifndef OPENMS_FORMAT_HANDLERS_MZMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZMLHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
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

#include <iostream>

namespace OpenMS
{
	namespace Internal
	{

		/**
			@brief XML handler for MzMLFile
			
			MapType has to be a MSExperiment or have the same interface.

			@note Do not use this class. It is only needed in MzMLFile.
		*/
		template <typename MapType>
		class MzMLHandler
			: public XMLHandler
		{
		 public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzMLHandler(MapType& exp, const String& filename, const String& version, ProgressLogger& logger)
				: XMLHandler(filename, version),
					exp_(&exp),
					cexp_(0),
					options_(),
					spec_(),
					data_(),
					default_array_length_(0),
					in_spectrum_list_(false),
					decoder_(),
					logger_(logger)
	  	{
			}

      /// Constructor for a read-only handler
      MzMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger)
				: XMLHandler(filename, version),
					exp_(0),
					cexp_(&exp),
					options_(),
					spec_(),
					data_(),
					default_array_length_(0),
					in_spectrum_list_(false),
					decoder_(),
					logger_(logger)
  		{
			}

      /// Destructor
      virtual ~MzMLHandler()
      {
      }
      //@}


			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, const unsigned int length);
			
			void setOptions(const PeakFileOptions& opt)
			{
				options_ = opt;
			}

		 protected:
      
      /// Peak type
      typedef typename MapType::PeakType PeakType;
      /// Spectrum type
			typedef MSSpectrum<PeakType, std::allocator<PeakType> > SpectrumType;
			
			/// Spectrum representation
			struct BinaryData
			{
				String base64;
				String precision;
				String name;
				UInt size;
				String compression;
				std::vector<Real> decoded_32;
				std::vector<DoubleReal> decoded_64;
			};
			
			/// map pointer for reading
			MapType* exp_;
			/// map pointer for writing
			const MapType* cexp_;
			
			///Options that can be set for loading/storing
			PeakFileOptions options_;
		
			/**@name temporary datastructures to hold parsed data */
			//@{
			/// The current spectrum
			SpectrumType spec_;
			/// The spectrum data
			std::vector<BinaryData> data_;
			/// The default number of peaks in the current spectrum
			UInt default_array_length_;
			/// flag that indicates that we're inside a spectum (in contrast to a chromatogram)
			bool in_spectrum_list_;
			//@}

			/// Decoder/Encoder for Base64-data in MzML
			Base64 decoder_;
			
			/// Progress logger
			const ProgressLogger& logger_;

			/// Fills the current spectrum with peaks and meta data
			void fillData_();			

			/// Handles CV terms
			void handleCV_(const String& parent_tag, const String& accession, const String& value);
		};



		//--------------------------------------------------------------------------------

		template <typename MapType>
		void MzMLHandler<MapType>::characters(const XMLCh* const chars, unsigned int /*length*/)
		{
			char* transcoded_chars = sm_.convert(chars);
			
			String& current_tag = open_tags_.back();
			
			if (current_tag == "binary" && in_spectrum_list_)
			{
				//chars may be split to several chunks => concatenate them
				data_.back().base64 += transcoded_chars;
			}
			else if (current_tag == "offset" || current_tag == "indexListOffset" || current_tag == "fileChecksum" || current_tag == "binary")
			{
				//do nothing for
				// - index
				// - checksum
				// - binary chromatogram data
			}
			else
			{
				String transcoded_chars2 = transcoded_chars;
				transcoded_chars2.trim();
				if (transcoded_chars2!="") warning(String("Unhandled character content in tag '") + current_tag + "': " + transcoded_chars2);
			}
		}

		template <typename MapType>
		void MzMLHandler<MapType>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
		{			
			static const XMLCh* s_count = xercesc::XMLString::transcode("count");
			static const XMLCh* s_defaultarraylength = xercesc::XMLString::transcode("defaultArrayLength");
			static const XMLCh* s_arraylength = xercesc::XMLString::transcode("arrayLength");
			static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			
			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			
			//determine the parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
			
			if (tag=="spectrum")
			{
				spec_ = SpectrumType();
				default_array_length_ = attributeAsInt_(attributes, s_defaultarraylength);
			}
			else if (tag=="spectrumList")
			{
				if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		  	UInt count = attributeAsInt_(attributes, s_count);
		  	exp_->reserve(count);
		  	logger_.startProgress(0,count,"loading mzML file");
				in_spectrum_list_ = true;
			}
			else if (tag=="binaryDataArrayList" && in_spectrum_list_)
			{
				data_.reserve(attributeAsInt_(attributes, s_count));
			}
			else if (tag=="binaryDataArray" && in_spectrum_list_)
			{
				data_.push_back(BinaryData());
				//set array length
				Int array_length = default_array_length_;
				optionalAttributeAsInt_(array_length, attributes, s_arraylength);
				data_.back().size = array_length;
			}
			else if (tag=="cvParam")
			{
				handleCV_(parent_tag, attributeAsString_(attributes, s_accession), attributeAsString_(attributes, s_value));
			}
			else if (tag=="referenceableParamGroup")
			{
				//TODO implement this and use it where it is allowed
			}
		}

		template <typename MapType>
		void MzMLHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static UInt scan_count=0;
			
			static const XMLCh* s_spectrum = xercesc::XMLString::transcode("spectrum");
			static const XMLCh* s_spectrumlist = xercesc::XMLString::transcode("spectrumList");
			static const XMLCh* s_mzml = xercesc::XMLString::transcode("mzML");
			
			open_tags_.pop_back();			
			
			if(equal_(qname,s_spectrum))
			{
				fillData_();
				exp_->push_back(spec_);
				logger_.setProgress(++scan_count);
				data_.clear();
				default_array_length_ = 0;
			}
			else if(equal_(qname,s_spectrumlist))
			{
				in_spectrum_list_ = false;
			}
			else if(equal_(qname,s_mzml))
			{
				logger_.endProgress();
				scan_count = 0;
			}
			
			sm_.clear();
		}

		template <typename MapType>
		void MzMLHandler<MapType>::fillData_()
		{
			//decode all base64 arrays
			for (UInt i=0; i<data_.size(); i++)
			{
				//remove whitespaces from binary data
				//this should not be necessary, but linebreaks inside the base64 data are unfortunately no exception
				data_[i].base64.removeWhitespaces();

				if (data_[i].precision=="64")
				{
					decoder_.decode(data_[i].base64, Base64::BYTEORDER_LITTLEENDIAN, data_[i].decoded_64);
				}
				else if (data_[i].precision=="32")
				{
					decoder_.decode(data_[i].base64, Base64::BYTEORDER_LITTLEENDIAN, data_[i].decoded_32);
				}
			}

			// this works only if MapType::PeakType is at least DPeak
			{
				//look up the precision and the index of the intensity and m/z array
				bool mz_precision_64 = true;
				bool int_precision_64 = true;
				Int mz_index = -1;
				Int int_index = -1;
				for (UInt i=0; i<data_.size(); i++)
				{
					if (data_[i].name=="mz")
					{
						mz_index = i;
						mz_precision_64 = (data_[i].precision=="64");
					}
					if (data_[i].name=="int")
					{
						int_index = i;
						int_precision_64 = (data_[i].precision=="64");
					}
				}
				
				//Abort if no m/z or intensity array is present
				if (int_index==-1 || mz_index==-1)
				{
					//if defaultArrayLength > 0 : warn that no m/z or int arrays is present
					if (default_array_length_!=0)
					{
						warning(String("The m/z or intensity array of spectrum ") + exp_->size() + " is missing and default_array_length_ is " + default_array_length_ + ".");
					}
					return;
				}
				
				//Warn if the decoded data has a differenct size than the the defaultArrayLength
				UInt mz_size = mz_precision_64 ? data_[mz_index].decoded_64.size() : data_[mz_index].decoded_32.size();
				if (default_array_length_!=mz_size)
				{
					warning(String("The base64-decoded m/z array of spectrum ") + exp_->size() + " has the size " + mz_size + ", but it should have the size " + default_array_length_ + " (defaultArrayLength).");
				}
				UInt int_size = int_precision_64 ? data_[int_index].decoded_64.size() : data_[int_index].decoded_32.size();
				if (default_array_length_!=int_size)
				{
					warning(String("The base64-decoded intensity array of spectrum ") + exp_->size() + " has the size " + int_size + ", but it should have the size " + default_array_length_ + " (defaultArrayLength).");
				}
				
				//create meta data arrays if necessary
				if (data_.size()>2)
				{
					//create meta data arrays and assign meta data
					spec_.getMetaDataArrays().resize(data_.size()-2);
					UInt meta_array_index = 0;
					for (UInt i=0; i<data_.size(); i++)
					{
						if (data_[i].name!="mz" && data_[i].name!="int")
						{
							spec_.getMetaDataArrays()[meta_array_index].setName(data_[i].name);
							spec_.getMetaDataArrays()[meta_array_index].reserve(data_[i].size);
							++meta_array_index;
						}
					}
				}
				
				//add the peaks and the meta data to the container (if they pass the restrictions)s
				spec_.reserve(default_array_length_);
				for (UInt n = 0 ; n < default_array_length_ ; n++)
				{
					DoubleReal mz = mz_precision_64 ? data_[mz_index].decoded_64[n] : data_[mz_index].decoded_32[n];
					DoubleReal intensity = int_precision_64 ? data_[int_index].decoded_64[n] : data_[int_index].decoded_32[n];
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(mz)))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(intensity))))
					{
						//add peak
						PeakType tmp;
						tmp.setIntensity(intensity);
						tmp.setPosition(mz);
						spec_.push_back(tmp);

						//add meta data
						UInt meta_array_index = 0;
						for (UInt i=0; i<data_.size(); i++)
						{
							if (n<data_[i].size && data_[i].name!="mz" && data_[i].name!="int")
							{
								DoubleReal value = (data_[i].precision=="64") ? data_[i].decoded_64[n] : data_[i].decoded_32[n];
								spec_.getMetaDataArrays()[meta_array_index].push_back(value);
								++meta_array_index;
							}
						}
					}
				}				
			}
		} //fillData_
		
		template <typename MapType>
		void MzMLHandler<MapType>::handleCV_(const String& parent_tag, const String& accession, const String& value)
		{
			//------------------------- binaryDataArray ----------------------------
			if (parent_tag=="binaryDataArray" && in_spectrum_list_)
			{
				//MS:1000518 ! binary data type
				if (accession=="MS:1000523") //64-bit float
				{
					data_.back().precision = "64";
				}
				else if (accession=="MS:1000522") //64-bit integer
				{
					throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
				else if (accession=="MS:1000521") //32-bit float
				{
					data_.back().precision = "32";
				}
				else if (accession=="MS:1000519") //32-bit integer
				{
					throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
				else if (accession=="MS:1000520") //16-bit float
				{
					throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
				//MS:1000513 ! binary data array
				else if (accession=="MS:1000514")//m/z array
				{
					data_.back().name = "mz";
				}
				else if (accession=="MS:1000515")//intensity array
				{
					data_.back().name = "int";
				}
				else if (accession=="MS:1000516")//charge array
				{
					data_.back().name = "charge";
				}
				else if (accession=="MS:1000517")//signal to noise array
				{
					data_.back().name = "signal to noise";
				}
				//MS:1000572 ! binary data compression type
				else if (accession=="MS:1000574")//zlib compression
				{
					data_.back().compression = "zlib";
					throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
				else if (accession=="MS:1000576")// no compression
				{
					data_.back().compression = "none";
				}
			}
			//------------------------- spectrum ----------------------------
			else if(parent_tag=="spectrum")
			{
				//MS:1000559 ! spectrum type
				if (accession=="MS:1000579") //MS1 spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::FULL);
				}
				else if (accession=="MS:1000580") //MSn spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::PRODUCT);
				}
				else if (accession=="MS:1000581") //CRM spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CRM);
				}
				else if (accession=="MS:1000582") //SIM spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::SIM);
				}
				else if (accession=="MS:1000583") //SRM spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
				}
				else if (accession=="MS:1000511") //ms level
				{
					//UNCLEAR Does this really belong here, or should it be under "spectrumDescription"
					spec_.setMSLevel(value.toInt());
				}
			}
			//------------------------- spectrumDescription ----------------------------
			else if(parent_tag=="spectrumDescription")
			{
				if (accession=="MS:1000127") //centroid mass spectrum
				{
					spec_.setType(SpectrumSettings::PEAKS);
				}
				else if (accession=="MS:1000128") //profile mass spectrum
				{
					spec_.setType(SpectrumSettings::RAWDATA);
				}
			}
			//------------------------- scan ----------------------------
			else if(parent_tag=="scan")
			{
				if (accession=="MS:1000011")//mass resolution
				{
					//EXTEND Currently only stored for each experiment, not for each spectrum => MassAnalyzer
				}
				else if (accession=="MS:1000015")//scan rate
				{
					//EXTEND Currently only stored for each experiment, not for each spectrum => MassAnalyzer
				}
				else if (accession=="MS:1000016")//scan time
				{
					//TODO Handle unit
					spec_.setRT(value.toDouble());
				}
				else if (accession=="MS:1000023")//isolation width
				{
					//EXTEND Currently only stored for each experiment, not for each spectrum => MassAnalyzer
				}
				else if (accession=="MS:1000092")//decreasing m/z scan
				{
					//EXTEND Currently only stored for each experiment, not for each spectrum => MassAnalyzer
				}
				else if (accession=="MS:1000093")//increasing m/z scan
				{
					//EXTEND Currently only stored for each experiment, not for each spectrum => MassAnalyzer
				}
				else if (accession=="MS:1000094")//scan law: exponential
				{
					//EXTEND Currently only stored for each experiment, not for each spectrum => MassAnalyzer
				}
				else if (accession=="MS:1000095")//scan law: linear
				{
					//EXTEND Currently only stored for each experiment, not for each spectrum => MassAnalyzer
				}
				else if (accession=="MS:1000096")//scan law: quadratic
				{
					//EXTEND Currently only stored for each experiment, not for each spectrum => MassAnalyzer
				}
				else if (accession=="MS:1000129")//negative scan
				{
					spec_.getInstrumentSettings().setPolarity(IonSource::NEGATIVE);
				}
				else if (accession=="MS:1000130")//positive scan
				{
					spec_.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
				}
			}
			//------------------------- scanWindow ----------------------------
			else if(parent_tag=="scanWindow")
			{
				//EXTEND parse and store more than one scan window. Currently only the last window is stored.
				if (accession=="MS:1000501") //scan m/z lower limit
				{
					spec_.getInstrumentSettings().setMzRangeStart(value.toDouble());
				}
				else if (accession=="MS:1000500") //scan m/z upper limit
				{
					spec_.getInstrumentSettings().setMzRangeStop(value.toDouble());
				}
			}
		}//handleCV_
		
	} // namespace Internal
} // namespace OpenMS

#endif
