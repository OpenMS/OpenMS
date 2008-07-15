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
					decoder_(),
					logger_(logger)
	  	{
			}

      /// Constructor for a read-only handler
      MzMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger)
				: XMLHandler(filename, version),
					exp_(0),
					cexp_(&exp),
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
			/// The number of peaks in the current spectrum
			UInt peak_count_;
			//@}

			/// Decoder/Encoder for Base64-data in MzML
			Base64 decoder_;
			
			/// Progress logger
			const ProgressLogger& logger_;

			/// fills the current spectrum with peaks and meta data
			void fillData_();			
		};



		//--------------------------------------------------------------------------------

		template <typename MapType>
		void MzMLHandler<MapType>::characters(const XMLCh* const chars, unsigned int /*length*/)
		{
			char* transcoded_chars = sm_.convert(chars);
			
			String& current_tag = open_tags_.back();
			
			if (current_tag == "binary")
			{
				//chars may be split to several chunks => concatenate them
				data_.back().base64 += transcoded_chars;
			}
			else if (current_tag == "offset" || current_tag == "indexListOffset" || current_tag == "fileChecksum")
			{
				//do nothing for index and checksum
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
			static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			
			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			
			//determine the parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
			
			if (tag=="spectrum")
			{
				spec_.clear();
				peak_count_ = attributeAsInt_(attributes, s_defaultarraylength); 
			}
			else if (tag=="spectrumList")
			{
				if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		  	UInt count = attributeAsInt_(attributes, s_count);
		  	exp_->reserve(count);
		  	logger_.startProgress(0,count,"loading mzML file");
			}
			else if (tag=="binaryDataArray")
			{
				data_.push_back(BinaryData());
			}
			else if (tag=="cvParam")
			{
				String accession = attributeAsString_(attributes, s_accession);
				String value = attributeAsString_(attributes, s_value);
				if (parent_tag=="binaryDataArray")
				{
					if (accession=="MS:1000523") //64-bit float
					{
						data_.back().precision = "64";
					}
					else if (accession=="MS:1000514")//m/z array
					{
						data_.back().name = "mz";
					}
					else if (accession=="MS:1000515")//intensity array
					{
						data_.back().name = "int";
					}
				}
				else if(parent_tag=="spectrum")
				{
					if (accession=="MS:1000511") //ms level
					{
						spec_.setMSLevel(attributeAsInt_(attributes, s_value));
					}
				}
			}
		}


		template <typename MapType>
		void MzMLHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static UInt scan_count=0;
			
			static const XMLCh* s_spectrum = xercesc::XMLString::transcode("spectrum");
			static const XMLCh* s_mzml = xercesc::XMLString::transcode("mzML");
			
			open_tags_.pop_back();			
			
			if(equal_(qname,s_spectrum))
			{
				fillData_();
				exp_->push_back(spec_);
				logger_.setProgress(++scan_count);
				data_.clear();
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
				else if (data_[i].precision=="16")
				{
					//TODO
				}
			}

			// this works only if MapType::PeakType is at least DPeak
			{
				//look up the precision and the index of the intensity and m/z array
				bool mz_precision_64 = true;
				bool int_precision_64 = true;
				UInt mz_index = 0;
				UInt int_index = 0;
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
								
				//push_back the peaks into the container				
				for (UInt n = 0 ; n < peak_count_ ; n++)
				{
					DoubleReal mz = mz_precision_64 ? data_[mz_index].decoded_64[n] : data_[mz_index].decoded_32[n];
					DoubleReal intensity = int_precision_64 ? data_[int_index].decoded_64[n] : data_[int_index].decoded_32[n];
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(mz)))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(intensity))))
					{
						PeakType tmp;
						tmp.setIntensity(intensity);
						tmp.setPosition(mz);
						spec_.push_back(tmp);
					}
				}
			}
			
		} //fillData_
		
	} // namespace Internal
} // namespace OpenMS

#endif
