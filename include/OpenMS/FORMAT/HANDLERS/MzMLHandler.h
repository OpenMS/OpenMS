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
			Do not use this class. It is only needed in MzMLFile.
			
			@todo Parse 16 bit encoded data? (Marc)
			@todo Implement everything else... (Marc)
			
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
					peak_count_(0),
					meta_id_descs_(),
					decoder_(),
					spec_write_counter_(1),
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
      MzMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger)
				: XMLHandler(filename, version),
					exp_(0),
					cexp_(&exp),
					peak_count_(0),
					meta_id_descs_(),
					decoder_(),
					spec_write_counter_(1),
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
			
			void setOptions(const PeakFileOptions& opt) { options_ = opt; }

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
			std::vector<String> binary_array_names_;
			//@}

			/// Decoder/Encoder for Base64-data in MzML
			Base64 decoder_;

			/// spectrum counter (needed because spectra without peaks are not written)
			UInt spec_write_counter_;
								
			/// Flag that indicates wether this spectrum should be skipped (due to options)
			bool skip_spectrum_;
			
			/// Progress logger
			const ProgressLogger& logger_;

			/// fills the current spectrum with peaks and meta data
			void fillData_();

			/** 
				@brief read attributes of MzML's cvParamType
	
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:1000001" name="@p name" value="@p value"/&gt;
				@p name and sometimes @p value are defined in the MzML ontology.
			*/
			void cvParam_(const String& name, const String& value);
			
		};



		//--------------------------------------------------------------------------------

		template <typename MapType>
		void MzMLHandler<MapType>::characters(const XMLCh* const chars, unsigned int /*length*/)
		{
			// skip current spectrum
			if (skip_spectrum_) return;
			
			char* transcoded_chars = sm_.convert(chars);
			
			String& current_tag = open_tags_.back();
			String& parent_tag = *(open_tags_.end()-2);
			
			if (current_tag == "binary")
			{
				//chars may be split to several chunks => concatenate them
				data_to_decode_.back() += transcoded_chars;
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
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");
			static const XMLCh* s_count = xercesc::XMLString::transcode("count");
			static const XMLCh* s_encodedlength = xercesc::XMLString::transcode("encodedLength");

			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			//std::cout << "Start: '" << tag << "'" << std::endl;
			
			//do nothing as until a new spectrum is reached
			if (tag!="spectrum" && skip_spectrum_) return;
			
			if (tag=="cvParam")
			{
				String accession = attributeAsString_(attributes, s_accession);
				String value = "";
				optionalAttributeAsString_(value, attributes, s_value);
				cvParam_(accession, value);
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
		  	logger_.startProgress(0,count,"loading mzML file");
		  	//std::cout << Date::now() << " done" << std::endl;
			}
			else if (tag=="binaryDataArray")
			{
				peak_count_ = attributeAsInt_(attributes, s_encodedlength);
				//spec_.getContainer().reserve(peak_count_); ?
				data_to_decode_.resize(data_to_decode_.size()+1);
			}
			//std::cout << "end startelement" << std::endl;
		}


		template <typename MapType>
		void MzMLHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static UInt scan_count = 0;
			
			static const XMLCh* s_spectrum = xercesc::XMLString::transcode("spectrum");
			static const XMLCh* s_mzml = xercesc::XMLString::transcode("mzML");
			
			open_tags_.pop_back();			
			//std::cout << "End: '" << sm_.convert(qname) << "'" << std::endl;
			
			if(equal_(qname,s_spectrum))
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
				binary_array_names_.clear();
				meta_id_descs_.clear();
			}
			else if(equal_(qname,s_mzml))
			{
				logger_.endProgress();
				scan_count = 0;
			}

			sm_.clear();
		}

		template <typename MapType>
		void MzMLHandler<MapType>::cvParam_(const String& accession, const String& value)
		{
			String error = "";

			String& parent_tag = *(open_tags_.end()-2);
			
			if (parent_tag == "binaryDataArray")
			{
				if (accession == "MS:1000514") // m/z binary array
				{
					binary_array_names_.push_back("m/z");
				}
				else if (accession == "MS:1000515") // intensity binary array
				{
					binary_array_names_.push_back("Intensity");
				}
				
				// what about meta data arrays?? (mailing list)

				else if (accession == "MS:1000521") // 32-bit float
				{
					precisions_.push_back("32");
				}
				else if (accession == "MS:1000523") // 64-bit float
				{
					precisions_.push_back("64");
				}
				
			}
			
//			else
//			{
//				warning(String("Unexpected cvParam: accession=\"") + accession + ", value=\"" + value + "\" in tag " + parent_tag);
//			}
//
//			if (error != "")
//			{
//				warning(String("Invalid cvParam: accession=\"") + accession + ", value=\"" + value + "\" in " + error);
//			}
			//std::cout << "End of MzMLHander::cvParam_" << std::endl;
		}

		template <typename MapType>
		void MzMLHandler<MapType>::fillData_()
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
					//std::cout << "nr. " << i << ": decoding as high-precision little endian" << std::endl;
					decoder_.decode(data_to_decode_[i], Base64::BYTEORDER_LITTLEENDIAN, decoded_double);
					// push_back the decoded double data - and an empty one into
					// the dingle-precision vector, so that we don't mess up the index
					//std::cout << "list size: " << decoded_double.size() << std::endl;
					decoded_double_list_.push_back(decoded_double);
					decoded_list_.push_back(std::vector<float>());
				}
				else if (precisions_[i]=="32") // precision 32 Bit
				{
					//std::cout << "nr. " << i << ": decoding as low-precision little endian" << std::endl;
					decoder_.decode(data_to_decode_[i], Base64::BYTEORDER_LITTLEENDIAN, decoded);
					//std::cout << "list size: " << decoded.size() << std::endl;
					decoded_list_.push_back(decoded);
					decoded_double_list_.push_back(std::vector<double>());
				}
				else // 16 bit
				{
					
				}
			}

			// this works only if MapType::PeakType is at least DPeak
			{
//				//reserve space for meta data arrays (peak count)
//				for (UInt i=0;i<spec_.getMetaDataArrays().size();++i)
//				{
//					spec_.getMetaDataArrays()[i].reserve(peak_count_);
//				}
				
				//push_back the peaks into the container				
				DoubleReal mz, intensity;
				for (UInt n = 0 ; n < peak_count_ ; n++)
				{
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(mz)))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(intensity))))
					{
						PeakType tmp;
						for (UInt i = 0; i < binary_array_names_.size(); i++)
						{
							if (binary_array_names_[i] == "m/z")
							{
								if (precisions_[i] == "64")
								{
									mz = decoded_double_list_[i][n];
								}
								else if (precisions_[i] == "32")
								{
									mz = decoded_list_[i][n];
								}
							}
							else if (binary_array_names_[i] == "Intensity")
							{
								if (precisions_[i] == "64")
								{
									intensity = decoded_double_list_[i][n];
								}
								else if (precisions_[i] == "32")
								{
									intensity = decoded_list_[i][n];
								}
							}
							else
							{
								// meta data arrays...
							}
						}
						tmp.setPosition(mz);
						tmp.setIntensity(intensity);
						
						spec_.push_back(tmp);
					}
				}
			}
		}
	} // namespace Internal
} // namespace OpenMS

#endif
