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

//TODO:
// - units
// - userParam
// - add to automatic tmp file validation in tests
// - TOPP tests for FileInfo, FileMerger, FileConverter
// - OpenMS test: 2 mal hintereinander laden
// - writing files

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
			/// Flag that indicates that we're inside a spectum (in contrast to a chromatogram)
			bool in_spectrum_list_;
			/// Id of the current list. Used for referencable param group, source file, sample, software, ...
			String current_id_;
			/// The referencable param groups: id => array (accession, value)
			Map<String, std::vector< std::pair<String,String> > > ref_param_;
			/// The source files: id => SourceFile
			Map<String, SourceFile> source_files_;
			/// The sample list: id => Sample
			Map<String, Sample> samples_;
			/// The software list: id => Software
			Map<String, Software> software_;
			/// The data processing lust: id => Software
			Map<String, Software> processing_;
			//@}

			/// Decoder/Encoder for Base64-data in MzML
			Base64 decoder_;
			
			/// Progress logger
			const ProgressLogger& logger_;

			/// Fills the current spectrum with peaks and meta data
			void fillData_();			

			/// Handles CV terms
			void handleCVParam_(const String& parent_tag, const String& accession, const String& value);

			/// Handles user terms
			void handleUserParam_(const String& parent_tag, const String& name, const String& type, const String& value);
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
			static const XMLCh* s_default_array_length = xercesc::XMLString::transcode("defaultArrayLength");
			static const XMLCh* s_array_length = xercesc::XMLString::transcode("arrayLength");
			static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");
			static const XMLCh* s_name = xercesc::XMLString::transcode("name");
			static const XMLCh* s_type = xercesc::XMLString::transcode("type");
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			static const XMLCh* s_id = xercesc::XMLString::transcode("id");
			static const XMLCh* s_ref = xercesc::XMLString::transcode("ref");
			static const XMLCh* s_number = xercesc::XMLString::transcode("number");
			static const XMLCh* s_version = xercesc::XMLString::transcode("version");
			static const XMLCh* s_location = xercesc::XMLString::transcode("location");
			static const XMLCh* s_sample_ref = xercesc::XMLString::transcode("sampleRef");
			static const XMLCh* s_software_ref = xercesc::XMLString::transcode("softwareRef");

			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			
			//determine parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);


			//determine the parent tag of the parent tag
			String parent_parent_tag;
			if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);
			
			if (tag=="spectrum")
			{
				spec_ = SpectrumType();
				default_array_length_ = attributeAsInt_(attributes, s_default_array_length);
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
				optionalAttributeAsInt_(array_length, attributes, s_array_length);
				data_.back().size = array_length;
			}
			else if (tag=="cvParam")
			{
				String value = "";
				optionalAttributeAsString_(value, attributes, s_value);
				handleCVParam_(parent_tag, attributeAsString_(attributes, s_accession), value);
			}
			else if (tag=="userParam")
			{
				String type = "";
				optionalAttributeAsString_(type, attributes, s_type);
				String value = "";
				optionalAttributeAsString_(value, attributes, s_value);
				handleUserParam_(parent_tag, attributeAsString_(attributes, s_name), type, value);
			}
			else if (tag=="referenceableParamGroup")
			{
				current_id_ = attributeAsString_(attributes, s_id);
			}
			else if (tag=="sourceFile")
			{
				current_id_ = attributeAsString_(attributes, s_id);
				source_files_[current_id_].setNameOfFile(attributeAsString_(attributes, s_name));
				source_files_[current_id_].setPathToFile(attributeAsString_(attributes, s_location));
			}
			else if (tag=="referenceableParamGroupRef")
			{
				//call handleCVParam_ with the parent tag for each parameter in the group
				String ref = attributeAsString_(attributes, s_ref);
				for (UInt i=0; i<ref_param_[ref].size(); ++i)
				{
					handleCVParam_(parent_tag,ref_param_[ref][i].first,ref_param_[ref][i].second);
				}
			}
			else if (tag=="acquisition")
			{
				Acquisition tmp;
				tmp.setNumber(attributeAsInt_(attributes, s_number));
				spec_.getAcquisitionInfo().push_back(tmp);
			}
			else if (tag=="mzML")
			{
				//check file version against schema version
				String file_version = attributeAsString_(attributes, s_version);
				if (file_version.toDouble()>version_.toDouble())
				{
					warning("The XML file (" + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
				}
			}
			else if (tag=="contact")
			{
				exp_->getContacts().push_back(ContactPerson());
			}
			//TODO Acquisition, Precursor and Spectrum can have a source file too
			else if(tag=="sourceFileRef" && parent_tag=="sourceFileRefList" && parent_parent_tag=="run")
			{
				//EXTEND Store more than one source file. Currently only the last file is stored (ExperimentalSettings)
				exp_->getSourceFile() = source_files_[ attributeAsString_(attributes, s_ref)];
			}
			else if (tag=="sample")
			{
				current_id_ = attributeAsString_(attributes, s_id);
				String name;
				if (optionalAttributeAsString_(name, attributes, s_name))
				{
					samples_[current_id_].setName(name);
				}
			}
			else if (tag=="run")
			{
				String sample_ref;
				if (optionalAttributeAsString_(sample_ref, attributes, s_sample_ref))
				{
					exp_->setSample(samples_[sample_ref]);
				}
			}
			else if (tag=="software")
			{
				current_id_ = attributeAsString_(attributes, s_id);
			}
			else if (tag=="softwareParam")
			{
				//Using an enum for software names is wrong in my (Marc) opinion, so we simply store the name as string
				software_[current_id_].setName(attributeAsString_(attributes, s_name));
				software_[current_id_].setVersion(attributeAsString_(attributes, s_version));
			}
			else if (tag=="dataProcessing")
			{
				current_id_ = attributeAsString_(attributes, s_id);
				processing_[current_id_] = software_[attributeAsString_(attributes, s_software_ref)];
			}
		}

		template <typename MapType>
		void MzMLHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static UInt scan_count=0;
			
			static const XMLCh* s_spectrum = xercesc::XMLString::transcode("spectrum");
			static const XMLCh* s_spectrum_list = xercesc::XMLString::transcode("spectrumList");
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
			else if(equal_(qname,s_spectrum_list))
			{
				in_spectrum_list_ = false;
			}
			else if(equal_(qname,s_mzml))
			{
				logger_.endProgress();
				scan_count = 0;
				ref_param_.clear();
				current_id_ = "";
				source_files_.clear();
				samples_.clear();
				software_.clear();
				processing_.clear();
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
		} //fillData_
		
		template <typename MapType>
		void MzMLHandler<MapType>::handleCVParam_(const String& parent_tag, const String& accession, const String& value)
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
					//TODO Does this really belong here, or should it be under "spectrumDescription"
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
			//------------------------- referenceableParamGroup ----------------------------
			else if(parent_tag=="referenceableParamGroup")
			{
				ref_param_[current_id_].push_back(std::make_pair(accession,value));
			}
			//------------------------- selectedIon ----------------------------
			else if(parent_tag=="selectedIon")
			{
				//EXTEND parse and store more than one precursor (isolationWindow,selectedIon,activation)
				if (accession=="MS:1000040") //m/z
				{
					spec_.getPrecursorPeak().getPosition()[0] = value.toDouble();
				}
				else if (accession=="MS:1000041") //charge state
				{
					spec_.getPrecursorPeak().setCharge(value.toInt());
				}
				else if (accession=="MS:1000042") //intensity
				{
					spec_.getPrecursorPeak().setIntensity(value.toDouble());
				}
				else if (accession=="MS:1000633") //possible charge state
				{
					//EXTEND store possible charge states as well
				}
			}
			//------------------------- activation ----------------------------
			else if(parent_tag=="activation")
			{
				if (accession=="MS:1000245") //charge stripping
				{
					//EXTEND 
				}
				else if (accession=="MS:1000246") //delayed extraction
				{
					//EXTEND 
				}
				else if (accession=="MS:1000045") //collision energy
				{
					//EXTEND 
				}
				else if (accession=="MS:1000412") //buffer gas
				{
					//EXTEND 
				}
				else if (accession=="MS:1000419") //collision gas
				{
					//EXTEND 
				}
				else if (accession=="MS:1000509") //activation energy
				{
					spec_.getPrecursor().setActivationEnergy(value.toDouble());
				}
				else if (accession=="MS:1000133") //collision-induced dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::CID);
				}
				else if (accession=="MS:1000134") //plasma desorption
				{
					spec_.getPrecursor().setActivationMethod(Precursor::PD);
				}
				else if (accession=="MS:1000135") //post-source decay
				{
					spec_.getPrecursor().setActivationMethod(Precursor::PSD);
				}
				else if (accession=="MS:1000136") //surface-induced dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::SID);
				}
				else if (accession=="MS:1000242") //blackbody infrared radiative dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::BIRD);
				}
				else if (accession=="MS:1000250") //electron capture dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::ECD);
				}
				else if (accession=="MS:1000262") //infrared multiphoton dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::IMD);
				}
				else if (accession=="MS:1000282") //sustained off-resonance irradiation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::SORI);
				}
				else if (accession=="MS:1000422") //high-energy collision-induced dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::HCID);
				}
				else if (accession=="MS:1000433") //low-energy collision-induced dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::LCID);
				}
				else if (accession=="MS:1000435") //photodissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::PHD);
				}
				else if (accession=="MS:1000598") //electron transfer dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::ETD);
				}
				else if (accession=="MS:1000599") //pulsed q dissociation
				{
					spec_.getPrecursor().setActivationMethod(Precursor::PQD);
				}
			}
			//------------------------- acquisitionList ----------------------------
			else if(parent_tag=="acquisitionList")
			{
				if (accession=="MS:1000571") //sum of spectra
				{
					spec_.getAcquisitionInfo().setMethodOfCombination("sum");
				}
				else if (accession=="MS:1000573") //median of spectra
				{
					spec_.getAcquisitionInfo().setMethodOfCombination("median");
				}
				else if (accession=="MS:1000575") //mean of spectra
				{
					spec_.getAcquisitionInfo().setMethodOfCombination("mean");
				}
			}
			//------------------------- acquisition ----------------------------
			else if (parent_tag=="acquisition")
			{
				//EXTEND Each acquisition can have all attributes like a scan (children of MS:1000503)
			}
			//------------------------- contact ----------------------------
			else if (parent_tag=="contact")
			{
				if (accession=="MS:1000586") //contact name
				{
					exp_->getContacts().back().setName(value);
				}
				else if (accession=="MS:1000587") //contact address
				{
					//EXTEND
				}
				else if (accession=="MS:1000588") //contact URL
				{
					//EXTEND
				}
				else if (accession=="MS:1000589") //contact email
				{
					exp_->getContacts().back().setEmail(value);
				}
				else if (accession=="MS:1000590") //contact organization
				{
					exp_->getContacts().back().setInstitution(value);
				}
			}
			//------------------------- sourceFile ----------------------------
			else if (parent_tag=="sourceFile")
			{
				if (accession=="MS:1000569") //SHA-1 checksum
				{
					source_files_[current_id_].setSha1(value);
				}
				else if (accession=="MS:1000568") //MD5 checksum
				{
					//EXTEND
				}
				else if (accession=="MS:1000561") //data file checksum type
				{
					//EXTEND
				}
				else if (accession=="MS:1000562") //wiff file
				{
					//EXTEND
				}
				else if (accession=="MS:1000563") //Xcalibur RAW file
				{
					//EXTEND
				}
				else if (accession=="MS:1000564") //mzData file
				{
					//EXTEND
				}
				else if (accession=="MS:1000565") //pkl file
				{
					//EXTEND
				}
				else if (accession=="MS:1000566") //mzXML file
				{
					//EXTEND
				}
				else if (accession=="MS:1000567") //yep file
				{
					//EXTEND
				}
				else if (accession=="MS:1000584") //mzML file
				{
					//EXTEND
				}
				else if (accession=="MS:1000613") //dta file
				{
					//EXTEND
				}
				else if (accession=="MS:1000614") //ProteinLynx Global Server mass spectrum XML file
				{
					//EXTEND
				}
				else if (accession=="MS:1000526") //MassLynx raw format
				{
					//EXTEND
				}
			}
			//------------------------- sample ----------------------------
			else if (parent_tag=="sample")
			{
				if (accession=="MS:1000004") //sample mass
				{
					samples_[current_id_].setMass(value.toDouble());
				}
				else if (accession=="MS:1000001") //sample number
				{
					samples_[current_id_].setNumber(value);
				}
				else if (accession=="MS:1000005") //sample volume
				{
					samples_[current_id_].setVolume(value.toDouble());
				}
				else if (accession=="MS:1000006") //sample concentration
				{
					samples_[current_id_].setConcentration(value.toDouble());
				}
				else if (accession=="MS:1000053") //sample batch
				{
					//EXTEND
				}
			}
		}//handleCVParam_



		template <typename MapType>
		void MzMLHandler<MapType>::handleUserParam_(const String& /*parent_tag*/, const String& /*name*/, const String& /*type*/, const String& /*value*/)
		{
			//TODO
		}//handleUserParam_
				
	} // namespace Internal
} // namespace OpenMS

#endif
