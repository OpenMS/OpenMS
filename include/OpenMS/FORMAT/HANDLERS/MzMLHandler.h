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

#ifndef OPENMS_FORMAT_HANDLERS_MZMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZMLHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

#include <sstream>
#include <iostream>

//TODO:
// - ExternalSpectrumID/ExternalNativeID (when back)
// - Check CV terms of spectrum and scan (when settled)
// - Check CV terms of spectrum type, scanning method and file content (when settled)
// - Check isolationWindow CV
// - DataProcessing of binaryDataArray
// - Sample: CVs for cellular compartement, source tissue and quality
// - units
// - scanSettingsList 
//
//TODO (PERHAPS):
// - DataProcessing of spectrum and chromatogram
// - InstrumentConfiguration of Acquisiton
//
//MISSING (AND NOT PLANNED):
// - more than one precursor per spectrum (warning if more than one)
// - more than one selected ion per precursor (warning if more than one)
// - scanWindowList for each acquisition separately

// xs:id/xs:idref prefix list
// - sf_ru : sourceFile (run)
// - sf_sp : sourceFile (spectrum)
// - sf_pr : sourceFile (precursor)
// - sf_ac : sourceFile (acquisition)
// - sa    : sample
// - ic    : instrumentConfiguration
// - so_dp : software (data processing)
// - so_in : software (instrument)
// - dp_ru : dataProcessing (run)
// - dp_sp : dataProcessing (spectrum)
// - dp_bi : dataProcessing (binary data array)
// - dp_ch : dataProcessing (chromatogram)

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
					logger_(logger),
					skip_spectrum_(false)
	  	{
	  		cv_.loadFromOBO("psi-ms", File::find("CV/psi-ms.obo"));
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
					logger_(logger),
					skip_spectrum_(false)
  		{
  			cv_.loadFromOBO("psi-ms", File::find("CV/psi-ms.obo"));
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
			
			//Docu in base class
			virtual void writeTo(std::ostream& os);
			
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
				MetaInfo meta;
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
			Map<String, std::vector< SemanticValidator::CVTerm > > ref_param_;
			/// The source files: id => SourceFile
			Map<String, SourceFile> source_files_;
			/// The sample list: id => Sample
			Map<String, Sample> samples_;
			/// The software list: id => Software
			Map<String, Software> software_;
			/// The data processing list: id => Instrument
			Map<String, Instrument> instruments_;
			/// The data processing list: id => Instrument
			Map<String, std::vector<DataProcessing> > processing_;
			//@}

			/// Decoder/Encoder for Base64-data in MzML
			Base64 decoder_;
			
			/// Progress logger
			const ProgressLogger& logger_;
			
			/// Flag that indicates whether this spectrum should be skipped (due to options)
			bool skip_spectrum_;
			
			///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
			ControlledVocabulary cv_;
			
			/// Fills the current spectrum with peaks and meta data
			void fillData_();			

			/// Handles CV terms
			void handleCVParam_(const String& parent_parent_parent_tag, const String& parent_parent_tag, const String& parent_tag, const String& accession, const String& name, const String& value);

			/// Handles user terms
			void handleUserParam_(const String& parent_tag, const String& name, const String& type, const String& value);
			
			/// Writes user terms
			void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent) const;
			
			/// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
			ControlledVocabulary::CVTerm getChildWithName_(const String& parent_accession, const String& name) const;
			
			/// Helper method that writes a software
			void writeSoftware_(std::ostream& os, const String& id, const Software& software);

			/// Helper method that writes a source file
			void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software);
		};

		//--------------------------------------------------------------------------------

		template <typename MapType>
		void MzMLHandler<MapType>::characters(const XMLCh* const chars, unsigned int /*length*/)
		{
			if (skip_spectrum_) return;
				
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
				if (transcoded_chars2!="") warning(LOAD, String("Unhandled character content in tag '") + current_tag + "': " + transcoded_chars2);
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
			static const XMLCh* s_spot_id = xercesc::XMLString::transcode("spotID");
			static const XMLCh* s_ref = xercesc::XMLString::transcode("ref");
			static const XMLCh* s_number = xercesc::XMLString::transcode("number");
			static const XMLCh* s_version = xercesc::XMLString::transcode("version");
			static const XMLCh* s_order = xercesc::XMLString::transcode("order");
			static const XMLCh* s_location = xercesc::XMLString::transcode("location");
			static const XMLCh* s_sample_ref = xercesc::XMLString::transcode("sampleRef");
			static const XMLCh* s_software_ref = xercesc::XMLString::transcode("softwareRef");
			static const XMLCh* s_source_file_ref = xercesc::XMLString::transcode("sourceFileRef");
			static const XMLCh* s_default_instrument_configuration_ref = xercesc::XMLString::transcode("defaultInstrumentConfigurationRef");
			static const XMLCh* s_default_data_processing_ref = xercesc::XMLString::transcode("defaultDataProcessingRef");
			static const XMLCh* s_start_time_stamp = xercesc::XMLString::transcode("startTimeStamp");
			//static const XMLCh* s_external_native_id = xercesc::XMLString::transcode("externalNativeID");
			//static const XMLCh* s_external_spectrum_id = xercesc::XMLString::transcode("externalSpectrumID");
			
			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			
			//determine parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
			String parent_parent_tag;
			if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);
			String parent_parent_parent_tag;
			if (open_tags_.size()>3) parent_parent_parent_tag = *(open_tags_.end()-4);

			//do nothing until a new spectrum is reached
			if (tag!="spectrum" && skip_spectrum_) return;

			if (tag=="spectrum")
			{
				//number of peaks
				spec_ = SpectrumType();
				default_array_length_ = attributeAsInt_(attributes, s_default_array_length);
				//spectrum source file
				String source_file_ref;
				if (optionalAttributeAsString_(source_file_ref, attributes, s_source_file_ref))
				{
					spec_.setSourceFile(source_files_[source_file_ref]);
				}
				//native id
				spec_.setNativeID(attributeAsString_(attributes, s_id));
				//maldi spot id
				String maldi_spot_id;
				if(optionalAttributeAsString_(maldi_spot_id,attributes, s_spot_id))
				{
					spec_.setMetaValue("maldi_spot_id",maldi_spot_id);
				}
			}
			else if (tag=="spectrumList")
			{
				//default data processing
				exp_->setDataProcessing(processing_[attributeAsString_(attributes, s_default_data_processing_ref)]);
				
				//Abort if we need meta data only
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
				handleCVParam_(parent_parent_parent_tag, parent_parent_tag, parent_tag, attributeAsString_(attributes, s_accession), attributeAsString_(attributes, s_name), value);
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
					handleCVParam_(parent_parent_parent_tag, parent_parent_tag, parent_tag,ref_param_[ref][i].accession,ref_param_[ref][i].name,ref_param_[ref][i].value);
				}
			}
			else if (tag=="scan")
			{
				//number
				Acquisition tmp;
				tmp.setNumber(attributeAsInt_(attributes, s_number));
				//source file => meta data
				String source_file_ref;
				if(optionalAttributeAsString_(source_file_ref, attributes, s_source_file_ref))
				{
					tmp.setMetaValue("source_file_name",source_files_[source_file_ref].getNameOfFile());
					tmp.setMetaValue("source_file_path",source_files_[source_file_ref].getPathToFile());
				}
//				//external native id => meta data
//				String external_native_id;
//				if(optionalAttributeAsString_(external_native_id, attributes, s_external_native_id))
//				{
//					tmp.setMetaValue("external_native_id",external_native_id);
//				}
//				//external spectrum id => meta data
//				String external_spectrum_id;
//				if(optionalAttributeAsString_(external_spectrum_id, attributes, s_external_spectrum_id))
//				{
//					tmp.setMetaValue("external_spectrum_id",external_spectrum_id);
//				}
				
				spec_.getAcquisitionInfo().push_back(tmp);
			}
			else if (tag=="mzML")
			{
				//check file version against schema version
				String file_version = attributeAsString_(attributes, s_version);
				DoubleReal double_version = 1.0;
				try
				{
					double_version = file_version.toDouble();
				}
				catch(...)
				{
					warning(LOAD, "Could not convert the mzML version string '" + file_version +"' to a double.");
				}
				if (double_version<1.1)
				{
					fatalError(LOAD, "MzML 1.0 is not supported!");
				}
				else if (double_version>version_.toDouble())
				{
					warning(LOAD, "The XML file (" + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
				}
				//handle file accession
				String accession;
				if (optionalAttributeAsString_(accession, attributes, s_accession))
				{
					exp_->setIdentifier(accession);
				}
				//handle file id
				String id;
				if (optionalAttributeAsString_(id, attributes, s_id))
				{
					exp_->setMetaValue("mzml_id",id);
				}
			}
			else if (tag=="contact")
			{
				exp_->getContacts().push_back(ContactPerson());
			}
			else if(tag=="sourceFileRef" && parent_tag=="sourceFileRefList" && parent_parent_tag=="run")
			{
				exp_->getSourceFiles().push_back(source_files_[ attributeAsString_(attributes, s_ref)]);
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
				//sample
				String sample_ref;
				if (optionalAttributeAsString_(sample_ref, attributes, s_sample_ref))
				{
					exp_->setSample(samples_[sample_ref]);
				}
				//instrument
				String instrument_ref = attributeAsString_(attributes, s_default_instrument_configuration_ref);
				exp_->setInstrument(instruments_[instrument_ref]);
				//start time
				String start_time;
				if (optionalAttributeAsString_(start_time, attributes, s_start_time_stamp))
				{
					exp_->setDateTime(asDateTime_(start_time));
				}
			}
			else if (tag=="software")
			{
				current_id_ = attributeAsString_(attributes, s_id);
				software_[current_id_].setVersion(attributeAsString_(attributes, s_version));
			}
			else if (tag=="dataProcessingRef")
			{
			  //"spectrum", "chromatogram" and "binaryDataArray" can have a DataProcessing too
			  //Currently this is not implemented as the mzML schema is messed up
			  //all dataProcessing is handled as global data processing of the whole file!
			  //See: end of dataProcessingList
			}
			else if (tag=="dataProcessing")
			{
				current_id_ = attributeAsString_(attributes, s_id);
			}
			else if (tag=="processingMethod")
			{
				DataProcessing dp;
				dp.setSoftware(software_[attributeAsString_(attributes, s_software_ref)]);
				processing_[current_id_].push_back(dp);
				//The order of processing methods is currently ignored
			}
			else if (tag=="instrumentConfiguration")
			{
				current_id_ = attributeAsString_(attributes, s_id);
			}
			else if (tag=="softwareRef")
			{
				//Set the software of the instrument
				instruments_[current_id_].setSoftware(software_[attributeAsString_(attributes, s_ref)]);
			}
			else if (tag=="source")
			{
				instruments_[current_id_].getIonSources().push_back(IonSource());
				instruments_[current_id_].getIonSources().back().setOrder(attributeAsInt_(attributes, s_order));
			}
			else if (tag=="analyzer")
			{
				instruments_[current_id_].getMassAnalyzers().push_back(MassAnalyzer());
				instruments_[current_id_].getMassAnalyzers().back().setOrder(attributeAsInt_(attributes, s_order));
			}
			else if (tag=="detector")
			{
				instruments_[current_id_].getIonDetectors().push_back(IonDetector());
				instruments_[current_id_].getIonDetectors().back().setOrder(attributeAsInt_(attributes, s_order));
			}
			else if (tag=="precursor")
			{
				//initialize
				spec_.getPrecursor() = Precursor();
				
				//source file => meta data
				String source_file_ref;
				if(optionalAttributeAsString_(source_file_ref, attributes, s_source_file_ref))
				{
					spec_.getPrecursor().setMetaValue("source_file_name",source_files_[source_file_ref].getNameOfFile());
					spec_.getPrecursor().setMetaValue("source_file_path",source_files_[source_file_ref].getPathToFile());
				}
//				//external native id => meta data
//				String external_native_id;
//				if(optionalAttributeAsString_(external_native_id, attributes, s_external_native_id))
//				{
//					spec_.getPrecursor().setMetaValue("external_native_id",external_native_id);
//				}
//				//external spectrum id => meta data
//				String external_spectrum_id;
//				if(optionalAttributeAsString_(external_spectrum_id, attributes, s_external_spectrum_id))
//				{
//					spec_.getPrecursor().setMetaValue("external_spectrum_id",external_spectrum_id);
//				}
			}
			else if (tag=="precursorList")
			{
				//Warn if more than one precursor is present
				if (attributeAsInt_(attributes, s_count)>1)
				{
					warning(LOAD, "OpenMS can only handle one precursor ion! Only the last precursor of each spectrum is loaded!");
				}
			}
			else if (tag=="selectedIon")
			{
				spec_.getPrecursorPeak() = typename SpectrumType::PrecursorPeakType();
			}
			else if (tag=="selectedIonList")
			{
				//Warn if more than one selected ion is present
				if (attributeAsInt_(attributes, s_count)>1)
				{
					warning(LOAD, "OpenMS can only handle one selection ion as precursor! Only the last ion is loaded!");
				}
			}
			else if (tag=="scanWindow")
			{
				spec_.getInstrumentSettings().getScanWindows().push_back(InstrumentSettings::ScanWindow());
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
				if (!skip_spectrum_)
				{
					fillData_();
					exp_->push_back(spec_);
				}
				skip_spectrum_ = false;
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
				instruments_.clear();
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
				if (data_[i].name=="m/z array")
				{
					mz_index = i;
					mz_precision_64 = (data_[i].precision=="64");
				}
				if (data_[i].name=="intensity array")
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
					warning(LOAD, String("The m/z or intensity array of spectrum ") + exp_->size() + " is missing and default_array_length_ is " + default_array_length_ + ".");
				}
				return;
			}
			
			//Warn if the decoded data has a differenct size than the the defaultArrayLength
			UInt mz_size = mz_precision_64 ? data_[mz_index].decoded_64.size() : data_[mz_index].decoded_32.size();
			if (default_array_length_!=mz_size)
			{
				warning(LOAD, String("The base64-decoded m/z array of spectrum ") + exp_->size() + " has the size " + mz_size + ", but it should have the size " + default_array_length_ + " (defaultArrayLength).");
			}
			UInt int_size = int_precision_64 ? data_[int_index].decoded_64.size() : data_[int_index].decoded_32.size();
			if (default_array_length_!=int_size)
			{
				warning(LOAD, String("The base64-decoded intensity array of spectrum ") + exp_->size() + " has the size " + int_size + ", but it should have the size " + default_array_length_ + " (defaultArrayLength).");
			}
			
			//create meta data arrays if necessary
			if (data_.size()>2)
			{
				//create meta data arrays and assign meta data
				spec_.getMetaDataArrays().resize(data_.size()-2);
				UInt meta_array_index = 0;
				for (UInt i=0; i<data_.size(); i++)
				{
					if (data_[i].name!="m/z array" && data_[i].name!="intensity array")
					{
						spec_.getMetaDataArrays()[meta_array_index].setName(data_[i].name);
						spec_.getMetaDataArrays()[meta_array_index].reserve(data_[i].size);
						//copy meta info into MetaInfoDescription
						std::vector<UInt> keys;
						data_[i].meta.getKeys(keys);
						for (UInt k=0;k<keys.size(); ++k)
						{
							spec_.getMetaDataArrays()[meta_array_index].setMetaValue(keys[k],data_[i].meta.getValue(keys[k]));
						}
						//go to next meta data array
						++meta_array_index;
					}
				}
			}
			
			//copy meta data from m/z and intensity binary
			//We don't have this as a separate location => store it in spectrum
			for (UInt i=0; i<data_.size(); i++)
			{
				if (data_[i].name=="m/z array" || data_[i].name=="intensity array")
				{
					std::vector<UInt> keys;
					data_[i].meta.getKeys(keys);
					for (UInt k=0;k<keys.size(); ++k)
					{
						spec_.setMetaValue(keys[k],data_[i].meta.getValue(keys[k]));
					}
				}
			}
			
			//add the peaks and the meta data to the container (if they pass the restrictions)
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
						if (n<data_[i].size && data_[i].name!="m/z array" && data_[i].name!="intensity array")
						{
							DoubleReal value = (data_[i].precision=="64") ? data_[i].decoded_64[n] : data_[i].decoded_32[n];
							spec_.getMetaDataArrays()[meta_array_index].push_back(value);
							++meta_array_index;
						}
					}
				}
			}				
		}
		
		template <typename MapType>
		void MzMLHandler<MapType>::handleCVParam_(const String& /*parent_parent_parent_tag*/, const String& /*parent_parent_tag*/, const String& parent_tag, const String& accession, const String& name, const String& value)
		{
			//Error checks of CV values
			if (cv_.exists(accession))
			{
				const ControlledVocabulary::CVTerm& term = cv_.getTerm(accession);
				//obsolete CV terms
				if (term.obsolete)
				{
					warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
				}
				//check if term name and parsed name match
				String parsed_name = name;
				parsed_name.trim();
				String correct_name = term.name;
				correct_name.trim();
				if (parsed_name!=correct_name)
				{
					warning(LOAD, String("Name of CV term not correct: '") + term.id + " - " + parsed_name + "' should be '" + correct_name + "'");
				}
				if (term.obsolete)
				{
					warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
				}
				//values used in wrong places and wrong value types
				if (value!="")
				{
					if (term.xref_type==ControlledVocabulary::CVTerm::NONE)
					{
						warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must not have a value. The value is '" + value + "'.");
					}
					else
					{
						switch(term.xref_type)
						{
							//string value can be anything
							case ControlledVocabulary::CVTerm::XSD_STRING:
								break;
							//int value => try casting
							case ControlledVocabulary::CVTerm::XSD_INTEGER:
							case ControlledVocabulary::CVTerm::XSD_NEGATIVE_INTEGER:
							case ControlledVocabulary::CVTerm::XSD_POSITIVE_INTEGER:
							case ControlledVocabulary::CVTerm::XSD_NON_NEGATIVE_INTEGER:
							case ControlledVocabulary::CVTerm::XSD_NON_POSITIVE_INTEGER:
								try
								{
									value.toInt();
								}
								catch(Exception::ConversionError&)
								{
									warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have an integer value. The value is '" + value + "'.");
									return;
								}
								break;
							//double value => try casting
							case ControlledVocabulary::CVTerm::XSD_DECIMAL:
								try
								{
									value.toDouble();
								}
								catch(Exception::ConversionError&)
								{
									warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have a floating-point value. The value is '" + value + "'.");
									return;
								}
								break;
							//date string => try conversion
							case ControlledVocabulary::CVTerm::XSD_DATE:
								try
								{
									DateTime tmp;
									tmp.set(value);
								}
								catch(Exception::ParseError&)
								{
									warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must be a valid date. The value is '" + value + "'.");
									return;
								}
								break;								
							default:
								warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' has the unknown value type '" + ControlledVocabulary::CVTerm::getXRefTypeName(term.xref_type) + "'.");
								break;
						}
					}
				}
				//no value, although there should be a numerical value
				else if (term.xref_type!=ControlledVocabulary::CVTerm::NONE && term.xref_type!=ControlledVocabulary::CVTerm::XSD_STRING)
				{
					warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' should have a numerical value. The value is '" + value + "'.");
					return;
	   		}
			}
			
			//------------------------- binaryDataArray ----------------------------
			if (parent_tag=="binaryDataArray")
			{
				if ( in_spectrum_list_)
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
					else if (accession=="MS:1000786") // non-standard binary data array (with name as value)
					{
						data_.back().name = value;
					}
					else if (cv_.isChildOf(accession,"MS:1000513")) //other array names as string
					{
						
						data_.back().name = cv_.getTerm(accession).name;
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
					else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
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
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::FULL);
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
				else if (accession=="MS:1000620") //PDA spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::PDA);
				}
				else if (accession=="MS:1000325") //constant neutral gain spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CNG);
				}
				else if (accession=="MS:1000326") //constant neutral loss spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CNL);
				}
				else if (accession=="MS:1000341") //precursor ion spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::PRECURSOR);
				}
				else if (accession=="MS:1000789") //enhanced multiply charged spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::EMC);
				}
				else if (accession=="MS:1000790") //time-delayed fragmentation spectrum
				{
					spec_.getInstrumentSettings().setScanMode(InstrumentSettings::TDF);
				}
				//representation
				else if (accession=="MS:1000127") //centroid mass spectrum
				{
					spec_.setType(SpectrumSettings::PEAKS);
				}
				else if (accession=="MS:1000128") //profile mass spectrum
				{
					spec_.setType(SpectrumSettings::RAWDATA);
				}
				else if (accession=="MS:1000527" || accession=="MS:1000528" || accession=="MS:1000504" || accession=="MS:1000505" || accession=="MS:1000285" )
				{
					//currently ignored
				}
				else if (accession=="MS:1000511") //ms level
				{
					spec_.setMSLevel(value.toInt());
					
					if (options_.hasMSLevels() && !options_.containsMSLevel(spec_.getMSLevel()))
					{
						skip_spectrum_ = true;
					}
				}
				//scan attributes
				else if (accession=="MS:1000011")//mass resolution
				{
					//No member => meta data
					spec_.setMetaValue("mass resolution",value);
				}
				else if (accession=="MS:1000015")//scan rate
				{
					//No member => meta data
					spec_.setMetaValue("scan rate",value);
				}
				else if (accession=="MS:1000016")//scan time
				{
					spec_.setRT(value.toDouble());
					
					if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRT())))
					{
						skip_spectrum_=true;
					}
				}
				else if (accession=="MS:1000023")//isolation width
				{
					//No member => meta data
					spec_.setMetaValue("isolation width",value);
				}
				else if (accession=="MS:1000512")//filter string
				{
					//No member => meta data
					spec_.setMetaValue("filter string",value);
				}
				else if (accession=="MS:1000616")//preset scan configuration
				{
					//No member => meta data
					spec_.setMetaValue("preset scan configuration",String("true"));
				}
				
				//scan direction
				else if (accession=="MS:1000092")//decreasing m/z scan
				{
					//No member => meta data
					spec_.setMetaValue("scan direction",String("decreasing"));
				}
				else if (accession=="MS:1000093")//increasing m/z scan
				{
					//No member => meta data
					spec_.setMetaValue("scan direction",String("increasing"));
				}
				
				//scan law
				else if (accession=="MS:1000094")//scan law: exponential
				{
					//No member => meta data
					spec_.setMetaValue("scan law",String("exponential"));
				}
				else if (accession=="MS:1000095")//scan law: linear
				{
					//No member => meta data
					spec_.setMetaValue("scan law",String("linear"));
				}
				else if (accession=="MS:1000096")//scan law: quadratic
				{
					//No member => meta data
					spec_.setMetaValue("scan law",String("quadratic"));
				}
				
				//scan polarity
				else if (accession=="MS:1000129")//negative scan
				{
					spec_.getInstrumentSettings().setPolarity(IonSource::NEGATIVE);
				}
				else if (accession=="MS:1000130")//positive scan
				{
					spec_.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			//------------------------- scanWindow ----------------------------
			else if(parent_tag=="scanWindow")
			{
				if (accession=="MS:1000501") //scan m/z lower limit
				{
					spec_.getInstrumentSettings().getScanWindows().back().begin = value.toDouble();
				}
				else if (accession=="MS:1000500") //scan m/z upper limit
				{
					spec_.getInstrumentSettings().getScanWindows().back().end = value.toDouble();
				}
				else if (accession=="MS:1000502") //dwell time
				{
					//Currently ignored
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			//------------------------- referenceableParamGroup ----------------------------
			else if(parent_tag=="referenceableParamGroup")
			{
				SemanticValidator::CVTerm term;
				term.accession = accession;
				term.name = name;
				term.value = value;
				ref_param_[current_id_].push_back(term);
			}
			//------------------------- selectedIon ----------------------------
			else if(parent_tag=="selectedIon")
			{
				if (accession=="MS:1000744") //selected m/z
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
					spec_.getPrecursorPeak().getPossibleChargeStates().push_back(value.toInt());
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			//------------------------- activation ----------------------------
			else if(parent_tag=="activation")
			{
				if (accession=="MS:1000245") //charge stripping
				{
					//No member => meta data
					spec_.getPrecursor().setMetaValue("charge_stripping",String("true"));
				}
				else if (accession=="MS:1000246") //delayed extraction
				{
					//No member => meta data
					spec_.getPrecursor().setMetaValue("delayed extraction",String("true"));
				}
				else if (accession=="MS:1000045") //collision energy
				{
					//No member => meta data
					spec_.getPrecursor().setMetaValue("collision energy",value); 
				}
				else if (accession=="MS:1000412") //buffer gas
				{
					//No member => meta data
					spec_.getPrecursor().setMetaValue("buffer gas",String("true"));
				}
				else if (accession=="MS:1000419") //collision gas
				{
					//No member => meta data
					spec_.getPrecursor().setMetaValue("collision gas",String("true"));
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
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			//------------------------- scanList ----------------------------
			else if(parent_tag=="scanList")
			{
				if (cv_.isChildOf(accession,"MS:1000570")) //method of combination as string
				{
					spec_.getAcquisitionInfo().setMethodOfCombination(cv_.getTerm(accession).name);
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			//------------------------- acquisition ----------------------------
			else if (parent_tag=="scan")
			{
				//scan attributes
				if (accession=="MS:1000011")//mass resolution
				{
					//No member => meta data
					spec_.getAcquisitionInfo().back().setMetaValue("mass resolution",value);
				}
				else if (accession=="MS:1000015")//scan rate
				{
					//No member => meta data
					spec_.getAcquisitionInfo().back().setMetaValue("scan rate",value);
				}
				else if (accession=="MS:1000016")//scan time
				{
					//No member => meta data
					spec_.getAcquisitionInfo().back().setMetaValue("scan time",value);
				}
				else if (accession=="MS:1000023")//isolation width
				{
					//No member => meta data
					spec_.getAcquisitionInfo().back().setMetaValue("isolation width",value);
				}
				else if (accession=="MS:1000511") //ms level
				{
					//No member => meta data
					spec_.getAcquisitionInfo().back().setMetaValue("ms level",value);
				}
				else if (accession=="MS:1000512")//filter string
				{
					//No member => meta data
					spec_.getAcquisitionInfo().back().setMetaValue("filter string",value);
				}
				else if (accession=="MS:1000616")//preset scan configuration
				{
					//No member => meta data
					spec_.getAcquisitionInfo().back().setMetaValue("preset scan configuration",String("true"));
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
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
					exp_->getContacts().back().setAddress(value);
				}
				else if (accession=="MS:1000588") //contact URL
				{
					exp_->getContacts().back().setURL(value);
				}
				else if (accession=="MS:1000589") //contact email
				{
					exp_->getContacts().back().setEmail(value);
				}
				else if (accession=="MS:1000590") //contact organization
				{
					exp_->getContacts().back().setInstitution(value);
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			//------------------------- sourceFile ----------------------------
			else if (parent_tag=="sourceFile")
			{
				if (accession=="MS:1000569") //SHA-1 checksum
				{
					source_files_[current_id_].setChecksum(value, SourceFile::SHA1);
				}
				else if (accession=="MS:1000568") //MD5 checksum
				{
					source_files_[current_id_].setChecksum(value, SourceFile::MD5);
				}
				else if (cv_.isChildOf(accession,"MS:1000560")) //file type as string
				{
					source_files_[current_id_].setFileType(cv_.getTerm(accession).name);
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
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
					//No member => meta data
					samples_[current_id_].setMetaValue("sample batch",String(value));
				}
				else if (accession=="MS:1000047") //emulsion
				{
					samples_[current_id_].setState(Sample::EMULSION);
				}
				else if (accession=="MS:1000048") //gas
				{
					samples_[current_id_].setState(Sample::GAS);
				}
				else if (accession=="MS:1000049") //liquid
				{
					samples_[current_id_].setState(Sample::LIQUID);
				}
				else if (accession=="MS:1000050") //solid
				{
					samples_[current_id_].setState(Sample::SOLID);
				}
				else if (accession=="MS:1000051") //solution
				{
					samples_[current_id_].setState(Sample::SOLUTION);
				}
				else if (accession=="MS:1000052") //suspension
				{
					samples_[current_id_].setState(Sample::SUSPENSION);
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			//------------------------- instrumentConfiguration ----------------------------
			else if (parent_tag=="instrumentConfiguration")
			{
				//instrument model
				if (cv_.isChildOf(accession,"MS:1000031")) //instrument name in value
				{
					instruments_[current_id_].setName(cv_.getTerm(accession).name);
				}
				//instrument attribute
				else if (accession=="MS:1000529") //instrument serial number
				{
					//No member => meta data
					instruments_[current_id_].setMetaValue("instrument serial number",String(value));
				}
				else if (accession=="MS:1000032") //customization
				{
					instruments_[current_id_].setCustomizations(value);
				}
				else if (accession=="MS:1000236") //transmission
				{
					//No member => metadata
					instruments_[current_id_].setMetaValue("transmission",value);
				}
				//ion optics type
				else if (accession=="MS:1000221") //magnetic deflection
				{
					instruments_[current_id_].setIonOptics(Instrument::MAGNETIC_DEFLECTION);
				}
				else if (accession=="MS:1000275") //collision quadrupole
				{
					instruments_[current_id_].setIonOptics(Instrument::COLLISION_QUADRUPOLE);
				}
				else if (accession=="MS:1000281") //selected ion flow tube
				{
					instruments_[current_id_].setIonOptics(Instrument::SELECTED_ION_FLOW_TUBE);
				}
				else if (accession=="MS:1000286") //time lag focusing
				{
					instruments_[current_id_].setIonOptics(Instrument::TIME_LAG_FOCUSING);
				}
				else if (accession=="MS:1000300") //reflectron
				{
					instruments_[current_id_].setIonOptics(Instrument::REFLECTRON);
				}
				else if (accession=="MS:1000307") //einzel lens
				{
					instruments_[current_id_].setIonOptics(Instrument::EINZEL_LENS);
				}
				else if (accession=="MS:1000309") //first stability region
				{
					instruments_[current_id_].setIonOptics(Instrument::FIRST_STABILITY_REGION);
				}
				else if (accession=="MS:1000310") //fringing field
				{
					instruments_[current_id_].setIonOptics(Instrument::FRINGING_FIELD);
				}
				else if (accession=="MS:1000311") //kinetic energy analyzer
				{
					instruments_[current_id_].setIonOptics(Instrument::KINETIC_ENERGY_ANALYZER);
				}
				else if (accession=="MS:1000320") //static field
				{
					instruments_[current_id_].setIonOptics(Instrument::STATIC_FIELD);
				}
				//ion optics attribute
				else if (accession=="MS:1000304") //accelerating voltage
				{
					//No member => metadata
					instruments_[current_id_].setMetaValue("accelerating voltage",value);
				}
				else if (accession=="MS:1000216") //field-free region
				{
					//No member => metadata
					instruments_[current_id_].setMetaValue("field-free region",String("true"));
				}
				else if (accession=="MS:1000308") //electric field strength
				{
					//No member => metadata
					instruments_[current_id_].setMetaValue("electric field strength",value);
				}
				else if (accession=="MS:1000319") //space charge effect
				{
					//No member => metadata
					instruments_[current_id_].setMetaValue("space charge effect",String("true"));
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			else if (parent_tag=="source")
			{
				//inlet type
				if (accession=="MS:1000055") //continuous flow fast atom bombardment
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::CONTINUOUSFLOWFASTATOMBOMBARDMENT);
				}
				else if (accession=="MS:1000056") //direct inlet
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::DIRECT);
				}
				else if (accession=="MS:1000057") //electrospray inlet
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::ELECTROSPRAYINLET);
				}
				else if (accession=="MS:1000058") //flow injection analysis
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::FLOWINJECTIONANALYSIS);
				}
				else if (accession=="MS:1000059") //inductively coupled plasma
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::INDUCTIVELYCOUPLEDPLASMA);
				}
				else if (accession=="MS:1000060") //infusion
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::INFUSION);
				}
				else if (accession=="MS:1000061") //jet separator
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::JETSEPARATOR);
				}
				else if (accession=="MS:1000062") //membrane separator
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::MEMBRANESEPARATOR);
				}
				else if (accession=="MS:1000063") //moving belt
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::MOVINGBELT);
				}
				else if (accession=="MS:1000064") //moving wire
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::MOVINGWIRE);
				}
				else if (accession=="MS:1000065") //open split
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::OPENSPLIT);
				}
				else if (accession=="MS:1000066") //particle beam
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::PARTICLEBEAM);
				}
				else if (accession=="MS:1000067") //reservoir
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::RESERVOIR);
				}
				else if (accession=="MS:1000068") //septum
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::SEPTUM);
				}
				else if (accession=="MS:1000069") //thermospray inlet
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::THERMOSPRAYINLET);
				}
				else if (accession=="MS:1000248") //direct insertion probe
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::BATCH);
				}
				else if (accession=="MS:1000249") //direct liquid introduction
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::CHROMATOGRAPHY);
				}
				else if (accession=="MS:1000396") //membrane inlet
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::MEMBRANE);
				}
				else if (accession=="MS:1000485") //nanospray inlet
				{
					instruments_[current_id_].getIonSources().back().setInletType(IonSource::NANOSPRAY);
				}
				//ionization type
				else if (accession=="MS:1000071") //chemical ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::CI);
				}
				else if (accession=="MS:1000073") //electrospray ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::ESI);
				}
				else if (accession=="MS:1000074") //fast atom bombardment ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FAB);
				}
				else if (accession=="MS:1000227") //multiphoton ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::MPI);
				}
				else if (accession=="MS:1000240") //atmospheric pressure ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::API);
				}
				else if (accession=="MS:1000247") //desorption ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::DI);
				}
				else if (accession=="MS:1000255") //flowing afterglow
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FA);
				}
				else if (accession=="MS:1000258") //field ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FII);
				}
				else if (accession=="MS:1000259") //glow discharge ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::GD_MS);
				}
				else if (accession=="MS:1000271") //Negative ion chemical ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::NICI);
				}
				else if (accession=="MS:1000272") //neutralization reionization mass spectrometry
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::NRMS);
				}
				else if (accession=="MS:1000273") //photoionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::PI);
				}
				else if (accession=="MS:1000274") //pyrolysis mass spectrometry
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::PYMS);
				}
				else if (accession=="MS:1000276") //resonance enhanced multiphoton ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::REMPI);
				}
				else if (accession=="MS:1000380") //adiabatic ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::AI);
				}
				else if (accession=="MS:1000381") //associative ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::ASI);
				}
				else if (accession=="MS:1000383") //autodetachment
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::AD);
				}
				else if (accession=="MS:1000384") //autoionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::AUI);
				}
				else if (accession=="MS:1000385") //charge exchange ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::CEI);
				}
				else if (accession=="MS:1000386") //chemi-ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::CHEMI);
				}
				else if (accession=="MS:1000388") //dissociative ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::DISSI);
				}
				else if (accession=="MS:1000389") //electron ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::EI);
				}
				else if (accession=="MS:1000395") //liquid secondary ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::LSI);
				}
				else if (accession=="MS:1000399") //penning ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::PEI);
				}
				else if (accession=="MS:1000400") //plasma desorption ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::PD);
				}
				else if (accession=="MS:1000402") //secondary ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SI);
				}
				else if (accession=="MS:1000403") //soft ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SOI);
				}
				else if (accession=="MS:1000404") //spark ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SPI);
				}
				else if (accession=="MS:1000406") //surface ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SUI);
				}
				else if (accession=="MS:1000407") //thermal ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::TI);
				}
				else if (accession=="MS:1000408") //vertical ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::VI);
				}
				else if (accession=="MS:1000446") //fast ion bombardment
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FIB);
				}
				else if (accession=="MS:1000070") //atmospheric pressure chemical ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::APCI);
				}
				else if (accession=="MS:1000239") //atmospheric pressure matrix-assisted laser desorption ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::AP_MALDI);
				}
				else if (accession=="MS:1000382") //atmospheric pressure photoionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::APPI);
				}
				else if (accession=="MS:1000075") //matrix-assisted laser desorption ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::MALDI);
				}
				else if (accession=="MS:1000257") //field desorption
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::FD);
				}
				else if (accession=="MS:1000387") //desorption/ionization on silicon
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SILI);
				}
				else if (accession=="MS:1000393") //laser desorption ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::LD);
				}
				else if (accession=="MS:1000405") //surface-assisted laser desorption ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SALDI);
				}
				else if (accession=="MS:1000397") //microelectrospray
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::MESI);
				}
				else if (accession=="MS:1000398") //nanoelectrospray
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::NESI);
				}
				else if (accession=="MS:1000278") //surface enhanced laser desorption ionization
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SELDI);
				}
				else if (accession=="MS:1000279") //surface enhanced neat desorption
				{
					instruments_[current_id_].getIonSources().back().setIonizationMethod(IonSource::SEND);
				}
				
				//source attribute
				else if (accession=="MS:1000392") //ionization efficiency
				{
					//No member => meta data
					instruments_[current_id_].getIonSources().back().setMetaValue("ionization efficiency",value);
				}
				else if (accession=="MS:1000486") //source potential
				{
					//No member => meta data
					instruments_[current_id_].getIonSources().back().setMetaValue("source potential",value);
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			else if (parent_tag=="analyzer")
			{
				//mass analyzer type
				if (accession=="MS:1000079") //fourier transform ion cyclotron resonance mass spectrometer
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::FOURIERTRANSFORM);
				}
				else if (accession=="MS:1000080") //magnetic sector
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::SECTOR);
				}
				else if (accession=="MS:1000081") //quadrupole
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::QUADRUPOLE);
				}
				else if (accession=="MS:1000084") //time-of-flight
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::TOF);
				}
				else if (accession=="MS:1000254") //electrostatic energy analyzer
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::ESA);
				}
				else if (accession=="MS:1000264") //ion trap
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::IT);
				}
				else if (accession=="MS:1000284") //stored waveform inverse fourier transform
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::SWIFT);
				}
				else if (accession=="MS:1000288") //cyclotron
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::CYCLOTRON);
				}
				else if (accession=="MS:1000484") //orbitrap
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::ORBITRAP);
				}
				else if (accession=="MS:1000078") //axial ejection linear ion trap
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::AXIALEJECTIONLINEARIONTRAP);
				}
				else if (accession=="MS:1000082") //quadrupole ion trap
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::PAULIONTRAP);
				}
				else if (accession=="MS:1000083") //radial ejection linear ion trap
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::RADIALEJECTIONLINEARIONTRAP);
				}
				else if (accession=="MS:1000291") //linear ion trap
				{
					instruments_[current_id_].getMassAnalyzers().back().setType(MassAnalyzer::LIT);
				}
				//mass analyzer attribute
				else if (accession=="MS:1000014") //accuracy
				{
					instruments_[current_id_].getMassAnalyzers().back().setAccuracy(value.toDouble());
				}
				else if (accession=="MS:1000022") //TOF Total Path Length
				{
					instruments_[current_id_].getMassAnalyzers().back().setTOFTotalPathLength(value.toDouble());
				}
				else if (accession=="MS:1000024") //final MS exponent
				{
					instruments_[current_id_].getMassAnalyzers().back().setFinalMSExponent(value.toInt());
				}
				else if (accession=="MS:1000025") //magnetic field strength
				{
					instruments_[current_id_].getMassAnalyzers().back().setMagneticFieldStrength(value.toDouble());
				}
				else if (accession=="MS:1000105") //reflectron off
				{
					instruments_[current_id_].getMassAnalyzers().back().setReflectronState(MassAnalyzer::OFF);
				}
				else if (accession=="MS:1000106") //reflectron on
				{
					instruments_[current_id_].getMassAnalyzers().back().setReflectronState(MassAnalyzer::ON);
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			else if (parent_tag=="detector")
			{
				//detector type
				if (accession=="MS:1000107") //channeltron
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::CHANNELTRON);
				}
				else if (accession=="MS:1000110") //daly detector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::DALYDETECTOR);
				}
				else if (accession=="MS:1000112") //faraday cup
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::FARADAYCUP);
				}
				else if (accession=="MS:1000114") //microchannel plate detector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::MICROCHANNELPLATEDETECTOR);
				}
				else if (accession=="MS:1000115") //multi-collector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::MULTICOLLECTOR);
				}
				else if (accession=="MS:1000116") //photomultiplier
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::PHOTOMULTIPLIER);
				}
				else if (accession=="MS:1000253") //electron multiplier
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::ELECTRONMULTIPLIER);
				}
				else if (accession=="MS:1000345") //array detector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::ARRAYDETECTOR);
				}
				else if (accession=="MS:1000346") //conversion dynode
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::CONVERSIONDYNODE);
				}
				else if (accession=="MS:1000347") //dynode
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::DYNODE);
				}
				else if (accession=="MS:1000348") //focal plane collector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::FOCALPLANECOLLECTOR);
				}
				else if (accession=="MS:1000349") //ion-to-photon detector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::IONTOPHOTONDETECTOR);
				}
				else if (accession=="MS:1000350") //point collector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::POINTCOLLECTOR);
				}
				else if (accession=="MS:1000351") //postacceleration detector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::POSTACCELERATIONDETECTOR);
				}
				else if (accession=="MS:1000621") //photodiode array detector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::PHOTODIODEARRAYDETECTOR);
				}
				else if (accession=="MS:1000624") //inductive detector
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::INDUCTIVEDETECTOR);
				}
				else if (accession=="MS:1000108") //conversion dynode electron multiplier
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::CONVERSIONDYNODEELECTRONMULTIPLIER);
				}
				else if (accession=="MS:1000109") //conversion dynode photomultiplier
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::CONVERSIONDYNODEPHOTOMULTIPLIER);
				}
				else if (accession=="MS:1000111") //electron multiplier tube
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::ELECTRONMULTIPLIERTUBE);
				}
				else if (accession=="MS:1000113") //focal plane array
				{
					instruments_[current_id_].getIonDetectors().back().setType(IonDetector::FOCALPLANEARRAY);
				}
				//detector attribute
				else if (accession=="MS:1000028") //detector resolution
				{
					instruments_[current_id_].getIonDetectors().back().setResolution(value.toDouble());
				}
				else if (accession=="MS:1000029") //sampling frequency
				{
					instruments_[current_id_].getIonDetectors().back().setADCSamplingFrequency(value.toDouble());
				}
				//dectector acquisition mode
				else if (accession=="MS:1000117") //analog-digital converter
				{
					instruments_[current_id_].getIonDetectors().back().setAcquisitionMode(IonDetector::ADC);
				}
				else if (accession=="MS:1000118") //pulse counting
				{
					instruments_[current_id_].getIonDetectors().back().setAcquisitionMode(IonDetector::PULSECOUNTING);
				}
				else if (accession=="MS:1000119") //time-digital converter
				{
					instruments_[current_id_].getIonDetectors().back().setAcquisitionMode(IonDetector::TDC);
				}
				else if (accession=="MS:1000120") //transient recorder
				{
					instruments_[current_id_].getIonDetectors().back().setAcquisitionMode(IonDetector::TRANSIENTRECORDER);
				} 
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			else if (parent_tag=="processingMethod")
			{
				//data processing parameter
				if (accession=="MS:1000629") //low intensity threshold
				{
					processing_[current_id_].back().setMetaValue("low_intensity_threshold",value.toDouble());
				}
				else if (accession=="MS:1000631") //high intensity threshold
				{
					processing_[current_id_].back().setMetaValue("high_intensity_threshold",value.toDouble());
				}
				else if (accession=="MS:1000787") //inclusive low intensity threshold
				{
					processing_[current_id_].back().setMetaValue("inclusive_low_intensity_threshold",value.toDouble());
				}
				else if (accession=="MS:1000788") //inclusive high intensity threshold
				{
					processing_[current_id_].back().setMetaValue("inclusive_high_intensity_threshold",value.toDouble());
				}
				else if (accession=="MS:1000747") //completion time
				{
					processing_[current_id_].back().setCompletionTime(asDateTime_(value));
				}
				//file format conversion
				else if (accession=="MS:1000544") //Conversion to mzML
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::CONVERSION_MZML);
				}
				else if (accession=="MS:1000545") //Conversion to mzXML
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::CONVERSION_MZXML);
				}
				else if (accession=="MS:1000546") //Conversion to mzData
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::CONVERSION_MZDATA);
				}
				else if (accession=="MS:1000741") //Conversion to DTA
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::CONVERSION_DTA);
				}
				//data processing action
				else if (accession=="MS:1000033") //deisotoping
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::DEISOTOPING);
				}
				else if (accession=="MS:1000034") //charge deconvolution
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::CHARGE_DECONVOLUTION);
				}
				else if (accession=="MS:1000035") //peak picking
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::PEAK_PICKING);
				}
				else if (accession=="MS:1000592" || cv_.isChildOf(accession,"MS:1000592")) //smoothing (or child terms, we make no difference
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::SMOOTHING);
				}
				else if (accession=="MS:1000778" || cv_.isChildOf(accession,"MS:1000778")) //charge state calculation (or child terms, we make no difference
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::CHARGE_CALCULATION);
				}
				else if (accession=="MS:1000780" || cv_.isChildOf(accession,"MS:1000780")) //precursor recalculation (or child terms, we make no difference
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::PRECURSOR_RECALCULATION);
				}
				else if (accession=="MS:1000593") //baseline reduction
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::BASELINE_REDUCTION);
				}
				else if (accession=="MS:1000594") //low intensity data point removal
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::LOW_INTENSITY_REMOVAL);
				}
				else if (accession=="MS:1000745") //retention time alignment
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::ALIGNMENT);
				}
				else if (accession=="MS:1000746") //high intensity data point removal
				{
					processing_[current_id_].back().getProcessingActions().insert(DataProcessing::HIGH_INTENSITY_REMOVAL);
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			else if (parent_tag=="fileContent")
			{
				if (accession=="MS:1000768") //Thermo nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::THERMO);
				}
				else if (accession=="MS:1000769") //Waters nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::WATERS);
				}
				else if (accession=="MS:1000770") //WIFF nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::WIFF);
				}
				else if (accession=="MS:1000771") //Bruker/Agilent YEP nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::BRUKER_AGILENT);
				}
				else if (accession=="MS:1000772") //Bruker BAF nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::BRUKER_BAF);
				}
				else if (accession=="MS:1000773") //Bruker FID nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::BRUKER_FID);
				}
				else if (accession=="MS:1000774") //multiple peak list nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::MULTIPLE_PEAK_LISTS);
				}
				else if (accession=="MS:1000775") //single peak list nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::SINGLE_PEAK_LIST);
				}
				else if (accession=="MS:1000776") //scan number only nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::SCAN_NUMBER);
				}
				else if (accession=="MS:1000777") //spectrum identifier nativeID format
				{
					exp_->setNativeIDType(ExperimentalSettings::SPECTRUM_IDENTIFIER);
				}
				else if (cv_.isChildOf(accession,"MS:1000524")) //data file content
				{
					//ignored
				}
				else if (cv_.isChildOf(accession,"MS:1000525")) //spectrum representation
				{
					//ignored
				}
				else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			}
			else if (parent_tag=="software")
			{
				//Using an enum for software names is not really practical in my (Marc) opinion
				// => we simply store the name as string
				software_[current_id_].setName(name);
			}
			else if (parent_tag=="chromatogram" || parent_tag=="target")
			{
				//allowed but, not needed
			}
			else warning(LOAD, String("Unhandled cvParam '") + accession + " in tag '" + parent_tag + "'.");
			
		}

		template <typename MapType>
		void MzMLHandler<MapType>::handleUserParam_(const String& parent_tag, const String& name, const String& type, const String& value)
		{
			//create a DataValue that contains the data in the right type
			DataValue data_value;
			//float type
			if (type=="xsd:double" || type=="xsd:float")
			{
				data_value = DataValue(value.toDouble());
			}
			//integer type
			else if (type=="xsd:byte" || type=="xsd:decimal" || type=="xsd:int" || type=="xsd:integer" || type=="xsd:long" || type=="xsd:negativeInteger" || type=="xsd:nonNegativeInteger" || type=="xsd:nonPositiveInteger" || type=="xsd:positiveInteger" || type=="xsd:short" || type=="xsd:unsignedByte" || type=="xsd:unsignedInt" || type=="xsd:unsignedLong" || type=="xsd:unsignedShort")
			{
				data_value = DataValue(value.toInt());
			}
			//everything else is treated as a string
			else
			{
				data_value = DataValue(value);
			}
			
			//find the right MetaInfoInterface
			if (parent_tag=="run")
			{
				exp_->setMetaValue(name,data_value);
			}
			else if (parent_tag=="instrumentConfiguration")
			{
				instruments_[current_id_].setMetaValue(name,data_value);
			}
			else if (parent_tag=="source")
			{
				instruments_[current_id_].getIonSources().back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="analyzer")
			{
				instruments_[current_id_].getMassAnalyzers().back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="detector")
			{
				instruments_[current_id_].getIonDetectors().back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="sample")
			{
				samples_[current_id_].setMetaValue(name,data_value);
			}
			else if (parent_tag=="software")
			{
				software_[current_id_].setMetaValue(name,data_value);
			}
			else if (parent_tag=="contact")
			{
				exp_->getContacts().back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="sourceFile")
			{
				source_files_[current_id_].setMetaValue(name,data_value);
			}
			else if (parent_tag=="binaryDataArray")
			{
				data_.back().meta.setValue(name,data_value);
			}
			else if (parent_tag=="spectrum")
			{
				spec_.setMetaValue(name,data_value);
			}
			else if (parent_tag=="scanList")
			{
				spec_.getAcquisitionInfo().setMetaValue(name,data_value);
			}
			else if (parent_tag=="scan")
			{
				spec_.getAcquisitionInfo().back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="isolationWindow")
			{
				//We don't have this as a separate location => store it in the precursor
				spec_.getPrecursor().setMetaValue(name,data_value);
			}
			else if (parent_tag=="selectedIon")
			{
				//We don't have this as a separate location => store it in the precursor
				spec_.getPrecursor().setMetaValue(name,data_value);
			}
			else if (parent_tag=="activation")
			{
				//We don't have this as a separate location => store it in the precursor
				spec_.getPrecursor().setMetaValue(name,data_value);
			}
			else if (parent_tag=="processingMethod")
			{
				processing_[current_id_].back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="fileContent")
			{
				//currently ignored
			}
			else warning(LOAD, String("Unhandled userParam '") + name + " in tag '" + parent_tag + "'.");
			
		}
	
		template <typename MapType>
		void MzMLHandler<MapType>::writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent) const
		{
			std::vector<String> keys;
			meta.getKeys(keys);
			
			for (UInt i = 0; i!=keys.size();++i)
			{
				os << String(indent,'\t') << "<userParam name=\"" << keys[i] << "\" type=\"";
				
				DataValue d = meta.getMetaValue(keys[i]);
				//determine type
				if (d.valueType()==DataValue::INT_VALUE)
				{
					os << "xsd:integer";
				}
				else if (d.valueType()==DataValue::DOUBLE_VALUE)
				{
					os << "xsd:double";
				}
				else //string or lists are converted to string
				{
					os << "xsd:string";
				}
				os << "\" value=\"" << (String)(d) << "\"/>" << std::endl;
			}
		}
		
		template <typename MapType>
		ControlledVocabulary::CVTerm MzMLHandler<MapType>::getChildWithName_(const String& parent_accession, const String& name) const
		{
			std::set<String> terms;
			cv_.getAllChildTerms(terms, parent_accession);
			for (std::set<String>::const_iterator it=terms.begin(); it!=terms.end(); ++it)
			{
				if (cv_.getTerm(*it).name==name)
				{
					return cv_.getTerm(*it);
				}
			}
			return ControlledVocabulary::CVTerm();
		}
		
		template <typename MapType>
		void MzMLHandler<MapType>::writeSoftware_(std::ostream& os, const String& id, const Software& software)
		{
			os  << "		<software id=\"" << id << "\" version=\"" << software.getVersion() << "\" >\n";
			ControlledVocabulary::CVTerm so_term = getChildWithName_("MS:1000531",software.getName());
			if (so_term.id!="")
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"" << so_term.id << "\" name=\"" << so_term.name << "\" />\n";
			}
			else //FORCED
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000752\" name=\"TOPP software\" />\n";
			}
			writeUserParam_(os, software, 3);
			os  << "		</software>\n";
		}

		template <typename MapType>
		void MzMLHandler<MapType>::writeSourceFile_(std::ostream& os, const String& id, const SourceFile& source_file)
		{
			os	<< "			<sourceFile id=\"" << id << "\" name=\"" << source_file.getNameOfFile() << "\" location=\"" << source_file.getPathToFile() << "\">\n";
			if (source_file.getChecksumType()==SourceFile::SHA1)
			{
				os  << "				<cvParam cvRef=\"MS\" accession=\"MS:1000569\" name=\"SHA-1\" value=\"" << source_file.getChecksum() << "\"/>\n";
			}
			else if (source_file.getChecksumType()==SourceFile::MD5)
			{
				os  << "				<cvParam cvRef=\"MS\" accession=\"MS:1000568\" name=\"MD5\" value=\"" << source_file.getChecksum() << "\"/>\n";
			}
			else //FORCED
			{
				os  << "				<cvParam cvRef=\"MS\" accession=\"MS:1000569\" name=\"SHA-1\" value=\"\"/>\n";
			}
			ControlledVocabulary::CVTerm sf_term = getChildWithName_("MS:1000560",source_file.getFileType());
			if (sf_term.id!="")
			{
				os  << "				<cvParam cvRef=\"MS\" accession=\"" << sf_term.id <<"\" name=\"" << sf_term.name << "\"/>\n";
			}
			else //FORCED
			{
				os  << "				<cvParam cvRef=\"MS\" accession=\"MS:1000564\" name=\"mzData file\"/>\n";
			}
			writeUserParam_(os, source_file, 4);
			os	<< "			</sourceFile>\n";
		}

		template <typename MapType>
		void MzMLHandler<MapType>::writeTo(std::ostream& os)
		{
			const MapType& exp = *(cexp_);
			logger_.startProgress(0,exp.size(),"storing mzML file");
			
			os	<< "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
					<< "<mzML xmlns=\"http://psi.hupo.org/schema_revision/mzML_1.0.0\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/schema_revision/mzML_1.0.0 mzML1.0.0.xsd\" accession=\"" << exp.getIdentifier() << "\" version=\"1.1\">\n";
			//--------------------------------------------------------------------------------------------
			// CV list
			//--------------------------------------------------------------------------------------------
			os	<< "	<cvList count=\"2\">\n"
					<< "		<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\"/>\n"
					<< "		<cv id=\"UO\" fullName=\"Unit Ontology\" URI=\"http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo\"/>\n"
					<< "	</cvList>\n";
			//--------------------------------------------------------------------------------------------
			// file content
			//--------------------------------------------------------------------------------------------
			os	<< "	<fileDescription>\n";
			os	<< "		<fileContent>\n";
			Map<InstrumentSettings::ScanMode, UInt> file_content;
			for (UInt i=0; i<exp.size(); ++i)
			{
				file_content[exp[i].getInstrumentSettings().getScanMode()]++;
			}
			if (file_content.has(InstrumentSettings::FULL))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::ZOOM))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::SIM))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000582\" name=\"SIM spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::SRM))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000583\" name=\"SRM spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::CRM))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000581\" name=\"CRM spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::PRECURSOR))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000341\" name=\"precursor ion spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::CNG))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000325\" name=\"constant neutral gain spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::CNL))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000326\" name=\"constant neutral loss spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::PDA))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000620\" name=\"PDA spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::EMC))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000789\" name=\"enhanced multiply charged spectrum\" />\n";
			}
			else if (file_content.has(InstrumentSettings::TDF))
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000789\" name=\"time-delayed fragmentation spectrum\" />\n";
			}
			else //FORCED
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" />\n";
			}
			//native ID format
			if (exp.getNativeIDType()==ExperimentalSettings::THERMO)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000768\" name=\"Thermo nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::WATERS)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000769\" name=\"Waters nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::WIFF)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000770\" name=\"WIFF nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::BRUKER_AGILENT)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000771\" name=\"Bruker/Agilent YEP nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::BRUKER_BAF)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000772\" name=\"Bruker BAF nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::BRUKER_FID)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000773\" name=\"Bruker FID nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::MULTIPLE_PEAK_LISTS)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000774\" name=\"multiple peak list nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::SINGLE_PEAK_LIST)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000775\" name=\"single peak list nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::SCAN_NUMBER)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000776\" name=\"scan number only nativeID format\" />\n";
			}
			else if (exp.getNativeIDType()==ExperimentalSettings::SPECTRUM_IDENTIFIER)
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000777\" name=\"spectrum identifier nativeID format\" />\n";
			}
			else //FORCED
			{
				os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000774\" name=\"multiple peak list nativeID format\" />\n";
			}

			os	<< "		</fileContent>\n";
			//--------------------------------------------------------------------------------------------
			// source file list
			//--------------------------------------------------------------------------------------------
			if (exp.getSourceFiles().size()!=0)
			{
				//find out how many spectra source files need to be written
				UInt sf_sp_count=0;
				for (UInt i=0; i<exp.size(); ++i)
				{
					if (exp[i].getSourceFile()!=SourceFile()) ++sf_sp_count;
				}
				os	<< "		<sourceFileList count=\"" << exp.getSourceFiles().size() + sf_sp_count << "\">\n";
				//write source files of run
				for (UInt i=0; i<exp.getSourceFiles().size(); ++i)
				{
					writeSourceFile_(os, String("sf_ru_")+i, exp.getSourceFiles()[i]);
				}
				//write source files of spectra
				for (UInt i=0; i<exp.size(); ++i)
				{
					if (exp[i].getSourceFile()!=SourceFile())
					{
						writeSourceFile_(os, String("sf_sp_")+i, exp[i].getSourceFile());
					}
				}
				os	<< "		</sourceFileList>\n";
			}
			//--------------------------------------------------------------------------------------------
			// contacts
			//--------------------------------------------------------------------------------------------
			for (UInt i=0; i<exp.getContacts().size(); ++i)
			{
				const ContactPerson& cp = exp.getContacts()[i];
				os  << "		<contact>\n";
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000586\" name=\"contact name\" value=\"" << cp.getLastName() << ", " << cp.getFirstName() << "\"/>\n";
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000590\" name=\"contact organization\" value=\"" << cp.getInstitution() << "\"/>\n";

				if (cp.getAddress()!="")
				{
					os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000587\" name=\"contact address\" value=\"" << cp.getAddress() << "\"/>\n";
				}
				if (cp.getURL()!="")
				{
					os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000588\" name=\"contact URL\" value=\"" << cp.getURL() << "\"/>\n";
				}
				if (cp.getEmail()!="")
				{
					os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000589\" name=\"contact email\" value=\"" << cp.getEmail() << "\"/>\n";
				}
				if (cp.getContactInfo()!="")
				{
					os  << "			<userParam name=\"contact_info\" type=\"xsd:string\" value=\"" << cp.getContactInfo() << "\"/>\n";
				}
				writeUserParam_(os, cp, 3);
				os << "		</contact>\n"; 
			}
			os	<< "	</fileDescription>\n";
			//--------------------------------------------------------------------------------------------
			// sample
			//--------------------------------------------------------------------------------------------
			const Sample& sa = exp.getSample();
			os  << "	<sampleList count=\"1\">\n";
			os  << "		<sample id=\"sa_0\" name=\"" << sa.getName() << "\">\n";
			if (sa.getNumber()!="")
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000001\" name=\"sample number\" value=\"" << sa.getNumber() << "\"/>\n";
			}
			os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000004\" name=\"sample mass\" value=\"" << sa.getMass() << "\"/>\n";
			os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000005\" name=\"sample volume\" value=\"" << sa.getVolume() << "\"/>\n";
			os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000006\" name=\"sample concentration\" value=\"" << sa.getConcentration() << "\"/>\n";
			if (sa.getComment()!="")
			{
				os  << "			<userParam name=\"comment\" type=\"xsd:string\" value=\"" << sa.getComment() << "\"/>\n";
			}
			if (sa.getState()==Sample::EMULSION)
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000047\" name=\"emulsion\" />\n";
			}
			else if (sa.getState()==Sample::GAS)
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000048\" name=\"gas\" />\n";
			}
			else if (sa.getState()==Sample::LIQUID)
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000049\" name=\"liquid\" />\n";
			}
			else if (sa.getState()==Sample::SOLID)
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000050\" name=\"solid\" />\n";
			}
			else if (sa.getState()==Sample::SOLUTION)
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000051\" name=\"solution\" />\n";
			}
			else if (sa.getState()==Sample::SUSPENSION)
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000052\" name=\"suspension\" />\n";
			}

			writeUserParam_(os, sa, 3);
			os  << "		</sample>\n";
			os  << "	</sampleList>\n";

			//--------------------------------------------------------------------------------------------
			// software
			//--------------------------------------------------------------------------------------------
			os  << "	<softwareList count=\"" << std::max((UInt)2,(UInt)exp.getDataProcessing().size()+1) << "\">\n";			
			//write instrument software
			writeSoftware_(os, "so_in_0", exp.getInstrument().getSoftware());
			//write data processing
			for (UInt i=0; i<exp.getDataProcessing().size(); ++i)
			{
				writeSoftware_(os, String("so_dp_") + i, exp.getDataProcessing()[i].getSoftware());
			}
			//FORCED - for DataProcessing
			if (exp.getDataProcessing().size()==0)
			{
				writeSoftware_(os, "so_dp_0", Software());
			}
			os  << "	</softwareList>\n";			

			//--------------------------------------------------------------------------------------------
			// instrument configuration (enclosing ion source, mass analyzer and detector)
			//--------------------------------------------------------------------------------------------
			const Instrument& in = exp.getInstrument();
			os  << "	<instrumentConfigurationList count=\"1\">\n";
			os  << "		<instrumentConfiguration id=\"ic_0\">\n";
			ControlledVocabulary::CVTerm in_term = getChildWithName_("MS:1000031",in.getName());
			if (in_term.id!="")
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"" << in_term.id <<"\" name=\"" << in_term.name << "\"/>\n";
			}
			else //FORCED
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000447\" name=\"LTQ\"/>\n";
			}
			
			if (in.getCustomizations()!="")
			{
				os  << "			<cvParam cvRef=\"MS\" accession=\"MS:1000032\" name=\"customization\" value=\"" << in.getCustomizations() << "\"/>\n";
			}

			//ion optics
			if (in.getIonOptics()==Instrument::MAGNETIC_DEFLECTION)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000221\" name=\"magnetic deflection\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::DELAYED_EXTRACTION)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000246\" name=\"delayed extraction\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::COLLISION_QUADRUPOLE)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000275\" name=\"collision quadrupole\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::SELECTED_ION_FLOW_TUBE)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000281\" name=\"selected ion flow tube\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::TIME_LAG_FOCUSING)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000286\" name=\"time lag focusing\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::REFLECTRON)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000300\" name=\"reflectron\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::EINZEL_LENS)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000307\" name=\"einzel lens\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::ELECTRIC_FIELD_STRENGTH)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000308\" name=\"electric field strength\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::FIRST_STABILITY_REGION)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000309\" name=\"first stability region\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::FRINGING_FIELD)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000310\" name=\"fringing field\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::KINETIC_ENERGY_ANALYZER)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000311\" name=\"kinetic energy analyzer\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::SPACE_CHARGE_EFFECT)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000319\" name=\"space charge effect\"/>\n";
			}
			else if (in.getIonOptics()==Instrument::STATIC_FIELD)
			{
				os << "			<cvParam cvRef=\"MS\" accession=\"MS:1000320\" name=\"static field\"/>\n";
			}

			writeUserParam_(os, in, 3);
			UInt component_count = in.getIonSources().size() + in.getMassAnalyzers().size() + in.getIonDetectors().size();
			if (component_count!=0)
			{
				os  << "			<componentList count=\"" << std::max((UInt)3,component_count) << "\">\n";
				//--------------------------------------------------------------------------------------------
				// ion source
				//--------------------------------------------------------------------------------------------
				for (UInt i=0; i<in.getIonSources().size(); ++i)
				{
					const IonSource& so = in.getIonSources()[i];
					os  << "				<source order=\"" << so.getOrder() << "\">\n";
					
					if(so.getInletType()==IonSource::CONTINUOUSFLOWFASTATOMBOMBARDMENT)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000055\" name=\"continuous flow fast atom bombardment\"/>\n";
					}
					else if(so.getInletType()==IonSource::DIRECT)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000056\" name=\"direct inlet\"/>\n";
					}
					else if(so.getInletType()==IonSource::ELECTROSPRAYINLET)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000057\" name=\"electrospray inlet\"/>\n";
					}
					else if(so.getInletType()==IonSource::FLOWINJECTIONANALYSIS)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000058\" name=\"flow injection analysis\"/>\n";
					}
					else if(so.getInletType()==IonSource::INDUCTIVELYCOUPLEDPLASMA)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000059\" name=\"inductively coupled plasma\"/>\n";
					}
					else if(so.getInletType()==IonSource::INFUSION)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000060\" name=\"infusion\"/>\n";
					}
					else if(so.getInletType()==IonSource::JETSEPARATOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000061\" name=\"jet separator\"/>\n";
					}
					else if(so.getInletType()==IonSource::MEMBRANESEPARATOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000062\" name=\"membrane separator\"/>\n";
					}
					else if(so.getInletType()==IonSource::MOVINGBELT)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000063\" name=\"moving belt\"/>\n";
					}
					else if(so.getInletType()==IonSource::MOVINGWIRE)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000064\" name=\"moving wire\"/>\n";
					}
					else if(so.getInletType()==IonSource::OPENSPLIT)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000065\" name=\"open split\"/>\n";
					}
					else if(so.getInletType()==IonSource::PARTICLEBEAM)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000066\" name=\"particle beam\"/>\n";
					}
					else if(so.getInletType()==IonSource::RESERVOIR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000067\" name=\"reservoir\"/>\n";
					}
					else if(so.getInletType()==IonSource::SEPTUM)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000068\" name=\"septum\"/>\n";
					}
					else if(so.getInletType()==IonSource::THERMOSPRAYINLET)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000069\" name=\"thermospray inlet\"/>\n";
					}
					else if(so.getInletType()==IonSource::BATCH)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000248\" name=\"direct insertion probe\"/>\n";
					}
					else if(so.getInletType()==IonSource::CHROMATOGRAPHY)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000249\" name=\"direct liquid introduction\"/>\n";
					}
					else if(so.getInletType()==IonSource::MEMBRANE)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000396\" name=\"membrane inlet\"/>\n";
					}
					else if(so.getInletType()==IonSource::NANOSPRAY)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000485\" name=\"nanospray inlet\"/>\n";	
					}
	
					if(so.getIonizationMethod()==IonSource::APCI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000070\" name=\"atmospheric pressure chemical ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::CI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000071\" name=\"chemical ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::ESI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000073\" name=\"electrospray ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::FAB)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000074\" name=\"fast atom bombardment ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::MALDI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000075\" name=\"matrix-assisted laser desorption ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::MPI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000227\" name=\"multiphoton ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::AP_MALDI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000239\" name=\"atmospheric pressure matrix-assisted laser desorption ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::API)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000240\" name=\"atmospheric pressure ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::DI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000247\" name=\"desorption ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::FA)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000255\" name=\"flowing afterglow\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::FD)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000257\" name=\"field desorption\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::FI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000258\" name=\"field ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::GD_MS)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000259\" name=\"glow discharge ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::NICI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000271\" name=\"Negative ion chemical ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::NRMS)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000272\" name=\"neutralization reionization mass spectrometry\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::PI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000273\" name=\"photoionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::PYMS)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000274\" name=\"pyrolysis mass spectrometry\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::REMPI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000276\" name=\"resonance enhanced multiphoton ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::SELDI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000278\" name=\"surface enhanced laser desorption ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::SEND)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000279\" name=\"surface enhanced neat desorption\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::AI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000380\" name=\"adiabatic ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::ASI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000381\" name=\"associative ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::APPI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000382\" name=\"atmospheric pressure photoionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::AD)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000383\" name=\"autodetachment\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::AUI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000384\" name=\"autoionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::CEI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000385\" name=\"charge exchange ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::CHEMI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000386\" name=\"chemi-ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::SILI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000387\" name=\"desorption/ionization on silicon\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::DISSI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000388\" name=\"dissociative ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::EI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000389\" name=\"electron ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::LD)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000393\" name=\"laser desorption ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::LSI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000395\" name=\"liquid secondary ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::MESI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000397\" name=\"microelectrospray\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::NESI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000398\" name=\"nanoelectrospray\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::PEI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000399\" name=\"penning ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::PD)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000400\" name=\"plasma desorption ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::SI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000402\" name=\"secondary ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::SOI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000403\" name=\"soft ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::SPI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000404\" name=\"spark ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::SALDI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000405\" name=\"surface-assisted laser desorption ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::SUI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000406\" name=\"surface ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::TI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000407\" name=\"thermal ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::VI)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000408\" name=\"vertical ionization\"/>\n";
					}
					else if(so.getIonizationMethod()==IonSource::FIB)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000446\" name=\"fast ion bombardment\"/>\n";
					}
						
					writeUserParam_(os, so, 5);
					os  << "				</source>\n";				
				}
				//FORCED
				if (component_count<3 && in.getIonSources().size()==0)
				{
					os  << "				<source order=\"1234\">\n";
					os  << "					<cvParam cvRef=\"MS\" accession=\"MS:1000446\" name=\"fast ion bombardment\"/>\n";
					os  << "					<userParam name=\"warning\" type=\"xsd:string\" value=\"invented ion source, to fulfill mzML schema\" />\n";
					os  << "				</source>\n";				
				}
				//--------------------------------------------------------------------------------------------
				// mass analyzer
				//--------------------------------------------------------------------------------------------
				for (UInt i=0; i<in.getMassAnalyzers().size(); ++i)
				{
					const MassAnalyzer& ma = in.getMassAnalyzers()[i];
					os  << "				<analyzer order=\"" << ma.getOrder() << "\">\n";
					
					os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000014\" name=\"accuracy\" value=\"" << ma.getAccuracy() << "\"/>\n";
					os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000022\" name=\"TOF Total Path Length\" value=\"" << ma.getTOFTotalPathLength() << "\"/>\n";
					os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000024\" name=\"final MS exponent\" value=\"" << ma.getFinalMSExponent() << "\"/>\n";
					os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000025\" name=\"magnetic field strength\" value=\"" << ma.getMagneticFieldStrength() << "\"/>\n";
					
					if (ma.getReflectronState()==MassAnalyzer::ON)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000106\" name=\"reflectron on\"/>\n";
						
					}
					else if (ma.getReflectronState()==MassAnalyzer::OFF)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000105\" name=\"reflectron off\"/>\n";
					}
	
					if (ma.getType()==MassAnalyzer::FOURIERTRANSFORM)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000079\" name=\"fourier transform ion cyclotron resonance mass spectrometer\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::SECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000080\" name=\"magnetic sector\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::QUADRUPOLE)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000081\" name=\"quadrupole\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::TOF)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000084\" name=\"time-of-flight\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::ESA)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000254\" name=\"electrostatic energy analyzer\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::IT)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000264\" name=\"ion trap\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::SWIFT)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000284\" name=\"stored waveform inverse fourier transform\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::CYCLOTRON)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000288\" name=\"cyclotron\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::ORBITRAP)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000484\" name=\"orbitrap\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::AXIALEJECTIONLINEARIONTRAP)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000078\" name=\"axial ejection linear ion trap\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::PAULIONTRAP)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000082\" name=\"quadrupole ion trap\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::RADIALEJECTIONLINEARIONTRAP)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000083\" name=\"radial ejection linear ion trap\"/>\n";
					}
					else if (ma.getType()==MassAnalyzer::LIT)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000291\" name=\"linear ion trap\"/>\n";
					}
					
					writeUserParam_(os, ma, 5);
					os  << "				</analyzer>\n";				
				}
				//FORCED
				if (component_count<3 && in.getMassAnalyzers().size()==0)
				{
					os  << "				<analyzer order=\"1234\">\n";
					os << "						<cvParam cvRef=\"MS\" accession=\"MS:1000288\" name=\"cyclotron\"/>\n";
					os  << "					<userParam name=\"warning\" type=\"xsd:string\" value=\"invented mass analyzer, to fulfill mzML schema\" />\n";
					os  << "				</analyzer>\n";				
				}
				//--------------------------------------------------------------------------------------------
				// ion detector
				//--------------------------------------------------------------------------------------------
				for (UInt i=0; i<in.getIonDetectors().size(); ++i)
				{
					const IonDetector& id = in.getIonDetectors()[i];
					os  << "				<detector order=\"" << id.getOrder() << "\">\n";
	
					os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000028\" name=\"detector resolution\" value=\"" << id.getResolution() << "\"/>\n";
					os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000029\" name=\"sampling frequency\" value=\"" << id.getADCSamplingFrequency() << "\"/>\n";
	
					if (id.getAcquisitionMode()==IonDetector::ADC)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000117\" name=\"analog-digital converter\"/>\n";
					}
					else if (id.getAcquisitionMode()==IonDetector::PULSECOUNTING)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000118\" name=\"pulse counting\"/>\n";
					}
					else if (id.getAcquisitionMode()==IonDetector::TDC)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000119\" name=\"time-digital converter\"/>\n";
					}
					else if (id.getAcquisitionMode()==IonDetector::TRANSIENTRECORDER)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000120\" name=\"transient recorder\"/>\n";
					}
	
					if (id.getType()==IonDetector::CHANNELTRON)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000107\" name=\"channeltron\"/>\n";
					}
					else if (id.getType()==IonDetector::DALYDETECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000110\" name=\"daly detector\"/>\n";
					}
					else if (id.getType()==IonDetector::FARADAYCUP)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000112\" name=\"faraday cup\"/>\n";
					}
					else if (id.getType()==IonDetector::MICROCHANNELPLATEDETECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000114\" name=\"microchannel plate detector\"/>\n";
					}
					else if (id.getType()==IonDetector::MULTICOLLECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000115\" name=\"multi-collector\"/>\n";
					}
					else if (id.getType()==IonDetector::PHOTOMULTIPLIER)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000116\" name=\"photomultiplier\"/>\n";
					}
					else if (id.getType()==IonDetector::ELECTRONMULTIPLIER)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000253\" name=\"electron multiplier\"/>\n";
					}
					else if (id.getType()==IonDetector::ARRAYDETECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000345\" name=\"array detector\"/>\n";
					}
					else if (id.getType()==IonDetector::CONVERSIONDYNODE)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000346\" name=\"conversion dynode\"/>\n";
					}
					else if (id.getType()==IonDetector::DYNODE)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000347\" name=\"dynode\"/>\n";
					}
					else if (id.getType()==IonDetector::FOCALPLANECOLLECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000348\" name=\"focal plane collector\"/>\n";
					}
					else if (id.getType()==IonDetector::IONTOPHOTONDETECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000349\" name=\"ion-to-photon detector\"/>\n";
					}
					else if (id.getType()==IonDetector::POINTCOLLECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000350\" name=\"point collector\"/>\n";
					}
					else if (id.getType()==IonDetector::POSTACCELERATIONDETECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000351\" name=\"postacceleration detector\"/>\n";
					}
					else if (id.getType()==IonDetector::PHOTODIODEARRAYDETECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000621\" name=\"photodiode array detector\"/>\n";
					}
					else if (id.getType()==IonDetector::INDUCTIVEDETECTOR)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000624\" name=\"inductive detector\"/>\n";
					}
					else if (id.getType()==IonDetector::CONVERSIONDYNODEELECTRONMULTIPLIER)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000108\" name=\"conversion dynode electron multiplier\"/>\n";
					}
					else if (id.getType()==IonDetector::CONVERSIONDYNODEPHOTOMULTIPLIER)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000109\" name=\"conversion dynode photomultiplier\"/>\n";
					}
					else if (id.getType()==IonDetector::ELECTRONMULTIPLIERTUBE)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000111\" name=\"electron multiplier tube\"/>\n";
					}
					else if (id.getType()==IonDetector::FOCALPLANEARRAY)
					{
						os << "					<cvParam cvRef=\"MS\" accession=\"MS:1000113\" name=\"focal plane array\"/>\n";
					}
	
					writeUserParam_(os, id, 5);
					os  << "				</detector>\n";				
				}
				//FORCED
				if (component_count<3 && in.getIonDetectors().size()==0)
				{
					os  << "				<detector order=\"1234\">\n";
					os  << "					<cvParam cvRef=\"MS\" accession=\"MS:1000107\" name=\"channeltron\"/>\n";
					os  << "					<userParam name=\"warning\" type=\"xsd:string\" value=\"invented ion detector, to fulfill mzML schema\" />\n";
					os  << "				</detector>\n";				
				}
				os  << "			</componentList>\n";
			}
			os  << "			<softwareRef ref=\"so_in_0\"/>\n";
			os  << "		</instrumentConfiguration>\n";
			os  << "	</instrumentConfigurationList>\n";

			//--------------------------------------------------------------------------------------------
			// data processing
			//--------------------------------------------------------------------------------------------
			os  << "	<dataProcessingList count=\"1\">\n";			
			os  << "		<dataProcessing id=\"dp_ru_0\">\n";
			for (UInt i=0; i<exp.getDataProcessing().size(); ++i)
			{
				const DataProcessing& dp = exp.getDataProcessing()[i];
				os  << "			<processingMethod order=\"0\" softwareRef=\"so_dp_" << i << "\">\n";
				if (dp.getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000034\" name=\"charge deconvolution\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::DEISOTOPING)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000033\" name=\"deisotoping\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::SMOOTHING)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000592\" name=\"smoothing\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::CHARGE_CALCULATION)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000778\" name=\"charge state calculation\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::PRECURSOR_RECALCULATION)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000780\" name=\"precursor recalculation\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::BASELINE_REDUCTION)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000593\" name=\"baseline reduction\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::PEAK_PICKING)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000035\" name=\"peak picking\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::FEATURE_FINDING)==1)
				{
					//no CV term for this
				}
				if (dp.getProcessingActions().count(DataProcessing::ALIGNMENT)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000745\" name=\"retention time alignment\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::LOW_INTENSITY_REMOVAL)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000594\" name=\"low intensity data point removal\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::HIGH_INTENSITY_REMOVAL)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000746\" name=\"high intensity data point removal\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::CONVERSION_MZDATA)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000546\" name=\"Conversion to mzData\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::CONVERSION_MZML)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::CONVERSION_MZXML)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000545\" name=\"Conversion to mzXML\"/>\n";
				}
				if (dp.getProcessingActions().count(DataProcessing::CONVERSION_DTA)==1)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000741\" name=\"Conversion to dta\"/>\n";
				}
				//data processing attribute
				if (dp.getCompletionTime().isValid())
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000747\" name=\"completion time\" value=\"" << dp.getCompletionTime().toString("yyyy-MM-dd+hh:mm").toStdString() << "\"/>\n";
				}
				
				writeUserParam_(os, dp, 4);
				os  << "			</processingMethod>\n";
			}
			//FORCED (also includes a forced software: so_dp_0)
			if (exp.getDataProcessing().size()==0)
			{
				os  << "			<processingMethod order=\"0\" softwareRef=\"so_dp_0\">\n";
				os  << "				<cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\"/>\n";
				os  << "				<userParam name=\"warning\" type=\"xsd:string\" value=\"invented data processing, to fulfill mzML schema\" />\n";
				os  << "			</processingMethod>\n";
			}
			os  << "		</dataProcessing>\n";
			os  << "	</dataProcessingList>\n";		
			//--------------------------------------------------------------------------------------------
			// acquisitionSettings
			//--------------------------------------------------------------------------------------------
			
			//--------------------------------------------------------------------------------------------
			// run
			//--------------------------------------------------------------------------------------------
			os  << "	<run id=\"ru_0\" defaultInstrumentConfigurationRef=\"ic_0\" sampleRef=\"sa_0\"";
			if (exp.getDateTime().isValid())
			{
				os << " startTimeStamp=\"" << exp.getDateTime().get().substitute(' ','T') << "\"";
			}
			os  << ">\n";
			
			writeUserParam_(os, exp, 2);
			
			if (exp.getSourceFiles().size()!=0)
			{
				os	<< "		<sourceFileRefList count=\"" << exp.getSourceFiles().size() << "\">\n";
				for (UInt i=0; i<exp.getSourceFiles().size(); ++i)
				{
					os	<< "			<sourceFileRef ref=\"sf_ru_" << i << "\"/>\n";
				}
				os	<< "		</sourceFileRefList>\n";
			}

			os	<< "		<spectrumList count=\"" << exp.size() << "\" defaultDataProcessingRef=\"dp_ru_0\">\n"; 
			//--------------------------------------------------------------------------------------------
			//spectrum
			//--------------------------------------------------------------------------------------------
			for (UInt s=0; s<exp.size(); ++s)
			{
				const SpectrumType& spec = exp[s];
				os	<< "			<spectrum id=\"" << spec.getNativeID() << "\" index=\"" << s << "\" defaultArrayLength=\"" << spec.size() << "\"";
				if (spec.getSourceFile()!=SourceFile())
				{
					os << " sourceFileRef=\"sf_sp_" << s << "\"";
				}
				os  << ">\n";
				
				if (spec.getType()==SpectrumSettings::PEAKS)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid mass spectrum\"/>\n";
				}
				else if (spec.getType()==SpectrumSettings::RAWDATA)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000128\" name=\"profile mass spectrum\"/>\n";
				}
				else //FORCED
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000128\" name=\"profile mass spectrum\"/>\n";
				}
				if (spec.getMSLevel()!=0)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"" << spec.getMSLevel() << "\"/>\n";	
				}
				if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::FULL)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::ZOOM)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::SIM)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000582\" name=\"SIM spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::SRM)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000583\" name=\"SRM spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::CRM)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000581\" name=\"CRM spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::PRECURSOR)
				{
					os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000341\" name=\"precursor ion spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::CNG)
				{
					os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000325\" name=\"constant neutral gain spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::CNL)
				{
					os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000326\" name=\"constant neutral loss spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::PDA)
				{
					os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000620\" name=\"PDA spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::EMC)
				{
					os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000789\" name=\"enhanced multiply charged spectrum\" />\n";
				}
				else if (spec.getInstrumentSettings().getScanMode()==InstrumentSettings::TDF)
				{
					os	<< "			<cvParam cvRef=\"MS\" accession=\"MS:1000789\" name=\"time-delayed fragmentation spectrum\" />\n";
				}
				
				if (spec.getInstrumentSettings().getPolarity()==IonSource::NEGATIVE)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000129\" name=\"negative scan\"/>\n";
				}
				else if (spec.getInstrumentSettings().getPolarity()==IonSource::POSITIVE)
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\"/>\n";
				}
				else //FORCED
				{
					os << "				<cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\"/>\n";
				}
				os  << "				<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan time\" value=\"" << spec.getRT() << "\"/>\n";
				writeUserParam_(os, spec, 5);
				//--------------------------------------------------------------------------------------------
				//scan list
				//--------------------------------------------------------------------------------------------
				os	<< "				<scanList count=\"" << std::max((UInt)1,(UInt)spec.getAcquisitionInfo().size()) << "\">\n";
				ControlledVocabulary::CVTerm ai_term = getChildWithName_("MS:1000570",spec.getAcquisitionInfo().getMethodOfCombination());
				if (ai_term.id!="")
				{
					os  << "					<cvParam cvRef=\"MS\" accession=\"" << ai_term.id <<"\" name=\"" << ai_term.name << "\"/>\n";
				}
				else
				{
					os  << "					<cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\"/>\n";
				}
				writeUserParam_(os, spec.getAcquisitionInfo(), 5);
				//--------------------------------------------------------------------------------------------
				//scan
				//--------------------------------------------------------------------------------------------
				for (UInt j=0; j<spec.getAcquisitionInfo().size(); ++j)
				{
					const Acquisition& ac = spec.getAcquisitionInfo()[j];
					os	<< "					<scan number=\"" << ac.getNumber() << "\">\n";
					//cvParam: all stored in userParam
					writeUserParam_(os, ac, 6);
					//scan windows
					if (j==0 && spec.getInstrumentSettings().getScanWindows().size()!=0)
					{
						os	<< "					<scanWindowList count=\"" << spec.getInstrumentSettings().getScanWindows().size() << "\">\n";
						for (UInt j=0; j<spec.getInstrumentSettings().getScanWindows().size(); ++j)
						{
							os	<< "						<scanWindow>\n";
							os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000501\" name=\"scan m/z lower limit\" value=\"" << spec.getInstrumentSettings().getScanWindows()[j].begin << "\"/>\n";
							os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000500\" name=\"scan m/z upper limit\" value=\"" << spec.getInstrumentSettings().getScanWindows()[j].end << "\"/>\n";
							os	<< "						</scanWindow>\n";
						}
						os	<< "						</scanWindowList>\n";
					}
					os	<< "					</scan>\n";
				}
				if (spec.getAcquisitionInfo().size()==0)
				{
					os	<< "					<scan number=\"0\">\n";
					//scan windows
					if (spec.getInstrumentSettings().getScanWindows().size()!=0)
					{
						os	<< "					<scanWindowList count=\"" << spec.getInstrumentSettings().getScanWindows().size() << "\">\n";
						for (UInt j=0; j<spec.getInstrumentSettings().getScanWindows().size(); ++j)
						{
							os	<< "						<scanWindow>\n";
							os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000501\" name=\"scan m/z lower limit\" value=\"" << spec.getInstrumentSettings().getScanWindows()[j].begin << "\"/>\n";
							os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000500\" name=\"scan m/z upper limit\" value=\"" << spec.getInstrumentSettings().getScanWindows()[j].end << "\"/>\n";
							os	<< "						</scanWindow>\n";
						}
						os	<< "						</scanWindowList>\n";
					}
					os	<< "					</scan>\n";
				}
				os	<< "				</scanList>\n";
				//--------------------------------------------------------------------------------------------
				//precursor list
				//--------------------------------------------------------------------------------------------
				if (spec.getPrecursor() != Precursor() || spec.getPrecursorPeak() != typename SpectrumType::PrecursorPeakType())
				{
					os	<< "				<precursorList count=\"1\">\n";
					os	<< "					<precursor>\n";
					//--------------------------------------------------------------------------------------------
					//isolation window
					//--------------------------------------------------------------------------------------------
					//cvParam: nothing to do here for now
					//userParam: no extra object for it => no user paramters
					
					//--------------------------------------------------------------------------------------------
					//selected ion list
					//--------------------------------------------------------------------------------------------
					os	<< "						<selectedIonList count=\"1\">\n";
					os	<< "							<selectedIon>\n";
					os  << "								<cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected m/z\" value=\"" << spec.getPrecursorPeak().getMZ() << "\"/>\n";
					os  << "								<cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"" << spec.getPrecursorPeak().getCharge() << "\"/>\n";
					os  << "								<cvParam cvRef=\"MS\" accession=\"MS:1000042\" name=\"intensity\" value=\"" << spec.getPrecursorPeak().getIntensity() << "\"/>\n";
					for (UInt j=0; j<spec.getPrecursorPeak().getPossibleChargeStates().size(); ++j)
					{
						os  << "								<cvParam cvRef=\"MS\" accession=\"MS:1000633\" name=\"possible charge state\" value=\"" << spec.getPrecursorPeak().getPossibleChargeStates()[j] << "\"/>\n";
					}
					//userParam: no extra object for it => no user paramters
					os	<< "							</selectedIon>\n";					
					os	<< "						</selectedIonList>\n";

					//--------------------------------------------------------------------------------------------
					//activation
					//--------------------------------------------------------------------------------------------
					os	<< "						<activation>\n";
					os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000509\" name=\"activation energy\" value=\"" << spec.getPrecursor().getActivationEnergy() << "\"/>\n";
					if (spec.getPrecursor().getActivationMethod()==Precursor::CID)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000133\" name=\"collision-induced dissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::PD)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000134\" name=\"plasma desorption\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::PSD)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000135\" name=\"post-source decay\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::SID)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000136\" name=\"surface-induced dissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::BIRD)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000242\" name=\"blackbody infrared radiative dissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::ECD)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000250\" name=\"electron capture dissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::IMD)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000262\" name=\"infrared multiphoton dissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::SORI)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000282\" name=\"sustained off-resonance irradiation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::HCID)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000422\" name=\"high-energy collision-induced dissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::LCID)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000433\" name=\"low-energy collision-induced dissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::PHD)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000435\" name=\"photodissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::ETD)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000598\" name=\"electron transfer dissociation\"/>\n";
					}
					else if (spec.getPrecursor().getActivationMethod()==Precursor::PQD)
					{
						os  << "							<cvParam cvRef=\"MS\" accession=\"MS:1000599\" name=\"pulsed q dissociation\"/>\n";
					}
					//as "precursor" has no own user param it's userParam is stored here
					writeUserParam_(os, spec.getPrecursor(), 7);
					os	<< "						</activation>\n";
					os	<< "					</precursor>\n";					
					os	<< "				</precursorList>\n";
				}
				writeUserParam_(os, spec.getInstrumentSettings(), 5);
				//--------------------------------------------------------------------------------------------
				//binary data array list
				//--------------------------------------------------------------------------------------------
				if (spec.size()!=0)
				{
					String encoded_string;
					std::vector<DoubleReal> data_to_encode;
					os	<< "				<binaryDataArrayList count=\"" << (spec.getMetaDataArrays().size()+2) << "\">\n";
					//write m/z array
					data_to_encode.resize(spec.size());
					for (UInt p=0; p<spec.size(); ++p) data_to_encode[p] = spec[p].getMZ();
					decoder_.encode(data_to_encode, Base64::BYTEORDER_LITTLEENDIAN, encoded_string);
					os	<< "					<binaryDataArray encodedLength=\"" << encoded_string.size() << "\">\n";
					os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"/>\n";
					os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\n";
					os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n";
					os	<< "						<binary>" << encoded_string << "</binary>\n";
					os	<< "					</binaryDataArray>\n";
					//write intensity array
					for (UInt p=0; p<spec.size(); ++p) data_to_encode[p] = spec[p].getIntensity();
					decoder_.encode(data_to_encode, Base64::BYTEORDER_LITTLEENDIAN, encoded_string);
					os	<< "					<binaryDataArray encodedLength=\"" << encoded_string.size() << "\">\n";
					os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>\n";
					os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\n";
					os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n";
					os	<< "						<binary>" << encoded_string << "</binary>\n";
					os	<< "					</binaryDataArray>\n";
					//write meta data array
					for (UInt m=0; m<spec.getMetaDataArrays().size(); ++m)
					{
						const typename SpectrumType::MetaDataArray& array = spec.getMetaDataArrays()[m];
						data_to_encode.resize(array.size());
						for (UInt p=0; p<array.size(); ++p) data_to_encode[p] = array[p];
						decoder_.encode(data_to_encode, Base64::BYTEORDER_LITTLEENDIAN, encoded_string);
						os	<< "					<binaryDataArray arrayLength=\"" << array.size() << "\" encodedLength=\"" << encoded_string.size() << "\">\n";
						os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\n";
						os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n";
						ControlledVocabulary::CVTerm bi_term = getChildWithName_("MS:1000513",array.getName());
						if (bi_term.id!="")
						{
							os  << "						<cvParam cvRef=\"MS\" accession=\"" << bi_term.id <<"\" name=\"" << bi_term.name << "\"/>\n";
						}
						else
						{
							os  << "						<cvParam cvRef=\"MS\" accession=\"MS:1000786\" name=\"non-standard data array\" value=\"" << array.getName() << "\"/>\n";
						}
						writeUserParam_(os, array, 8);
						os	<< "						<binary>" << encoded_string << "</binary>\n";
						os	<< "					</binaryDataArray>\n";
					}
					os	<< "				</binaryDataArrayList>\n";
				}
				
				os	<< "			</spectrum>\n";
			}
			os	<< "		</spectrumList>\n";
			os  << "	</run>\n";		
			
			os	<< "</mzML>";
			logger_.endProgress();
		}

	} // namespace Internal
} // namespace OpenMS

#endif
