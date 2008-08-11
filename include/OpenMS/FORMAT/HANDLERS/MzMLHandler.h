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
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <sstream>

#include <iostream>

//TODO:
// - writing files, add to automatic tmp file validation in tests
// - units

//EXTEND:
// - Handle new enum Types when reading/writing mzXML, mzData and netCDF
// - chromatograms
// - acquisitionSettings
// - isolationWindow

//OUR MODEL:
// - remove ExperimentalSettings::type_ ?

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
					logger_(logger)
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
			Map<String, std::vector< std::pair<String,String> > > ref_param_;
			/// The source files: id => SourceFile
			Map<String, SourceFile> source_files_;
			/// The sample list: id => Sample
			Map<String, Sample> samples_;
			/// The software list: id => Software
			Map<String, Software> software_;
			/// The data processing list: id => Instrument
			Map<String, Instrument> instruments_;
			//@}

			/// Decoder/Encoder for Base64-data in MzML
			Base64 decoder_;
			
			/// Progress logger
			const ProgressLogger& logger_;
			
			///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
			ControlledVocabulary cv_;
			
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
			static const XMLCh* s_source_file_ref = xercesc::XMLString::transcode("sourceFileRef");
			//static const XMLCh* s_order = xercesc::XMLString::transcode("order");
			static const XMLCh* s_default_instrument_configuration_ref = xercesc::XMLString::transcode("defaultInstrumentConfigurationRef");
			
			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			
			//std::cout << "TAG: " << tag << std::endl;
			
			//determine parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);


			//determine the parent tag of the parent tag
			String parent_parent_tag;
			if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);
			
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
				DoubleReal double_version = 1.0;
				try
				{
					double_version = file_version.toDouble();
				}
				catch(...)
				{
					warning("Could not convert the mzML version string '" + file_version +"' to a double.");
				}
				if (double_version>version_.toDouble())
				{
					warning("The XML file (" + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
				}
			}
			else if (tag=="contact")
			{
				exp_->getContacts().push_back(ContactPerson());
			}
			//EXTEND "acquisition", "precursor" and "acquisition settings" can have a source file too
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
				//sample
				String sample_ref;
				if (optionalAttributeAsString_(sample_ref, attributes, s_sample_ref))
				{
					exp_->setSample(samples_[sample_ref]);
				}
				//instrument
				String instrument_ref = attributeAsString_(attributes, s_default_instrument_configuration_ref);
				exp_->setInstrument(instruments_[instrument_ref]);
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
			  //EXTEND the software should not be set here directly. It is determined through defaultInstrumentConfiguration.
			  //       As we do not have Software in Instrument yet, this hack is used...
			  //EXTEND "spectrum" and "binaryDataArray" also have a DataProcessingRef. What do we do with it?
				current_id_ = attributeAsString_(attributes, s_id);
				exp_->setSoftware(software_[attributeAsString_(attributes, s_software_ref)]);
			}

			else if (tag=="processingMethod")
			{
			  //EXTEND the processing should not be set here directly.
			  //       But where is it definded for the whole run? Ask the PSI people!
			  //       As we do not have Software in Instrument yet, this hack is used...
				//EXTEND Allow more then one processing step. Currently only the last one is stored
				//EXTEND Add order
				
			}
			else if (tag=="instrumentConfiguration")
			{
				current_id_ = attributeAsString_(attributes, s_id);
			}
			else if (tag=="softwareRef")
			{
				//EXTEND Add software to Instrument class
				//String ref = attributeAsString_(attributes, s_ref);
				//...
			}
			else if (tag=="source")
			{
				//EXTEND Allow several ion sources
				//EXTEND Add order to instrument components
				//Int order = attributeAsInt_(attributes, s_order);
				//...
			}
			else if (tag=="analyzer")
			{
				//EXTEND Add order to instrument components
				//Int order = attributeAsInt_(attributes, s_order);
				instruments_[current_id_].getMassAnalyzers().push_back(MassAnalyzer());
				
			}
			else if (tag=="detector")
			{
				//EXTEND Allow several detectors
				//EXTEND Add order to instrument components
				//Int order = attributeAsInt_(attributes, s_order);
				//...
			}
		}

		template <typename MapType>
		void MzMLHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static UInt scan_count=0;
			
			static const XMLCh* s_spectrum = xercesc::XMLString::transcode("spectrum");
			static const XMLCh* s_spectrum_list = xercesc::XMLString::transcode("spectrumList");
			static const XMLCh* s_mzml = xercesc::XMLString::transcode("mzML");

			//std::cout << "/TAG: " << open_tags_.back() << std::endl;
						
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
				instruments_.clear();
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
			//Warn when using obsolete CV terms
			if (cv_.exists(accession) && cv_.getTerm(accession).obsolete)
			{
				warning(String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'");
			}
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
				//EXTEND? Each acquisition can have all attributes like a scan (children of MS:1000503)
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
			//------------------------- instrumentConfiguration ----------------------------
			else if (parent_tag=="instrumentConfiguration")
			{
				//instrument model
				if (cv_.isChildOf(accession,"MS:1000031")) //is this a child of "instrument name"?
				{
					instruments_[current_id_].setName(cv_.getTerm(accession).name);
				}
				//instrument attribute
				else if (accession=="MS:1000529") //instrument serial number
				{
					//EXTEND
				}
				else if (accession=="MS:1000032") //customization
				{
					instruments_[current_id_].setCustomizations(value);
				}
				else if (accession=="MS:1000236") //transmission
				{
					//EXTEND
				}
				//ion optics type
				else if (accession=="MS:1000221") //magnetic deflection
				{
					//EXTEND
				}
				else if (accession=="MS:1000246") //delayed extraction
				{
					//EXTEND
				}
				else if (accession=="MS:1000275") //collision quadrupole
				{
					//EXTEND
				}
				else if (accession=="MS:1000281") //selected ion flow tube
				{
					//EXTEND
				}
				else if (accession=="MS:1000286") //time lag focusing
				{
					//EXTEND
				}
				else if (accession=="MS:1000300") //reflectron
				{
					//EXTEND
				}
				else if (accession=="MS:1000304") //accelerating voltage
				{
					//EXTEND
				}
				else if (accession=="MS:1000307") //einzel lens
				{
					//EXTEND
				}
				else if (accession=="MS:1000309") //first stability region
				{
					//EXTEND
				}
				else if (accession=="MS:1000310") //fringing field
				{
					//EXTEND
				}
				else if (accession=="MS:1000311") //kinetic energy analyzer
				{
					//EXTEND
				}
				else if (accession=="MS:1000320") //static field
				{
					//EXTEND
				}
				//ion optics attribute
				else if (accession=="MS:1000216") //field-free region
				{
					//EXTEND
				}
				else if (accession=="MS:1000308") //electric field strength
				{
					//EXTEND
				}
				else if (accession=="MS:1000319") //space charge effect
				{
					//EXTEND
				}
			}
			else if (parent_tag=="source")
			{
				//inlet type
				if (accession=="MS:1000055") //continuous flow fast atom bombardment
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::CONTINUOUSFLOWFASTATOMBOMBARDMENT);
				}
				else if (accession=="MS:1000056") //direct inlet
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::DIRECT);
				}
				else if (accession=="MS:1000057") //electrospray inlet
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::ELECTROSPRAYINLET);
				}
				else if (accession=="MS:1000058") //flow injection analysis
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::FLOWINJECTIONANALYSIS);
				}
				else if (accession=="MS:1000059") //inductively coupled plasma
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::INDUCTIVELYCOUPLEDPLASMA);
				}
				else if (accession=="MS:1000060") //infusion
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::INFUSION);
				}
				else if (accession=="MS:1000061") //jet separator
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::JETSEPARATOR);
				}
				else if (accession=="MS:1000062") //membrane separator
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::MEMBRANESEPARATOR);
				}
				else if (accession=="MS:1000063") //moving belt
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::MOVINGBELT);
				}
				else if (accession=="MS:1000064") //moving wire
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::MOVINGWIRE);
				}
				else if (accession=="MS:1000065") //open split
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::OPENSPLIT);
				}
				else if (accession=="MS:1000066") //particle beam
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::PARTICLEBEAM);
				}
				else if (accession=="MS:1000067") //reservoir
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::RESERVOIR);
				}
				else if (accession=="MS:1000068") //septum
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::SEPTUM);
				}
				else if (accession=="MS:1000069") //thermospray inlet
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::THERMOSPRAYINLET);
				}
				else if (accession=="MS:1000248") //direct insertion probe
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::BATCH);
				}
				else if (accession=="MS:1000249") //direct liquid introduction
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::CHROMATOGRAPHY);
				}
				else if (accession=="MS:1000396") //membrane inlet
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::MEMBRANE);
				}
				else if (accession=="MS:1000485") //nanospray inlet
				{
					instruments_[current_id_].getIonSource().setInletType(IonSource::NANOSPRAY);
				}
				//ionization type
				else if (accession=="MS:1000071") //chemical ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::CI);
				}
				else if (accession=="MS:1000073") //electrospray ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::ESI);
				}
				else if (accession=="MS:1000074") //fast atom bombardment ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::FAB);
				}
				else if (accession=="MS:1000227") //multiphoton ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::MPI);
				}
				else if (accession=="MS:1000240") //atmospheric pressure ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::API);
				}
				else if (accession=="MS:1000247") //desorption ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::DI);
				}
				else if (accession=="MS:1000255") //flowing afterglow
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::FA);
				}
				else if (accession=="MS:1000258") //field ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::FII);
				}
				else if (accession=="MS:1000259") //glow discharge ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::GD_MS);
				}
				else if (accession=="MS:1000271") //Negative ion chemical ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::NICI);
				}
				else if (accession=="MS:1000272") //neutralization reionization mass spectrometry
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::NRMS);
				}
				else if (accession=="MS:1000273") //photoionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::PI);
				}
				else if (accession=="MS:1000274") //pyrolysis mass spectrometry
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::PYMS);
				}
				else if (accession=="MS:1000276") //resonance enhanced multiphoton ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::REMPI);
				}
				else if (accession=="MS:1000380") //adiabatic ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::AI);
				}
				else if (accession=="MS:1000381") //associative ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::ASI);
				}
				else if (accession=="MS:1000383") //autodetachment
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::AD);
				}
				else if (accession=="MS:1000384") //autoionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::AUI);
				}
				else if (accession=="MS:1000385") //charge exchange ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::CEI);
				}
				else if (accession=="MS:1000386") //chemi-ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::CHEMI);
				}
				else if (accession=="MS:1000388") //dissociative ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::DISSI);
				}
				else if (accession=="MS:1000389") //electron ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::EI);
				}
				else if (accession=="MS:1000395") //liquid secondary ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::LSI);
				}
				else if (accession=="MS:1000399") //penning ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::PEI);
				}
				else if (accession=="MS:1000400") //plasma desorption ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::PD);
				}
				else if (accession=="MS:1000402") //secondary ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::SI);
				}
				else if (accession=="MS:1000403") //soft ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::SOI);
				}
				else if (accession=="MS:1000404") //spark ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::SPI);
				}
				else if (accession=="MS:1000406") //surface ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::SUI);
				}
				else if (accession=="MS:1000407") //thermal ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::TI);
				}
				else if (accession=="MS:1000408") //vertical ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::VI);
				}
				else if (accession=="MS:1000446") //fast ion bombardment
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::FIB);
				}
				else if (accession=="MS:1000070") //atmospheric pressure chemical ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::APCI);
				}
				else if (accession=="MS:1000239") //atmospheric pressure matrix-assisted laser desorption ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::AP_MALDI);
				}
				else if (accession=="MS:1000382") //atmospheric pressure photoionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::APPI);
				}
				else if (accession=="MS:1000075") //matrix-assisted laser desorption ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::MALDI);
				}
				else if (accession=="MS:1000257") //field desorption
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::FD);
				}
				else if (accession=="MS:1000387") //desorption/ionization on silicon
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::SILI);
				}
				else if (accession=="MS:1000393") //laser desorption ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::LD);
				}
				else if (accession=="MS:1000405") //surface-assisted laser desorption ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::SALDI);
				}
				else if (accession=="MS:1000397") //microelectrospray
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::MESI);
				}
				else if (accession=="MS:1000398") //nanoelectrospray
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::NESI);
				}
				else if (accession=="MS:1000278") //surface enhanced laser desorption ionization
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::SELDI);
				}
				else if (accession=="MS:1000279") //surface enhanced neat desorption
				{
					instruments_[current_id_].getIonSource().setIonizationMethod(IonSource::SEND);
				}
				//source attribute
				else if (accession=="MS:1000392") //ionization efficiency
				{
					//EXTEND
				}
				else if (accession=="MS:1000486") //source potential
				{
					//EXTEND
				}
				else if (accession=="MS:1000552") //maldi spot identifier
				{
					//EXTEND
				}
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
			}
			else if (parent_tag=="detector")
			{
				//detector type
				if (accession=="MS:1000107") //channeltron
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::CHANNELTRON);
				}
				else if (accession=="MS:1000110") //daly detector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::DALYDETECTOR);
				}
				else if (accession=="MS:1000112") //faraday cup
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::FARADAYCUP);
				}
				else if (accession=="MS:1000114") //microchannel plate detector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::MICROCHANNELPLATEDETECTOR);
				}
				else if (accession=="MS:1000115") //multi-collector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::MULTICOLLECTOR);
				}
				else if (accession=="MS:1000116") //photomultiplier
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::PHOTOMULTIPLIER);
				}
				else if (accession=="MS:1000253") //electron multiplier
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::ELECTRONMULTIPLIER);
				}
				else if (accession=="MS:1000345") //array detector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::ARRAYDETECTOR);
				}
				else if (accession=="MS:1000346") //conversion dynode
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::CONVERSIONDYNODE);
				}
				else if (accession=="MS:1000347") //dynode
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::DYNODE);
				}
				else if (accession=="MS:1000348") //focal plane collector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::FOCALPLANECOLLECTOR);
				}
				else if (accession=="MS:1000349") //ion-to-photon detector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::IONTOPHOTONDETECTOR);
				}
				else if (accession=="MS:1000350") //point collector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::POINTCOLLECTOR);
				}
				else if (accession=="MS:1000351") //postacceleration detector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::POSTACCELERATIONDETECTOR);
				}
				else if (accession=="MS:1000621") //photodiode array detector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::PHOTODIODEARRAYDETECTOR);
				}
				else if (accession=="MS:1000624") //inductive detector
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::INDUCTIVEDETECTOR);
				}
				else if (accession=="MS:1000108") //conversion dynode electron multiplier
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::CONVERSIONDYNODEELECTRONMULTIPLIER);
				}
				else if (accession=="MS:1000109") //conversion dynode photomultiplier
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::CONVERSIONDYNODEPHOTOMULTIPLIER);
				}
				else if (accession=="MS:1000111") //electron multiplier tube
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::ELECTRONMULTIPLIERTUBE);
				}
				else if (accession=="MS:1000113") //focal plane array
				{
					instruments_[current_id_].getIonDetector().setType(IonDetector::FOCALPLANEARRAY);
				}
				//detector attribute
				else if (accession=="MS:1000028") //detector resolution
				{
					instruments_[current_id_].getIonDetector().setResolution(value.toDouble());
				}
				else if (accession=="MS:1000029") //sampling frequency
				{
					instruments_[current_id_].getIonDetector().setADCSamplingFrequency(value.toDouble());
				}
				//dectector acquisition mode
				else if (accession=="MS:1000117") //analog-digital converter
				{
					instruments_[current_id_].getIonDetector().setAcquisitionMode(IonDetector::ADC);
				}
				else if (accession=="MS:1000118") //pulse counting
				{
					instruments_[current_id_].getIonDetector().setAcquisitionMode(IonDetector::PULSECOUNTING);
				}
				else if (accession=="MS:1000119") //time-digital converter
				{
					instruments_[current_id_].getIonDetector().setAcquisitionMode(IonDetector::TDC);
				}
				else if (accession=="MS:1000120") //transient recorder
				{
					instruments_[current_id_].getIonDetector().setAcquisitionMode(IonDetector::TRANSIENTRECORDER);
				} 
			}
			else if (parent_tag=="processingMethod")
			{
				//data processing parameter
				if (accession=="MS:1000629") //low intensity threshold
				{
					exp_->getProcessingMethod().setIntensityCutoff(value.toDouble());
				}
				else if (accession=="MS:1000631") //high intensity threshold
				{
					//EXTEND
				}
				//file format conversion
				else if (accession=="MS:1000544") //Conversion to mzML
				{
					//EXTEND
				}
				else if (accession=="MS:1000545") //Conversion to mzXML
				{
					//EXTEND
				}
				else if (accession=="MS:1000546") //Conversion to mzData
				{
					//EXTEND
				}
				//data processing action
				else if (accession=="MS:1000033") //deisotoping
				{
					exp_->getProcessingMethod().setDeisotoping(true);
				}
				else if (accession=="MS:1000034") //charge deconvolution
				{
					exp_->getProcessingMethod().setChargeDeconvolution(true);
				}
				else if (accession=="MS:1000035") //peak picking
				{
					//EXTEND
				}
				else if (accession=="MS:1000592") //smoothing
				{
					//EXTEND
				}
				else if (accession=="MS:1000593") //baseline reduction
				{
					//EXTEND
				}
				else if (accession=="MS:1000594") //low intensity data point removal
				{
					//EXTEND
				}
			}
		}//handleCVParam_



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
				instruments_[current_id_].getIonSource().setMetaValue(name,data_value);
			}
			else if (parent_tag=="analyzer")
			{
				instruments_[current_id_].getMassAnalyzers().back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="detector")
			{
				instruments_[current_id_].getIonDetector().setMetaValue(name,data_value);
			}
			else if (parent_tag=="sample")
			{
				samples_[current_id_].setMetaValue(name,data_value);
			}
			else if (parent_tag=="contact")
			{
				exp_->getContacts().back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="sourceFile")
			{
				//EXTEND Derive SourceFile from MetaInfoInterface
			}
			else if (parent_tag=="spectrum")
			{
				spec_.setMetaValue(name,data_value);
			}
			else if (parent_tag=="binaryDataArray")
			{
				if (data_.back().name=="mz" || data_.back().name=="int")
				{
					warning(String("Unhandled userParam in m/z or intensity binaryDataArray (name: '") + name + "' value: '" + value + "')");
				}
				else
				{
					data_.back().meta.setValue(name,data_value);
				}
			}
			else if (parent_tag=="spectrumDescription")
			{
				//EXTEND? Where should we put this?
			}
			else if (parent_tag=="scan")
			{
				spec_.getInstrumentSettings().setMetaValue(name,data_value);
			}
			else if (parent_tag=="acquisitionList")
			{
				//EXTEND Derive AcquisitionInfo from MetaDataInterface
			}
			else if (parent_tag=="acquisition")
			{
				spec_.getAcquisitionInfo().back().setMetaValue(name,data_value);
			}
			else if (parent_tag=="isolationWindow")
			{
				//EXTEND? Where should we put this?
			}
			else if (parent_tag=="selectedIon")
			{
				//EXTEND? Where should we put this?
			}
			else if (parent_tag=="activation")
			{
				//EXTEND? Where should we put this?
			}
			else if (parent_tag=="processingMethod")
			{
				exp_->getProcessingMethod().setMetaValue(name,data_value);
			}
		}//handleUserParam_

		template <typename MapType>
		void MzMLHandler<MapType>::writeTo(std::ostream& /*os*/)
		{
			throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}

	} // namespace Internal
} // namespace OpenMS

#endif
