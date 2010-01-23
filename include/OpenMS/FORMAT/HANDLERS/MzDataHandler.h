// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <sstream>

namespace OpenMS
{
	namespace Internal
	{
		/**
			@brief XML handler for MzDataFile
			
			MapType has to be a MSExperiment or have the same interface.
			Do not use this class. It is only needed in MzDataFile.
			
			@improvement Add implementation and tests of 'supDataArray' to store IntegerDataArray and StringDataArray of MSSpectrum (Hiwi)
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
					decoder_(),
					skip_spectrum_(false),
					logger_(logger)
	  	{
				init_();				
			}

      /// Constructor for a read-only handler
      MzDataHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger)
				: XMLHandler(filename, version),
					exp_(0),
					cexp_(&exp),
					peak_count_(0),
					meta_id_descs_(),
					decoder_(),
					skip_spectrum_(false),
					logger_(logger)
  		{
				init_();				
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
      virtual void characters(const XMLCh* const chars, const XMLSize_t length);

  		/// Writes the contents to a stream
			void writeTo(std::ostream& os);
		
			///Sets the options
			void setOptions(const PeakFileOptions& options)
			{ 
				options_ = options; 
			}
		
		private:
			void init_()
			{
				cv_terms_.resize(19);
				// SampleState
				String(";Solid;Liquid;Gas;Solution;Emulsion;Suspension").split(';',cv_terms_[0]);
				// IonizationMode
				String(";PositiveIonMode;NegativeIonMode").split(';',cv_terms_[1]);
				// ResolutionMethod
				String(";FWHM;TenPercentValley;Baseline").split(';',cv_terms_[2]);
				// ResolutionType
				String(";Constant;Proportional").split(';',cv_terms_[3]);
				// ScanFunction
				// is no longer used cv_terms_[4] is empty now
				// ScanDirection
				String(";Up;Down").split(';',cv_terms_[5]);
				// ScanLaw
				String(";Exponential;Linear;Quadratic").split(';',cv_terms_[6]);
				// PeakProcessing
				String(";CentroidMassSpectrum;ContinuumMassSpectrum").split(';',cv_terms_[7]);
				// ReflectronState
				String(";On;Off;None").split(';',cv_terms_[8]);
				// AcquisitionMode
				String(";PulseCounting;ADC;TDC;TransientRecorder").split(';',cv_terms_[9]);
				// IonizationType
				String(";ESI;EI;CI;FAB;TSP;LD;FD;FI;PD;SI;TI;API;ISI;CID;CAD;HN;APCI;APPI;ICP").split(';',cv_terms_[10]);
				// InletType
				String(";Direct;Batch;Chromatography;ParticleBeam;MembraneSeparator;OpenSplit;JetSeparator;Septum;Reservoir;MovingBelt;MovingWire;FlowInjectionAnalysis;ElectrosprayInlet;ThermosprayInlet;Infusion;ContinuousFlowFastAtomBombardment;InductivelyCoupledPlasma").split(';',cv_terms_[11]);
				// TandemScanningMethod
				// is no longer used cv_terms_[12] is empty now
				// DetectorType
				String(";EM;Photomultiplier;FocalPlaneArray;FaradayCup;ConversionDynodeElectronMultiplier;ConversionDynodePhotomultiplier;Multi-Collector;ChannelElectronMultiplier").split(';',cv_terms_[13]);
				// AnalyzerType
				String(";Quadrupole;PaulIonTrap;RadialEjectionLinearIonTrap;AxialEjectionLinearIonTrap;TOF;Sector;FourierTransform;IonStorage").split(';',cv_terms_[14]);
				// EnergyUnits
				// is no longer used cv_terms_[15] is empty now
				// ScanMode
				// is no longer used cv_terms_[16] is empty now
				// Polarity
				// is no longer used cv_terms_[17] is empty now
				// ActivationMethod
				String("CID;PSD;PD;SID").split(';',cv_terms_[18]);			
			}


		 protected:
      
      /// Peak type
      typedef typename MapType::PeakType PeakType;
      /// Spectrum type
			typedef MSSpectrum<PeakType> SpectrumType;
			
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
			//@}

			/// Decoder/Encoder for Base64-data in MzData
			Base64 decoder_;

			/// Flag that indicates whether this spectrum should be skipped (due to options)
			bool skip_spectrum_;
			
			/// Progress logger
			const ProgressLogger& logger_;

			/// fills the current spectrum with peaks and meta data
			void fillData_();

			///@name cvParam and userParam handling methods (for mzData and FeatureXML)
			//@{
			/**  
				@brief write cvParam containing strings to stream
				
				@p value string value
				@p acc accession number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
			*/
			inline void writeCVS_(std::ostream& os, DoubleReal value, const String& acc, const String& name, UInt indent=4) const
			{
				if (value!=0.0)
				{
					os << String(indent,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:" << acc << "\" name=\"" << name << "\" value=\"" << value << "\"/>\n";
				}
			}
			/**  
				@brief write cvParam containing strings to stream
				
				@p value string value
				@p acc accession number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
			*/
			inline void writeCVS_(std::ostream& os, const String& value, const String& acc, const String& name, UInt indent=4) const
			{
				if (value!="")
				{
					os << String(indent,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:" << acc << "\" name=\"" << name << "\" value=\"" << value << "\"/>\n";
				}
			}

			/**  
				@brief write cvParam element to stream

				@p os Output stream
				@p value enumeration value	
				@p map index if the terms in cv_terms_
				@p acc accession number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value=""/&gt;
			*/
			inline void writeCVS_(std::ostream& os, UInt value, UInt map, const String& acc, const String& name, UInt indent=4)
			{
				//abort when receiving a wrong map index
				if (map>=cv_terms_.size())
				{
					warning(STORE, String("Cannot find map '") + map + "' needed to write CV term '" + name + "' with accession '" + acc + "'.");
					return;
				}
				//abort when receiving a wrong term index
				if (value>=cv_terms_[map].size())
				{
					warning(STORE, String("Cannot find value '") + value + "' needed to write CV term '" + name + "' with accession '" + acc + "'.");
					return;
				}
				writeCVS_(os, cv_terms_[map][value], acc, name, indent);
			}
			
			///Writing the MetaInfo as UserParam to the file
			inline void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent=4)
			{
				std::vector<String> keys;
				meta.getKeys(keys);
				for (std::vector<String>::const_iterator it = keys.begin(); it!=keys.end(); ++it)
				{
					if ( (*it)[0] != '#')  // internally used meta info start with '#'
					{
						os << String(indent,'\t') << "<userParam name=\"" << *it << "\" value=\"" << meta.getMetaValue(*it) << "\"/>\n";
					}
				}
			}
			/** 
				@brief read attributes of MzData's cvParamType

				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:1000001" name="@p name" value="@p value"/&gt;
				@p name and sometimes @p value are defined in the MzData ontology.
			*/
			void cvParam_(const String& name, const String& value);
			//@}

			/**
				@brief write binary data to stream (first one)
			
				The @p name and @p id are only used if the @p tag is @em supDataArrayBinary or @em supDataArray.
			*/
			inline void writeBinary_(std::ostream& os, Size size, const String& tag, const String& name="", SignedSize id=-1)
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

				String str;				
				decoder_.encode(data_to_encode_, Base64::BYTEORDER_LITTLEENDIAN, str);
				data_to_encode_.clear();
				os << "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\""
					 << size << "\">"
					 << str
					 << "</data>\n\t\t\t</" << tag << ">\n";
			}
			
			//Data processing auxilary variable
			DataProcessing data_processing_;
			
		};

		//--------------------------------------------------------------------------------

		template <typename MapType>
		void MzDataHandler<MapType>::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
		{
			// skip current spectrum
			if (skip_spectrum_) return;
			
			char* transcoded_chars = sm_.convert(chars);
			
			//current tag
			const String& current_tag = open_tags_.back();
			
			//determine the parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
			
			
			if (current_tag=="sampleName")
			{
				exp_->getSample().setName( sm_.convert(chars) );
			}
			else if (current_tag=="instrumentName")
			{
				exp_->getInstrument().setName(sm_.convert(chars));
			}
			else if (current_tag=="version")
			{
				data_processing_.getSoftware().setVersion( sm_.convert(chars) );
			}
			else if (current_tag=="institution")
			{
				exp_->getContacts().back().setInstitution( sm_.convert(chars) );
			}
			else if (current_tag=="contactInfo")
			{
				exp_->getContacts().back().setContactInfo( sm_.convert(chars) );
			}
			else if (current_tag=="name" && parent_tag=="contact")
			{
				exp_->getContacts().back().setName(sm_.convert(chars));
			}
			else if (current_tag=="name" && parent_tag=="software")
			{
				data_processing_.getSoftware().setName( sm_.convert(chars) );
			}
			else if (current_tag=="comments" && parent_tag=="software")
			{
				data_processing_.getSoftware().setMetaValue("comment", String(sm_.convert(chars)) );
			}
			else if (current_tag == "comments" && parent_tag=="spectrumDesc")
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
				spec_.getFloatDataArrays().back().setName(transcoded_chars);
			}
			else if (current_tag=="nameOfFile" && parent_tag == "sourceFile")
			{
				exp_->getSourceFiles().back().setNameOfFile( sm_.convert(chars) );
			}
			else if (current_tag == "nameOfFile" && parent_tag == "supSourceFile")
			{
				//ignored
			}
			else if (current_tag=="pathToFile" && parent_tag == "sourceFile")
			{
				exp_->getSourceFiles().back().setPathToFile( sm_.convert(chars) );
			}
			else if (current_tag == "pathToFile" && parent_tag == "supSourceFile")
			{
				//ignored
			}
			else if (current_tag=="fileType" && parent_tag == "sourceFile")
			{
				exp_->getSourceFiles().back().setFileType( sm_.convert(chars) );
			}
			else if (current_tag == "fileType" && parent_tag == "supSourceFile")
			{
				//ignored
			}
			else
			{
				String trimmed_transcoded_chars = transcoded_chars;
				trimmed_transcoded_chars.trim();
				if (trimmed_transcoded_chars!="")
				{
					warning(LOAD, String("Unhandled character content in tag '") + current_tag + "': " + trimmed_transcoded_chars);
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
			static const XMLCh* s_accessionnumber = xercesc::XMLString::transcode("accessionNumber");

			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);
			//std::cout << "Start: '" << tag << "'" << std::endl;
			
			//determine the parent tag
			String parent_tag;
			if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
							
			//do nothing until a new spectrum is reached
			if (tag!="spectrum" && skip_spectrum_) return;


			// Do something depending on the tag
			if (tag=="sourceFile")
			{
				exp_->getSourceFiles().push_back(SourceFile());
			}
			if (tag=="contact")
			{
				exp_->getContacts().resize(exp_->getContacts().size()+1);
			}
			else if (tag=="source")
			{
				exp_->getInstrument().getIonSources().resize(1);
			}
			else if (tag=="detector")
			{
				exp_->getInstrument().getIonDetectors().resize(1);
			}
			else if (tag=="analyzer")
			{
				exp_->getInstrument().getMassAnalyzers().resize(exp_->getInstrument().getMassAnalyzers().size()+1);
			}
			else if (tag=="software")
			{
				data_processing_ = DataProcessing();
				if (attributes.getIndex(sm_.convert("completionTime"))!=-1)
				{
					data_processing_.setCompletionTime( asDateTime_(sm_.convert(attributes.getValue(sm_.convert("completionTime")))) );
				}
			}
			else if (tag=="precursor")
			{
				spec_.getPrecursors().push_back(Precursor());
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
				if (optionalAttributeAsString_(comment, attributes, s_comment))
				{
					meta_id_descs_.back().second.setMetaValue("comment",comment);
				}
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
					spec_.getPrecursors().back().setMetaValue(name, value);
				}
				else if (parent_tag=="activation")
				{
					spec_.getPrecursors().back().setMetaValue(name, value);
				}
				else if (parent_tag=="supDataDesc")
				{
					meta_id_descs_.back().second.setMetaValue(name, value);
				}
				else if (parent_tag=="detector")
				{
					exp_->getInstrument().getIonDetectors().back().setMetaValue(name,value);
				}
				else if (parent_tag=="source")
				{
					exp_->getInstrument().getIonSources().back().setMetaValue(name,value);
				}
				else if (parent_tag=="sampleDescription")
				{
					exp_->getSample().setMetaValue(name,value);
				}
				else if (parent_tag=="analyzer")
				{
					exp_->getInstrument().getMassAnalyzers().back().setMetaValue(name,value);
				}
				else if (parent_tag=="additional")
				{
					exp_->getInstrument().setMetaValue(name,value);
				}
				else if (parent_tag=="processingMethod")
				{
					data_processing_.setMetaValue(name,value);			
				}
				else
				{
					warning(LOAD, "Invalid userParam: name=\"" + name + ", value=\"" + value + "\"");
				}
			}
			else if (tag=="supDataArrayBinary")
			{
				
				//create FloatDataArray
				typename MapType::SpectrumType::FloatDataArray mda;
				//Assign the right MetaInfoDescription ("supDesc" tag)
				String id = attributeAsString_(attributes, s_id);
				for (Size i=0;i<meta_id_descs_.size(); ++i)
				{
					if (meta_id_descs_[i].first==id)
					{
						mda.MetaInfoDescription::operator=(meta_id_descs_[i].second);
						break;
					}
				}
				//append FloatDataArray
				spec_.getFloatDataArrays().push_back(mda);
			}
			else if (tag=="spectrum")
			{
				spec_ = SpectrumType();
				spec_.setNativeID(String("spectrum=") + attributeAsString_(attributes, s_id));
				spec_.getDataProcessing().push_back(data_processing_);
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
			else if (tag=="mzData")
			{
				//handle file id
				exp_->setIdentifier(attributeAsString_(attributes, s_accessionnumber));
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
					warning(LOAD, String("Invalid spectrum type '") + tmp_type + "'.");
				}
				
				spec_.getAcquisitionInfo().setMethodOfCombination(attributeAsString_(attributes, s_methodofcombination));
			}
			else if (tag=="acquisition")
			{
				spec_.getAcquisitionInfo().insert(spec_.getAcquisitionInfo().end(), Acquisition());
				spec_.getAcquisitionInfo().back().setIdentifier(attributeAsString_(attributes, s_acqnumber));
			}
			else if (tag=="spectrumInstrument" || tag=="acqInstrument")
			{
				spec_.setMSLevel(attributeAsInt_(attributes, s_mslevel));
				ScanWindow window;
				optionalAttributeAsDouble_(window.begin, attributes, s_mzrangestart);
				optionalAttributeAsDouble_(window.end, attributes, s_mzrangestop);
				if (window.begin!=0.0 || window.end!=0.0)
				{
					spec_.getInstrumentSettings().getScanWindows().push_back(window);
				}
				
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
					spec_.reserve(peak_count_);
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
			else if (tag=="arrayName" && parent_tag=="supDataArrayBinary")
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
			
			static const XMLCh* s_spectrum = xercesc::XMLString::transcode("spectrum");
			static const XMLCh* s_mzdata = xercesc::XMLString::transcode("mzData");
			
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
		void MzDataHandler<MapType>::fillData_()
		{
			std::vector<Real> decoded;
			std::vector<DoubleReal> decoded_double;
			
			// data_to_decode is an encoded spectrum, represented as
			// vector of base64-encoded strings:
			// Each string represents one property (e.g. mzData) and decodes
			// to a vector of property values - one value for every peak in the spectrum.
			for (Size i=0; i<data_to_decode_.size(); ++i)
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

			// this works only if MapType::PeakType is a Peak1D or derived from it
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
				for (Size i=0;i<spec_.getFloatDataArrays().size();++i)
				{
					spec_.getFloatDataArrays()[i].reserve(peak_count_);
				}
				
				//push_back the peaks into the container				
				for (Size n = 0 ; n < peak_count_ ; ++n)
				{
					DoubleReal mz = mz_precision_64 ? decoded_double_list_[0][n] : decoded_list_[0][n];
					DoubleReal intensity = int_precision_64 ? decoded_double_list_[1][n] : decoded_list_[1][n];
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(mz)))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(intensity))))
					{
						PeakType tmp;
						tmp.setIntensity(intensity);
						tmp.setMZ(mz);
						spec_.push_back(tmp);
						//load data from meta data arrays
						for (Size i=0;i<spec_.getFloatDataArrays().size();++i)
						{
							spec_.getFloatDataArrays()[i].push_back(precisions_[2+i]=="64" ? decoded_double_list_[2+i][n] : decoded_list_[2+i][n]);
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
				 << "<mzData version=\"1.05\" accessionNumber=\"" << cexp_->getIdentifier() << "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"http://psidev.sourceforge.net/ms/xml/mzdata/mzdata.xsd\">\n";
			
			//---------------------------------------------------------------------------------------------------
			//DESCRIPTION
			const Sample& sm = cexp_->getSample();
			os << "\t<description>\n"
				 << "\t\t<admin>\n"
				 << "\t\t\t<sampleName>"
				 << sm.getName()
				 << "</sampleName>\n";

			if( sm.getNumber()!="" || sm.getState() || sm.getMass() || sm.getVolume() || sm.getConcentration()	|| !sm.isMetaEmpty())
			{
				os << "\t\t\t<sampleDescription>\n";
				writeCVS_(os, sm.getNumber(), "1000001", "SampleNumber");
				writeCVS_(os, sm.getState(), 0, "1000003", "SampleState");
				writeCVS_(os, sm.getMass(), "1000004", "SampleMass");
				writeCVS_(os, sm.getVolume(), "1000005", "SampleVolume");
				writeCVS_(os, sm.getConcentration(), "1000006", "SampleConcentration");
				writeUserParam_(os, cexp_->getSample());
				os << "\t\t\t</sampleDescription>\n";
			}

			if (cexp_->getSourceFiles().size()>=1)
			{
				os << "\t\t\t<sourceFile>\n"
					 << "\t\t\t\t<nameOfFile>" << cexp_->getSourceFiles()[0].getNameOfFile() << "</nameOfFile>\n"
					 << "\t\t\t\t<pathToFile>" << cexp_->getSourceFiles()[0].getPathToFile() << "</pathToFile>\n";
				if (cexp_->getSourceFiles()[0].getFileType()!="")
					os << "\t\t\t\t<fileType>" << cexp_->getSourceFiles()[0].getFileType() << "</fileType>\n";
				os << "\t\t\t</sourceFile>\n";
			}
			if (cexp_->getSourceFiles().size()>1)
			{
				warning(STORE, "The MzData format can store only one source file. Only the first one is stored!");
			}

			for (Size i=0; i < cexp_->getContacts().size(); ++i)
			{
				os << "\t\t\t<contact>\n"
					 << "\t\t\t\t<name>" << cexp_->getContacts()[i].getFirstName() << " " << cexp_->getContacts()[i].getLastName() << "</name>\n"
					 << "\t\t\t\t<institution>" << cexp_->getContacts()[i].getInstitution() << "</institution>\n";
				if (cexp_->getContacts()[i].getContactInfo()!="")
					os << "\t\t\t\t<contactInfo>" << cexp_->getContacts()[i].getContactInfo() << "</contactInfo>\n";
				os << "\t\t\t</contact>\n";
			}
			//no contacts given => add empty entry as there must be a contact entry
			if (cexp_->getContacts().size()==0)
			{
				os << "\t\t\t<contact>\n"
					 << "\t\t\t\t<name></name>\n"
					 << "\t\t\t\t<institution></institution>\n";
				os << "\t\t\t</contact>\n";
			}
			
			os << "\t\t</admin>\n";
			const Instrument& inst = cexp_->getInstrument();
			os << "\t\t<instrument>\n"
				 << "\t\t\t<instrumentName>" << inst.getName() << "</instrumentName>\n"
				 << "\t\t\t<source>\n";
			if (inst.getIonSources().size()>=1)
			{
				writeCVS_(os, inst.getIonSources()[0].getInletType(),11, "1000007", "InletType");
				writeCVS_(os, inst.getIonSources()[0].getIonizationMethod(), 10, "1000008","IonizationType");
				writeCVS_(os, inst.getIonSources()[0].getPolarity(), 1, "1000009", "IonizationMode");
				writeUserParam_(os, inst.getIonSources()[0]);
			}
			if (inst.getIonSources().size()>1)
			{
				warning(STORE, "The MzData format can store only one ion source. Only the first one is stored!");
			}
			os << "\t\t\t</source>\n";
						
			//no analyzer given => add empty entry as there must be one entry
			if (inst.getMassAnalyzers().size() == 0)
			{
				os << "\t\t\t<analyzerList count=\"1\">\n"
				   << "\t\t\t\t<analyzer>\n"
				   << "\t\t\t\t</analyzer>\n";
			}
			else
			{
				os << "\t\t\t<analyzerList count=\"" << inst.getMassAnalyzers().size() << "\">\n";
				for (Size i=0; i<inst.getMassAnalyzers().size(); ++i)
				{
					os << "\t\t\t\t<analyzer>\n";
					const MassAnalyzer& ana = inst.getMassAnalyzers()[i];
					writeCVS_(os, ana.getType(), 14, "1000010", "AnalyzerType",5);
					writeCVS_(os, ana.getResolution(), "1000011", "MassResolution",5);
					writeCVS_(os, ana.getResolutionMethod(), 2,"1000012", "ResolutionMethod",5);
					writeCVS_(os, ana.getResolutionType(), 3, "1000013", "ResolutionType",5);
					writeCVS_(os, ana.getAccuracy(), "1000014", "Accuracy",5);
					writeCVS_(os, ana.getScanRate(), "1000015", "ScanRate",5);
					writeCVS_(os, ana.getScanTime(), "1000016", "ScanTime",5);
					writeCVS_(os, ana.getScanDirection(), 5,	"1000018", "ScanDirection",5);
					writeCVS_(os, ana.getScanLaw(), 6, "1000019", "ScanLaw",5);
					writeCVS_(os, ana.getReflectronState(), 8, "1000021", "ReflectronState",5);
					writeCVS_(os, ana.getTOFTotalPathLength(), "1000022", "TOFTotalPathLength",5);
					writeCVS_(os, ana.getIsolationWidth(), "1000023", "IsolationWidth",5);
					writeCVS_(os, ana.getFinalMSExponent(), "1000024", "FinalMSExponent",5);
					writeCVS_(os, ana.getMagneticFieldStrength(), "1000025", "MagneticFieldStrength",5);
					writeUserParam_(os, ana, 5);
					os << "\t\t\t\t</analyzer>\n";
				}
			}
			os << "\t\t\t</analyzerList>\n";
			
			os << "\t\t\t<detector>\n";
			if (inst.getIonDetectors().size()>=1)
			{
				writeCVS_(os, inst.getIonDetectors()[0].getType(), 13, "1000026", "DetectorType");
				writeCVS_(os, inst.getIonDetectors()[0].getAcquisitionMode(), 9, "1000027", "DetectorAcquisitionMode");
				writeCVS_(os, inst.getIonDetectors()[0].getResolution(), "1000028", "DetectorResolution");
				writeCVS_(os, inst.getIonDetectors()[0].getADCSamplingFrequency(), "1000029", "SamplingFrequency");
				writeUserParam_(os, inst.getIonDetectors()[0]);
			}
			if (inst.getIonDetectors().size()>1)
			{
				warning(STORE, "The MzData format can store only one ion detector. Only the first one is stored!");
			}
			os << "\t\t\t</detector>\n";
			if (inst.getVendor()!="" || inst.getModel()!="" || inst.getCustomizations()!="")
			{
				os << "\t\t\t<additional>\n";
				writeCVS_(os, inst.getVendor(), "1000030", "Vendor");
				writeCVS_(os, inst.getModel(), "1000031", "Model");
				writeCVS_(os, inst.getCustomizations(), "1000032", "Customization");
				writeUserParam_(os, inst);
				os << "\t\t\t</additional>\n";
			}
			os << "\t\t</instrument>\n";
			
			//the data processing information of the first spectrum is used for the whole file
			if (cexp_->size()==0 || (*cexp_)[0].getDataProcessing().size()==0)
			{
				os << "\t\t<dataProcessing>\n"
					 << "\t\t\t<software>\n"
					 << "\t\t\t\t<name></name>\n"
					 << "\t\t\t\t<version></version>\n"
					 << "\t\t\t</software>\n"
					 << "\t\t</dataProcessing>\n";
			}
			else
			{
				const DataProcessing& data_processing = (*cexp_)[0].getDataProcessing()[0];
				os << "\t\t<dataProcessing>\n"
					 << "\t\t\t<software";
				if (data_processing.getCompletionTime()!=DateTime())
				{
					os << " completionTime=\"" << data_processing.getCompletionTime().get().substitute(' ','T') << "\"";
				}
				os << ">\n"
					 << "\t\t\t\t<name>" << data_processing.getSoftware().getName() << "</name>\n"
					 << "\t\t\t\t<version>" << data_processing.getSoftware().getVersion() << "</version>\n";
				os << "\t\t\t</software>\n"
					 << "\t\t\t<processingMethod>\n";
				if(data_processing.getProcessingActions().count(DataProcessing::DEISOTOPING)==1)
				{
					os << "\t\t\t\t<cvParam cvLabel=\"psi\" name=\"Deisotoping\" accession=\"PSI:1000033\" />\n";
				}
				if(data_processing.getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION)==1)
				{
					os << "\t\t\t\t<cvParam cvLabel=\"psi\" name=\"ChargeDeconvolution\" accession=\"PSI:1000034\" />\n";
				}
				if(data_processing.getProcessingActions().count(DataProcessing::PEAK_PICKING)==1)
				{
					os << "\t\t\t\t<cvParam cvLabel=\"psi\" name=\"Centroid Mass Spectrum\" accession=\"PSI:1000127\"/>\n";
				}
				writeUserParam_(os, data_processing);
				os << "\t\t\t</processingMethod>\n"
					 << "\t\t</dataProcessing>\n";
			}
			os << "\t</description>\n";

			//---------------------------------------------------------------------------------------------------
			//ACTUAL DATA
			if (cexp_->size()!=0)
			{
				//check if the nativeID of all spectra are numbers or numbers prefixed with 'spectrum='
				//If not we need to renumber all spectra.
				bool all_numbers = true;
				bool all_empty = true;
				bool all_prefixed_numbers = true;
				for (Size s=0; s<cexp_->size(); s++)
				{
					String native_id = (*cexp_)[s].getNativeID();
					if (!native_id.hasPrefix("spectrum="))
					{
						all_prefixed_numbers = false;
					}
					else
					{
						native_id = native_id.substr(9);
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
					warning(STORE, "Not all spectrum native IDs are numbers or correctly prefixed with 'spectrum='. The spectra are renumbered and the native IDs are lost!");
				}
				//Map to store the last spectrum ID for each MS level (needed to find precursor spectra)
				Map<Int,Size> level_id; 
				
				os << "\t<spectrumList count=\"" << cexp_->size() << "\">\n";
				for (Size s=0; s<cexp_->size(); ++s)
				{
					logger_.setProgress(s);
					const SpectrumType& spec = (*cexp_)[s];
					
					Size spectrum_id = s+1;
					if (all_prefixed_numbers)
					{
						spectrum_id = spec.getNativeID().substr(9).toInt();
					}
					else if (all_numbers)
					{
						spectrum_id = spec.getNativeID().toInt();
					}
					os << "\t\t<spectrum id=\"" << spectrum_id << "\">\n"
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
						else
						{
							warning(STORE, "Spectrum type is unknown, assuming 'discrete'");
							os << "discrete";
						}
	
						os << "\" methodOfCombination=\"" << spec.getAcquisitionInfo().getMethodOfCombination() << "\""
						   << " count=\"" << spec.getAcquisitionInfo().size() << "\">\n";
						for (Size i=0; i<spec.getAcquisitionInfo().size(); ++i)
						{
							const Acquisition& ac = spec.getAcquisitionInfo()[i];
							Int acq_number = 0;
							try
							{
								if (ac.getIdentifier() != "")
								{
									acq_number =  ac.getIdentifier().toInt();
								}
							}
							catch(...)
							{
								warning(STORE, String("Could not convert acquisition identifier '") + ac.getIdentifier() + "' to an integer. Using '0' instead!");
								acq_number = 0;
							}
							os << "\t\t\t\t\t\t<acquisition acqNumber=\"" << acq_number << "\">\n";
							writeUserParam_(os, ac, 7);
							os << "\t\t\t\t\t\t</acquisition>\n";
						}
						os << "\t\t\t\t\t</acqSpecification>\n";
					}
	
					const InstrumentSettings& iset = spec.getInstrumentSettings();
					os << "\t\t\t\t\t<spectrumInstrument msLevel=\"" << spec.getMSLevel() << "\"";
					level_id[spec.getMSLevel()] = spectrum_id;
					
					if (iset.getScanWindows().size() > 0)
					{
						os << " mzRangeStart=\"" << iset.getScanWindows()[0].begin << "\" mzRangeStop=\"" << iset.getScanWindows()[0].end << "\"";
					}
					if (iset.getScanWindows().size() > 1)
					{
						warning(STORE, "The MzData format can store only one scan window for each scan. Only the first one is stored!");
					}
					os << ">\n";
	
					//scan mode
					switch(iset.getScanMode())
					{
						case InstrumentSettings::UNKNOWN:
							//do nothing here
							break;
						case InstrumentSettings::MASSSPECTRUM:
							if (iset.getZoomScan())
							{
								os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"Zoom\"/>\n";
							}
							else
							{
								os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"MassScan\"/>\n";
							}
							break;
						case InstrumentSettings::SIM:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"SelectedIonDetection\"/>\n";
							break;
						case InstrumentSettings::SRM:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"SelectedReactionMonitoring\"/>\n";
							break;
						case InstrumentSettings::CRM:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"ConsecutiveReactionMonitoring\"/>\n";
							break;
						case InstrumentSettings::CNG:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"ConstantNeutralGainScan\"/>\n";
							break;
						case InstrumentSettings::CNL:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"ConstantNeutralLossScan\"/>\n";
							break;
						case InstrumentSettings::PRECURSOR:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"PrecursorIonScan\"/>\n";
							break;
						case InstrumentSettings::ABSORBTION:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"PhotodiodeArrayDetector\"/>\n";
							break;
						case InstrumentSettings::EMC:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"EnhancedMultiplyChargedScan\"/>\n";
							break;
						case InstrumentSettings::TDF:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"TimeDelayedFragmentationScan\"/>\n";
							break;
						default:
							os << "						<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"MassScan\"/>\n";
							warning(STORE, String("Scan mode '") + InstrumentSettings::NamesOfScanMode[iset.getScanMode()] + "' not supported by mzData. Using 'MassScan' scan mode!");
					}
					
				 	//scan polarity
					if (spec.getInstrumentSettings().getPolarity()==IonSource::POSITIVE)
					{
						os << String(6,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:1000037\" name=\"Polarity\" value=\"Positive\"/>\n";
					}
					else if (spec.getInstrumentSettings().getPolarity()==IonSource::NEGATIVE)
					{
						os << String(6,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:1000037\" name=\"Polarity\" value=\"Negative\"/>\n";
					}
					
					//Retiontion time already in TimeInSeconds
					writeCVS_(os, spec.getRT(), "1000039", "TimeInSeconds",6);
					writeUserParam_(os, spec.getInstrumentSettings(), 6);
					os 	<< "\t\t\t\t\t</spectrumInstrument>\n\t\t\t\t</spectrumSettings>\n";
	
					if (spec.getPrecursors().size()!=0)
					{
						Int precursor_ms_level = spec.getMSLevel()-1;
						SignedSize precursor_id = -1;
						if (level_id.has(precursor_ms_level))
						{
							precursor_id = level_id[precursor_ms_level];
						}
						os << "\t\t\t\t<precursorList count=\"" << spec.getPrecursors().size() << "\">\n";
						for (Size i=0; i< spec.getPrecursors().size(); ++i)
						{
							const Precursor& precursor = spec.getPrecursors()[i];
							os << "\t\t\t\t\t<precursor msLevel=\"" << precursor_ms_level << "\" spectrumRef=\"" << precursor_id << "\">\n";
							os << "\t\t\t\t\t\t<ionSelection>\n";
							if (precursor != Precursor())
							{
								writeCVS_(os, precursor.getMZ(), "1000040", "MassToChargeRatio",7);
								writeCVS_(os, precursor.getCharge(), "1000041", "ChargeState",7);
								writeCVS_(os, precursor.getIntensity(), "1000042", "Intensity",7);
								os << "\t\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000043\" name=\"IntensityUnit\" value=\"NumberOfCounts\"/>\n";
								writeUserParam_(os, precursor, 7);
							}
							os << "\t\t\t\t\t\t</ionSelection>\n";
							os << "\t\t\t\t\t\t<activation>\n";
							if (precursor != Precursor())
							{
								if (precursor.getActivationMethods().size()>0)
								{
									writeCVS_(os, *(precursor.getActivationMethods().begin()), 18, "1000044", "ActivationMethod",7);
								}
								writeCVS_(os, precursor.getActivationEnergy(), "1000045", "CollisionEnergy",7);
								os << "\t\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000046\" name=\"EnergyUnit\" value=\"eV\"/>\n";
							}
							os << "\t\t\t\t\t\t</activation>\n";
							os << "\t\t\t\t\t</precursor>\n";
						}
						os << "\t\t\t\t</precursorList>\n";
					}
					os << "\t\t\t</spectrumDesc>\n";
	
					// write the supplementary data?
					if (options_.getWriteSupplementalData())
					{
						//write meta data array descriptions
						for (Size i=0; i<spec.getFloatDataArrays().size(); ++i)
						{
							const MetaInfoDescription& desc = spec.getFloatDataArrays()[i];
							os << "\t\t\t<supDesc supDataArrayRef=\"" << (i+1) << "\">\n";
							if (!desc.isMetaEmpty())
							{
								os << "\t\t\t\t<supDataDesc>\n";
								writeUserParam_(os, desc, 5);
								os << "\t\t\t\t</supDataDesc>\n";
							}
							os << "\t\t\t</supDesc>\n";
						}
					}
					
					//write m/z and intensity arrays
					data_to_encode_.clear();
					for (Size i=0; i<spec.size(); i++)
					{
						data_to_encode_.push_back(spec[i].getPosition()[0]);
					}
					
					writeBinary_(os,spec.size(),"mzArrayBinary");
	
					// intensity
					data_to_encode_.clear();
					for (Size i=0; i<spec.size(); i++)
					{
						data_to_encode_.push_back(spec[i].getIntensity());
					}
					
					writeBinary_(os,spec.size(),"intenArrayBinary");
	
					// write the supplementary data?
					if (options_.getWriteSupplementalData())
					{
						//write supplemental data arrays
						for (Size i=0; i<spec.getFloatDataArrays().size(); ++i)
						{
							const typename MapType::SpectrumType::FloatDataArray& mda = spec.getFloatDataArrays()[i];
							//check if spectrum and meta data array have the same length
							if (mda.size()!=spec.size())
							{
								error(LOAD, String("Length of meta data array (index:'")+i+"' name:'"+mda.getName()+"') differs from spectrum length. meta data array: " + mda.size() + " / spectrum: " + spec.size() +" .");
							}
							//encode meta data array
							data_to_encode_.clear();
							for (Size j=0; j<mda.size(); j++)
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
					if (value=="Zoom")
					{
						spec_.getInstrumentSettings().setZoomScan(true);
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
					}
					else if (value=="MassScan")
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
					}
					else if (value=="SelectedIonDetection")
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::SIM);
					}
					else if (value=="SelectedReactionMonitoring")
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
					}
					else if (value=="ConsecutiveReactionMonitoring")
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CRM);
					}
					else if (value=="ConstantNeutralGainScan")
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CNG);
					}
					else if (value=="ConstantNeutralLossScan")
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::CNL);
					}
					else if (value=="ProductIonScan")
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
						spec_.setMSLevel(2);
					}
					else if (value=="PrecursorIonScan")
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::PRECURSOR);
					}
					else if (value=="EnhancedResolutionScan")
					{
						spec_.getInstrumentSettings().setZoomScan(true);
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
					}
					else
					{
						spec_.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
						warning(LOAD, String("Unknown scan mode '") + value + "'. Assuming full scan");
					}
				}
				else if (accession=="PSI:1000038") //Time in minutes
				{
					spec_.setRT(asDouble_(value)*60); //Minutes to seconds
					if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRT())))
					{
						skip_spectrum_=true;
					}
				}
				else if (accession=="PSI:1000039") //Time in seconds
				{
					spec_.setRT(asDouble_(value));
					if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRT())))
					{
						skip_spectrum_=true;
					}
				}
				else if (accession=="PSI:1000037") //Polarity
				{
					if (value=="Positive" || value=="positive" || value=="+") //be flexible here, actually only the first one is correct
					{
						spec_.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
					}
					else if (value=="Negative" || value=="negative" || value=="-") //be flexible here, actually only the first one is correct
					{
						spec_.getInstrumentSettings().setPolarity(IonSource::NEGATIVE);
					}
					else
					{
						warning(LOAD, String("Invalid scan polarity (PSI:1000037) detected: \"") + value + "\". Valid are 'Positive' or 'Negative'.");
					}
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
					spec_.getPrecursors().back().setMZ(asDouble_(value));			
				}
				else if (accession=="PSI:1000041") //Charge
				{
					if (spec_.getPrecursors().back().getCharge() != 0)
					{
						warning(LOAD, String("Multiple precursor charges detected, expected only one! Ignoring this charge settings! accession=\"") + accession + "\", value=\"" + value + "\"");
						spec_.getPrecursors().back().setCharge(0);
					}
					else
					{
						spec_.getPrecursors().back().setCharge(asInt_(value));
					}
				}
				else if (accession=="PSI:1000042") //Intensity
				{
					spec_.getPrecursors().back().setIntensity(asDouble_(value));		
				}
				else if (accession=="PSI:1000043") //Intensity unit
				{
					//ignored
				}
				else
				{
					error = "PrecursorList.Precursor.IonSelection.UserParam";
				}
			}
			else if (parent_tag=="activation") 
			{
				if (accession=="PSI:1000044") //activationmethod
				{
					spec_.getPrecursors().back().getActivationMethods().insert((Precursor::ActivationMethod)cvStringToEnum_(18, value,"activation method"));
				}
				else if (accession=="PSI:1000045") //Energy
				{
					spec_.getPrecursors().back().setActivationEnergy( asDouble_(value) );
				}
				else if (accession=="PSI:1000046") //Energy unit
				{
					//ignored - we assume electronvolt
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
			else if (parent_tag=="detector")
			{
				if (accession=="PSI:1000026")
				{
					exp_->getInstrument().getIonDetectors().back().setType( (IonDetector::Type)cvStringToEnum_(13,value, "detector type") );
				}
				else if (accession=="PSI:1000028")
				{
					exp_->getInstrument().getIonDetectors().back().setResolution( asDouble_(value) );
				}
				else if (accession=="PSI:1000029")
				{
					exp_->getInstrument().getIonDetectors().back().setADCSamplingFrequency( asDouble_(value) );
				}
				else if (accession=="PSI:1000027")
				{
					exp_->getInstrument().getIonDetectors().back().setAcquisitionMode((IonDetector::AcquisitionMode)cvStringToEnum_(9, value, "acquisition mode") );
				}
				else
				{
					error = "Description.Instrument.Detector.UserParam";
				}
			}
			else if (parent_tag=="source")
			{
				if (accession=="PSI:1000008")
				{
					exp_->getInstrument().getIonSources().back().setIonizationMethod( (IonSource::IonizationMethod)cvStringToEnum_(10, value, "ion source") );
				}
				else if (accession=="PSI:1000007")
				{
					exp_->getInstrument().getIonSources().back().setInletType( (IonSource::InletType)cvStringToEnum_(11, value,"inlet type") );
				}
				else if (accession=="PSI:1000009")
				{
					exp_->getInstrument().getIonSources().back().setPolarity( (IonSource::Polarity)cvStringToEnum_(1, value,"polarity") );
				}
				else 
				{
					error = "Description.Instrument.Source.UserParam";
				}
			}
			else if (parent_tag=="sampleDescription")
			{
				if (accession=="PSI:1000001")
				{
					exp_->getSample().setNumber( value );
				}
				else if (accession=="PSI:1000003")
				{
					exp_->getSample().setState( (Sample::SampleState)cvStringToEnum_(0, value, "sample state") );
				}
				else if (accession=="PSI:1000004")
				{
					exp_->getSample().setMass( asDouble_(value) );
				}
				else if (accession=="PSI:1000005")
				{
					exp_->getSample().setVolume( asDouble_(value) );
				}
				else if (accession=="PSI:1000006")
				{
					exp_->getSample().setConcentration( asDouble_(value) );
				}
				else 
				{
					error = "Description.Admin.SampleDescription.UserParam";
				}
			}
			else if (parent_tag=="analyzer")
			{
				if (accession=="PSI:1000010")
				{
					exp_->getInstrument().getMassAnalyzers().back().setType( (MassAnalyzer::AnalyzerType)cvStringToEnum_(14, value,"analyzer type"));
				}
				else if (accession=="PSI:1000011")
				{
					exp_->getInstrument().getMassAnalyzers().back().setResolution( asDouble_(value) );
				}
				else if (accession=="PSI:1000012")
				{
					exp_->getInstrument().getMassAnalyzers().back().setResolutionMethod( (MassAnalyzer::ResolutionMethod)cvStringToEnum_(2, value,"resolution method"));
				}
				else if (accession=="PSI:1000013")
				{
					exp_->getInstrument().getMassAnalyzers().back().setResolutionType( (MassAnalyzer::ResolutionType)cvStringToEnum_(3, value, "resolution type"));
				}
				else if (accession=="PSI:1000014")
				{
					exp_->getInstrument().getMassAnalyzers().back().setAccuracy( asDouble_(value) );
				}
				else if (accession=="PSI:1000015")
				{
					exp_->getInstrument().getMassAnalyzers().back().setScanRate( asDouble_(value) );
				}
				else if (accession=="PSI:1000016")
				{
					exp_->getInstrument().getMassAnalyzers().back().setScanTime( asDouble_(value) );
				}
				else if (accession=="PSI:1000018")
				{
					exp_->getInstrument().getMassAnalyzers().back().setScanDirection( (MassAnalyzer::ScanDirection)cvStringToEnum_(5, value, "scan direction"));
				}
				else if (accession=="PSI:1000019")
				{
					exp_->getInstrument().getMassAnalyzers().back().setScanLaw( (MassAnalyzer::ScanLaw)cvStringToEnum_(6, value, "scan law"));
				}
				else if (accession=="PSI:1000020")
				{
					// ignored
				}
				else if (accession=="PSI:1000021")
				{
					exp_->getInstrument().getMassAnalyzers().back().setReflectronState( (MassAnalyzer::ReflectronState)cvStringToEnum_(8, value, "reflectron state"));
				}
				else if (accession=="PSI:1000022")
				{
					exp_->getInstrument().getMassAnalyzers().back().setTOFTotalPathLength( asDouble_(value) );
				}
				else if (accession=="PSI:1000023")
				{
					exp_->getInstrument().getMassAnalyzers().back().setIsolationWidth( asDouble_(value) );
				}
				else if (accession=="PSI:1000024")
				{
					exp_->getInstrument().getMassAnalyzers().back().setFinalMSExponent( asInt_(value) );
				}
				else if (accession=="PSI:1000025")
				{
					exp_->getInstrument().getMassAnalyzers().back().setMagneticFieldStrength( asDouble_(value) );
				}
				else if (accession=="PSI:1000017")
				{
					//ignored
				}
				else 
				{
					error = "AnalyzerList.Analyzer.UserParam";
				}
			}
			else if (parent_tag=="additional")
			{
				if (accession=="PSI:1000030")
				{
					exp_->getInstrument().setVendor(value);
				}
				else if (accession=="PSI:1000031")
				{
					exp_->getInstrument().setModel(value);
				}
				else if (accession=="PSI:1000032")
				{
					exp_->getInstrument().setCustomizations(value);
				}
				else 
				{
					error = "Description.Instrument.Additional";
				}
			}
			else if (parent_tag=="processingMethod")
			{
				if (accession=="PSI:1000033")
				{
					data_processing_.getProcessingActions().insert(DataProcessing::DEISOTOPING);
				}
				else if (accession=="PSI:1000034")
				{
					data_processing_.getProcessingActions().insert(DataProcessing::CHARGE_DECONVOLUTION);
				}
				else if (accession=="PSI:1000127")
				{
					data_processing_.getProcessingActions().insert(DataProcessing::PEAK_PICKING);
				}
				else if (accession=="PSI:1000035")
				{
					//ignored
				}
				else 
				{
					error = "DataProcessing.DataProcessing.UserParam";
				}
			}
			else
			{
				warning(LOAD, String("Unexpected cvParam: accession=\"") + accession + "\" value=\"" + value + "\" in tag " + parent_tag);
			}

			if (error != "")
			{
				warning(LOAD, String("Invalid cvParam: accession=\"") + accession + "\" value=\"" + value + "\" in " + error);
			}
			//std::cout << "End of MzDataHander::cvParam_" << std::endl;
		}
		
	} // namespace Internal

} // namespace OpenMS

#endif
