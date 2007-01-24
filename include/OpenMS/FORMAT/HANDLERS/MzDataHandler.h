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

#ifndef OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataExpSettHandler.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DPeak.h>

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
			
			@todo Softly abort parsing after metadata, if only metadata should be read (Thomas S.)
		*/
		template <typename MapType>
		class MzDataHandler
			: public SchemaHandler
		{
		 public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzDataHandler(MapType& exp, const String& filename)
				: SchemaHandler(TAG_NUM,MAP_NUM,filename), // number of tags, number of maps
					exp_(&exp),
					cexp_(0),
					peak_count_(0),
					prec_(0),	acq_(0),
					meta_id_(),
					exp_sett_(),
					decoder_(2),
					spec_write_counter_(1)
	  	{
				fillMaps_(Schemes::MzData[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}

      /// Constructor for a read-only handler
      MzDataHandler(const MapType& exp, const String& filename)
				: SchemaHandler(TAG_NUM,MAP_NUM,filename), // number of tags, number of maps
					exp_(0),
					cexp_(&exp),
					peak_count_(0),
					prec_(0),	acq_(0),
					meta_id_(),
					exp_sett_(),
					decoder_(2),
					spec_write_counter_(1)
  		{
				fillMaps_(Schemes::MzData[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}

      /// Destructor
      virtual ~MzDataHandler(){}
      //@}


			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, const unsigned int length);

  		/// Writes the contents to a stream
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
			enum Tags { TAGNULL, MZDATA, DESCRIPTION, SPECTRUMLIST, SPECTRUM,
								SPECTRUMDESC, SPECTRUMSETTINGS,
								ACQSPEC, ACQUISITION,	SPECTRUMINSTRUMENT, PRECURSORLIST,
								IONSELECTION, ACTIVATION, PRECURSOR, SUPDATADESC,
								SUPDESC, SUPSRCFILE, DATA, INTENARRAYBINARY, MZARRAYBINARY,
								CVPARAM, USERPARAM, ACQINSTRUMENT, ACQSETTINGS, ACQDESC, CVLOOKUP,
								SUPARRAYBINARY, SUPARRAY, ARRAYNAME, COMMENTS,
								NAMEOFFILE,	PATHTOFILE, FILETYPE, TAG_NUM};


			/** @brief indices for attributes used by mzData

			If you add attributes, also add them to XMLSchemes.h
			*/
			enum Attributes { ATTNULL, NAME, VALUE, ID, COUNT, SPECTRUMTYPE, METHOD_OF_COMBINATION,
			                 ACQNUMBER, MSLEVEL, MZRANGE_START, MZRANGE_STOP,
			                 SUP_DATA_ARRAY_REF, ATT_PRECISION, ATT_ENDIAN, LENGTH, VERSION, ATT_NUM};
			
			/** @brief indices for ontology terms used by mzData

			If you add terms, also add them to XMLSchemes.h
			*/
			enum Ontology { ONTNULL, SCANMODE, POLARITY, TIMEMIN, TIMESEC,
										MZ_ONT, CHARGESTATE, INTENSITY, IUNITS, METHOD, ENERGY, EUNITS};

			/** @brief indices for enum2str-maps used by mzData

			Used to access enum2str_().
			If you add maps, also add them to XMLSchemes.h.
			Add no elements to the enum after MAP_NUM.
			Each map corresponds to a string in XMLSchemes.h.
			*/
			enum MapTypes {	PRECISION, ENDIAN, EUNITSMAP,
				SCANMODEMAP, POLARITYMAP, ACTMETHODMAP, ONTOLOGYMAP, TAGMAP, ATTMAP, MAP_NUM};

			/// Possible precisions for Base64 data encoding
			enum Precision { UNKNOWN_PRECISION, REAL, DOUBLE};

			/// Possible endian-types for Base64 data encoding
			enum Endian { UNKNOWN_ENDIAN, LITTLE, BIG};

			typedef typename MapType::SpectrumType SpectrumType;
			typedef typename MapType::PeakType PeakType;

			PeakFileOptions options_;
		
			/**@name temporary datastructures to hold parsed data */
			//@{
				
			Size peak_count_;
			SpectrumType spec_;
			Precursor* prec_;
			Acquisition* acq_;
			String meta_id_;
			std::vector<String> data_;
			std::vector<String> array_name_;
			std::vector<Precision> precisions_;
			std::vector<Endian> endians_;
			//@}

			/// stream to collect experimental settings
			std::stringstream exp_sett_;

			/// Decoder/Encoder for Base64-data in MzData
			std::vector<Base64> decoder_;

			/// spectrum counter (spectra without peaks are not written)
			UnsignedInt spec_write_counter_;

			/// fills the experiment with peaks
			void fillData_();

			/** @brief read attributes of MzData's cvParamType

			Example:
			&lt;cvParam cvLabel="psi" accession="PSI:1000001" name="@p name" value="@p value"/&gt;
			@p name and sometimes @p value are defined in the MzData ontology.
			*/
			void cvParam_(const String& name, const String& value);

			/** @brief read attributes of MzData's userParamType

			Example:
			&lt;userParam name="@p name" value="@p value"/&gt;
			@p name and @p value are stored as MetaValues
			*/
			void userParam_(const String& name, const String& value);

			/// write binary data to stream using the first decoder_ (previously filled)
			inline void writeBinary_(std::ostream& os, Size size, const String& tag, const String& desc="", int id=-1)
			{
				os 	<< "\t\t\t<" << tag;
				if (id>=0)
				{
					os << " id=\"" << id << "\"";
				}
				os << ">\n";
				if (desc!="")
				{
					os << "\t\t\t\t<arrayName>" << desc << "</arrayName>\n";
				}
				os << "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\""
					 << size << "\">"
					 << decoder_[0].encodeFloat()
					 << "</data>\n\t\t\t</" << tag << ">\n";
			}

			inline double getDatum_(const std::vector<void*>& ptrs,  UnsignedInt member, UnsignedInt index)
			{
				if (precisions_[member]==DOUBLE)
				{
					return static_cast<double*>(ptrs[member])[index];
				}
				else
				{
					return static_cast<float*>(ptrs[member])[index];
				}
			}

			/**
				 @brief Write supplemental data for derived classes of DPeak, e.g. for
				 picked peaks.  Default is to do nothing.
			*/
			template <typename ContainerType>
			void writeDerivedPeakSupplementalData_( std::ostream& /* os */ , ContainerType const & /* container */ )
			{
			}

			/**
				 @brief Read supplemental data for derived classes of DPeak, e.g. for
				 picked peaks.  Default is to do nothing.
			*/
			template <typename PeakType>
			void readPeakSupplementalData_( std::vector<void*>& /*data*/, PeakType& /*peak*/, Size /*n*/)
			{
			}

		};

		//--------------------------------------------------------------------------------

		template <typename MapType>
		void MzDataHandler<MapType>::characters(const XMLCh* const chars, const unsigned int /*length*/)
		{
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << xercesc::XMLString::transcode(chars);
				return;
			}

			// find the tag that the parser is in right now
 			for (Size i=0; i<is_parser_in_tag_.size(); i++)
 			{
				if (is_parser_in_tag_[i])
				{
					switch(i) 
					{
	  				case COMMENTS:		// <comment> is child of more than one other tags
							if (is_parser_in_tag_[ACQDESC])
							{
								spec_.setComment( xercesc::XMLString::transcode(chars) );
							}
							else
							{
								warning(String("Unhandled tag \"comments\" with content:") + xercesc::XMLString::transcode(chars));
							}
							break;
						case DATA:
							if (options_.getMetadataOnly()) break;
							
							data_.push_back(xercesc::XMLString::transcode(chars));		// store characters for later
							if (is_parser_in_tag_[MZARRAYBINARY]) array_name_.push_back("mz");
							if (is_parser_in_tag_[INTENARRAYBINARY]) array_name_.push_back("intens");
							break;
					  case ARRAYNAME:
							array_name_.push_back(xercesc::XMLString::transcode(chars));
							if (spec_.getMetaInfoDescriptions().find(meta_id_) != spec_.getMetaInfoDescriptions().end())
							{
								spec_.getMetaInfoDescriptions()[meta_id_].setName(xercesc::XMLString::transcode(chars));
							}
							break;
						case NAMEOFFILE: 	// <nameOfFile> is child of more than one other tags
							if (is_parser_in_tag_[SUPSRCFILE])
							{
								spec_.getMetaInfoDescriptions()[meta_id_].getSourceFile().setNameOfFile( xercesc::XMLString::transcode(chars) );
							}
							else
							{
								warning(String("Unhandled tag \"nameOfFile\" with content: ") + xercesc::XMLString::transcode(chars));
							}
							break;
						case PATHTOFILE: // <pathOfFile> is child of more than one other tags
							if (is_parser_in_tag_[SUPSRCFILE])
							{
								spec_.getMetaInfoDescriptions()[meta_id_].getSourceFile().setPathToFile( xercesc::XMLString::transcode(chars) );
							}
							else
							{
								warning(String("Unhandled tag \"pathToFile\" with content: ") + xercesc::XMLString::transcode(chars));
							}
							break;
						case FILETYPE: // <fileType> is child of more than one other tags
							if (is_parser_in_tag_[SUPSRCFILE])
							{
								spec_.getMetaInfoDescriptions()[meta_id_].getSourceFile().setFileType( xercesc::XMLString::transcode(chars) );
							}
							else
							{
								warning(String("Unhandled tag \"fileType\" with content: ") + xercesc::XMLString::transcode(chars));
							}
							break;	
					}
				}
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
		{
			
// 			std::cout << "Start: '" << xercesc::XMLString::transcode(qname) << "'" << std::endl;
			
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << '<' << xercesc::XMLString::transcode(qname);
				Size n=attributes.getLength();
				for (Size i=0; i<n; ++i)
				{
					exp_sett_ << ' ' << xercesc::XMLString::transcode(attributes.getQName(i)) << "=\""	<< xercesc::XMLString::transcode(attributes.getValue(i)) << '\"';
				}
				exp_sett_ << '>';
				return;
			}

			int tag = enterTag(qname, attributes);

			// Do something depending on the tag
			String tmp_type;
			switch(tag) 
			{
				case DESCRIPTION: 
					exp_sett_ << '<' << xercesc::XMLString::transcode(qname) << '>'; 
					break;
				case CVPARAM: 
				case USERPARAM:
				{
					String name = getAttributeAsString_(NAME);
					String value = getAttributeAsString_(VALUE);
					
					if (name == "")
					{
						error("missing required attribute 'name'");
					}
					else if  (value == "")
					{
						error("missing required attribute 'value'");
					}
					else 
					{
						if (tag == USERPARAM)
						{
							userParam_(name, value);
						}
						else
						{
							cvParam_(name, value);
						}
					}
					break;
			  }
				case SUPARRAYBINARY:
					meta_id_ = getAttributeAsString_(ID);
					break;
				case SPECTRUM:
// 					exp_->push_back(SpectrumType());
					spec_ = SpectrumType();
					break;
			  case SPECTRUMLIST:
			  	//std::cout << Date::now() << " Reserving space for spectra" << std::endl;
			  	exp_->reserve( asSignedInt_(getAttributeAsString_(COUNT)) );
			  	//std::cout << Date::now() << " done" << std::endl;
			  	break;
				case ACQSPEC:
					tmp_type = getAttributeAsString_(SPECTRUMTYPE);
					if  (tmp_type == "CentroidMassSpectrum")
					{
						spec_.setType(SpectrumSettings::PEAKS);
					}
					else if (tmp_type == "ContinuumMassSpectrum")
					{
						spec_.setType(SpectrumSettings::RAWDATA);
					}
					else
					{
						spec_.setType(SpectrumSettings::UNKNOWN);
					}
					
					spec_.getAcquisitionInfo().setMethodOfCombination(getAttributeAsString_(METHOD_OF_COMBINATION));
					break;
				case ACQUISITION:
					{
						spec_.getAcquisitionInfo().insert(spec_.getAcquisitionInfo().end(), Acquisition());
						acq_ = &(spec_.getAcquisitionInfo().back());
						acq_->setNumber(asSignedInt_(getAttributeAsString_(ACQNUMBER)));
					}	
					break;
				case SPECTRUMINSTRUMENT:
				case ACQINSTRUMENT:
				{
					spec_.setMSLevel(asSignedInt_(getAttributeAsString_(MSLEVEL)));
					String start = getAttributeAsString_(MZRANGE_START);
					String stop = getAttributeAsString_(MZRANGE_STOP);
					
					if  (start != "")
					{
						spec_.getInstrumentSettings().setMzRangeStart(asDouble_(start));
					}
					if  (stop != "")
					{
						spec_.getInstrumentSettings().setMzRangeStop(asDouble_(stop));
					}
					
					if (options_.hasMSLevels())
					{
						if (!options_.containsMSLevel(spec_.getMSLevel()))
						{
							// HACK: skip the top 4 tags: spectrum, spectrumDesc, SpectrumSettings and spectrumInstrument
							if (is_parser_in_tag_[SPECTRUM] && is_parser_in_tag_[SPECTRUMDESC] && is_parser_in_tag_[SPECTRUMSETTINGS])
							{
								for (int i = 0; i != 4; ++i) skip_tag_.pop();
								for (int i = 0; i != 4; ++i) skip_tag_.push(true);
							}
						}
					}
					
					break;
				}
				case PRECURSOR:
					prec_ = &(spec_.getPrecursor());
					//UNHANDLED: "spectrumRef";
					break;
				case SUPDESC:
					meta_id_ = getAttributeAsString_(SUP_DATA_ARRAY_REF);
					break;
				case DATA:
					if (options_.getMetadataOnly()) break;
					
					// store precision for later
					precisions_.push_back((Precision)str2enum_(PRECISION, getAttributeAsString_(ATT_PRECISION)));
					endians_.push_back((Endian)str2enum_(ENDIAN, getAttributeAsString_(ATT_ENDIAN)));
					if (is_parser_in_tag_[MZARRAYBINARY])
					{
						peak_count_ = asSignedInt_(getAttributeAsString_(LENGTH));
						//std::cout << Date::now() << " Reserving space for peaks" << std::endl;
						spec_.getContainer().reserve(peak_count_);
					}
					break;
				case MZDATA:
					{
						String s = getAttributeAsString_(VERSION);
						for (UnsignedInt index=0; index<Schemes::MzData_num; ++index)
						{
							if (s!=String(schema_) && s.hasSubstring(Schemes::MzData[index][0]))
							{
								schema_ = index;
								// refill maps with older schema
								for (Size i=0; i<str2enum_array_.size(); i++)	str2enum_array_[i].clear();
								for (Size i=0; i<enum2str_array_.size(); i++)	enum2str_array_[i].clear();
								fillMaps_(Schemes::MzData[schema_]);
							}
						}
					}
					break;
			}
		}


		template <typename MapType>
		void MzDataHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			
// 			std::cout << "End: '" << xercesc::XMLString::transcode(qname) << "'" << std::endl;
			
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << "</" << xercesc::XMLString::transcode(qname) << ">\n";
				if (xercesc::XMLString::transcode(qname) != enum2str_(TAGMAP,DESCRIPTION))	
				{
					return;
				}
			}
			
			bool skip = skip_tag_.top();

			int tag = leaveTag(qname);

			// Do something depending on the tag
			switch(tag) {
			case DESCRIPTION:
				// delegate control to ExperimentalSettings handler
				{
					// initialize parser
					xercesc::XMLPlatformUtils::Initialize();
					xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
					parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
					parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);

					MzDataExpSettHandler handler( exp_->getExperimentalSettings(),file_); // *((ExperimentalSettings*)
					handler.resetErrors();
					parser->setContentHandler(&handler);
					parser->setErrorHandler(&handler);
					
					String tmp(exp_sett_.str().c_str());
// 					std::cout << tmp << std::endl;
					xercesc::MemBufInputSource source((const XMLByte*)(tmp.c_str()), tmp.size(), "dummy");
	      	parser->parse(source);
	      	delete(parser);
				}
				break;
			case SPECTRUM:
				if (!skip)
				{
					fillData_();
					exp_->push_back(spec_);
				}
				
				data_.clear();
				array_name_.clear();
				precisions_.clear();
				endians_.clear();
				break;
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::userParam_(const String& name, const String& value)
		{
			if(is_parser_in_tag_[SPECTRUMINSTRUMENT] || is_parser_in_tag_[ACQINSTRUMENT])
			{
				setAddInfo_(spec_.getInstrumentSettings(), name, value, "SpectrumSettings.SpectrumInstrument.UserParam");
			}
			else if(is_parser_in_tag_[ACQUISITION])
			{
				setAddInfo_(*acq_, name, value, "SpectrumSettings.AcqSpecification.Acquisition.UserParam");
			}
			else if (is_parser_in_tag_[IONSELECTION])
			{
				setAddInfo_(spec_.getPrecursorPeak(), name, value, "PrecursorList.Precursor.IonSelection.UserParam");
			}
			else if (is_parser_in_tag_[ACTIVATION])
			{
				setAddInfo_(*prec_, name, value, "PrecursorList.Precursor.Activation.UserParam");
			}
			else if (is_parser_in_tag_[SUPDATADESC])
			{
				setAddInfo_(spec_.getMetaInfoDescriptions()[meta_id_], name, value, "Spectrum.SupDesc.SupDataDesc.UserParam");
			}
			else
			{
				warning("Invalid userParam: name=\"" + name + ", value=\"" + value + "\"");
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::cvParam_(const String& name, const String& value)
		{
			int ont = str2enum_(ONTOLOGYMAP, name, "cvParam element"); // index of current ontology term
			std::string error = "";
			if(is_parser_in_tag_[SPECTRUMINSTRUMENT] || is_parser_in_tag_[ACQINSTRUMENT])
			{
				InstrumentSettings& sett = spec_.getInstrumentSettings();
				bool skip = false;
				
				switch (ont)
				{
					case SCANMODE:
						sett.setScanMode( (InstrumentSettings::ScanMode)str2enum_(SCANMODEMAP, value) );
						break;
					case TIMEMIN:
						spec_.setRetentionTime(asFloat_(value)*60); //Minutes to seconds
						if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRetentionTime())))
						{
							skip = true;
						}
						break;
					case TIMESEC:
						spec_.setRetentionTime(asFloat_(value));
						if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRetentionTime())))
						{
							skip = true;
						}
						break;
					case POLARITY:
						sett.setPolarity( (IonSource::Polarity)str2enum_(POLARITYMAP, value) );
						break;
				  default:      
				  	error = "SpectrumDescription.SpectrumSettings.SpectrumInstrument";
				}
				
				if (skip)
				{
					// HACK: skip the top five tags: spectrum, spectrumDesc, spectrumSettings, {spectrum,acq}Instrument and cvParam
					if (is_parser_in_tag_[SPECTRUM] && is_parser_in_tag_[SPECTRUMDESC] && is_parser_in_tag_[SPECTRUMSETTINGS])
					{
						for (int i = 0; i != 5; ++i) skip_tag_.pop();
						for (int i = 0; i != 5; ++i) skip_tag_.push(true);
						
						return;
					}
				}
			}
			else if (is_parser_in_tag_[IONSELECTION]) 
			{
				switch (ont)
				{
					case MZ_ONT:
						spec_.getPrecursorPeak().getPosition()[0] = asFloat_(value);
						break;
					case CHARGESTATE:
						spec_.getPrecursorPeak().setCharge(asSignedInt_(value));
						break;
					case INTENSITY:
						spec_.getPrecursorPeak().getIntensity() = asFloat_(value);
						break;
					case IUNITS:
						setAddInfo_(spec_.getPrecursorPeak(),"#IntensityUnits", value, "Precursor.IonSelection.IntensityUnits");
						break;
					default:
						error = "PrecursorList.Precursor.IonSelection.UserParam";
				}
			}
			else if (is_parser_in_tag_[ACTIVATION]) 
			{
				switch (ont)
				{
					case METHOD:
						prec_->setActivationMethod((Precursor::ActivationMethod)str2enum_(ACTMETHODMAP, value));
						break;
					case ENERGY: 
						prec_->setActivationEnergy( asFloat_(value) );
						break;
					case EUNITS:
						prec_->setActivationEnergyUnit((Precursor::EnergyUnits)str2enum_(EUNITSMAP, value));
						break;
					default:
						error = "PrecursorList.Precursor.Activation.UserParam";
				}
			}
			else
			{
				warning(String("Invalid cvParam: name=\"") + name + ", value=\"" + value + "\"");
			}

			if (error != "")
			{
				warning(String("Invalid cvParam: name=\"") + name + ", value=\"" + value + "\" in " + error);
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::fillData_()
		{
			//std::cout << "reading scan" << std::endl;
			if (data_.size() > decoder_.size()) // not enough decoder
				decoder_.resize(data_.size());

			std::vector<void*> ptrs(data_.size(),0);		//pointers to data of each decoder
			for (Size i=0; i<data_.size(); i++)
			{
				if (precisions_[i]==DOUBLE)	// precision 64 Bit
				{
					if (endians_[i]==BIG)
					{
						ptrs[i] = static_cast<void*>( decoder_[i].decodeDoubleCorrected(data_[i].c_str(), data_[i].length()));
					}
					else
					{
						ptrs[i] = static_cast<void*>( decoder_[i].decodeDouble(data_[i].c_str(), data_[i].length()));
					}
				}
				else
				{											// precision 32 Bit
					if (endians_[i]==BIG)
					{
						ptrs[i] = static_cast<void*>(decoder_[i].decodeFloatCorrected(data_[i].c_str(), data_[i].length()));
					}
					else
					{
						ptrs[i] = static_cast<void*>(decoder_[i].decodeFloat(data_[i].c_str(), data_[i].length()));
					}
				}
			}

			const int MZ = 0;
			const int INTENS = 1;

			// this works only if MapType::PeakType is at leat DRawDataPoint
			{
				//push_back the peaks into the container
				for (Size n = 0 ; n < peak_count_ ; n++)
				{
					double mz = getDatum_(ptrs,MZ,n);
					double intensity = getDatum_(ptrs,INTENS,n);
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(mz)))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(intensity))))
					{
						spec_.insert(spec_.end(), PeakType());
						spec_.back().getIntensity() = intensity;
						spec_.back().getPosition()[0] = mz;
						//read supplemental data for derived classes (do nothing for DPeak)
						readPeakSupplementalData_(ptrs,spec_.back(),n);
					}
				}
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::writeTo(std::ostream& os)
		{
			os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
				 << "<mzData version=\"1.05\" accessionNumber=\"OpenMS:\">\n";

			// delegate control to ExperimentalSettings handler
			Internal::MzDataExpSettHandler handler( cexp_->getExperimentalSettings(),"");
			handler.writeTo(os);

			//determine how many spectra there are (count only those with peaks)
			UnsignedInt count_tmp_  = 0;
			for (UnsignedInt s=0; s<cexp_->size(); s++)
			{
				const SpectrumType& spec = (*cexp_)[s];
				if (spec.size()!=0) ++count_tmp_;
			}

			os << "\t<spectrumList count=\"" << count_tmp_ << "\">\n";

			int spectrum_ref = -1;
			for (UnsignedInt s=0; s<cexp_->size(); s++)
			{
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
						os << "CentroidMassSpectrum";
					}
					else if (spec.getType()==SpectrumSettings::RAWDATA)
					{
						os << "ContinuumMassSpectrum";
					}

					os << "\" methodOfCombination=\""
						 << spec.getAcquisitionInfo().getMethodOfCombination() << "\" count=\""
						 << spec.getAcquisitionInfo().size() << "\">\n";
					for (Size i=0; i<spec.getAcquisitionInfo().size(); ++i)
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

				writeCVS_(os, spec.getInstrumentSettings().getScanMode(), SCANMODEMAP,
									"1000036", "ScanMode",6);
				writeCVS_(os, spec.getInstrumentSettings().getPolarity(), POLARITYMAP,
									"1000037", "Polarity",6);
				//Retiontion time already in TimeInSeconds
				//writeCVS_(os, spec.getRetentionTime()/60, "1000038", "TimeInMinutes",6);
				writeCVS_(os, spec.getRetentionTime(), "1000039", "TimeInSeconds",6);
				writeUserParam_(os, spec.getInstrumentSettings(), 6);
				os 	<< "\t\t\t\t\t</spectrumInstrument>\n\t\t\t\t</spectrumSettings>\n";

				typedef typename SpectrumType::PrecursorPeakType PrecursorPeak;
				if (spec.getPrecursorPeak() != PrecursorPeak()
						|| spec.getPrecursor() != Precursor())
				{
					os	<< "\t\t\t\t<precursorList count=\"1\">\n"
							<< "\t\t\t\t\t<precursor msLevel=\"2\" spectrumRef=\""
							<< spectrum_ref << "\">\n";
					if (spec.getPrecursorPeak() != PrecursorPeak())
					{
						const PrecursorPeak& peak = spec.getPrecursorPeak();
						os << "\t\t\t\t\t\t<ionSelection>\n";
						writeCVS_(os, peak.getPosition()[0], "1000040", "MassToChargeRatio",7);
						writeCVS_(os, peak.getCharge(), "1000041", "ChargeState",7);
						writeCVS_(os, peak.getIntensity(), "1000042", "Intensity",7);
						if (peak.metaValueExists("#IntensityUnits"))
						{
							writeCVS_(os, String(peak.getMetaValue("#IntensityUnits")),
												"1000043", "IntensityUnits",7);
						}
						writeUserParam_(os, peak, 7);
						os << "\t\t\t\t\t\t</ionSelection>\n";
						os << "\t\t\t\t\t\t<activation>\n";
					}
					if (spec.getPrecursor() != Precursor())
					{
						const Precursor& prec = spec.getPrecursor();
						writeCVS_(os, prec.getActivationMethod(), ACTMETHODMAP, "1000044", "Method",7);
						writeCVS_(os, prec.getActivationEnergy(), "1000045", "CollisionEnergy",7);
						writeCVS_(os, prec.getActivationEnergyUnit(), EUNITSMAP,"1000046", "EnergyUnits",7);
						writeUserParam_(os, prec,7);
					}
						os << "\t\t\t\t\t\t</activation>\n";
					os << "\t\t\t\t\t</precursor>\n"
						 << "\t\t\t\t</precursorList>\n";
				}
				os << "\t\t\t</spectrumDesc>\n";

				typedef const std::map<String,MetaInfoDescription> Map;
				if (spec.getMetaInfoDescriptions().size()>0)
				{
					for (Map::const_iterator it = spec.getMetaInfoDescriptions().begin();
							 it != spec.getMetaInfoDescriptions().end(); ++it)
					{
						os << "\t\t\t<supDesc supDataArrayRef=\"" << it->first << "\">\n";
						if (!it->second.isMetaEmpty())
						{
							os << "\t\t\t\t<supDataDesc>\n";
							writeUserParam_(os, it->second, 5);
							os << "\t\t\t\t</supDataDesc>\n";
						}
						if (it->second.getSourceFile()!=SourceFile())
						{
							os << "\t\t\t\t<supSourceFile>\n"
					 				<< "\t\t\t\t\t<nameOfFile>" << it->second.getSourceFile().getNameOfFile()
									<< "</nameOfFile>\n"
					 				<< "\t\t\t\t\t<pathToFile>" << it->second.getSourceFile().getPathToFile()
									<< "</pathToFile>\n";
							if (it->second.getSourceFile().getFileType()!="")	os << "\t\t\t\t\t<fileType>"
								<< it->second.getSourceFile().getFileType()	<< "</fileType>\n";
							os << "\t\t\t\t</supSourceFile>\n";
						}
						os << "\t\t\t</supDesc>\n";
					}
				}

				// m/z
				float* tmp = decoder_[0].getFloatBuffer(spec.size());
				for (UnsignedInt i=0; i<spec.size(); i++)
					tmp[i] = spec.getContainer()[i].getPosition()[0];
				writeBinary_(os,spec.size(),"mzArrayBinary");

				// intensity
				for (UnsignedInt i=0; i<spec.size(); i++)
					tmp[i] = spec.getContainer()[i].getIntensity();
				writeBinary_(os,spec.size(),"intenArrayBinary");

				// write the supplementary data for picked peaks (is a no-op otherwise)
				if (options_.getWriteSupplementalData())
				{
					this->writeDerivedPeakSupplementalData_(os, spec.getContainer());
				}

				os <<"\t\t</spectrum>\n";
			}
			os << "\t</spectrumList>\n</mzData>\n";
		}


		/**
			 @brief Partial specialization that writes supplemental data for picked peaks.

			 @note Partial specialization must be placed in .C file
		*/
		template <>
		template <>
		void MzDataHandler <MSExperiment<DPickedPeak<1,KernelTraits> > >::writeDerivedPeakSupplementalData_ < DPeakArray<1, DPickedPeak<1,KernelTraits> > >( std::ostream& os, DPeakArray<1, DPickedPeak<1,KernelTraits> > const & container);

		/**
			 @brief Partial specialization that reads supplemental data for picked peaks.

			 @note Partial specialization must be placed in .C file
		*/
		template <>
		template <>
		void MzDataHandler <MSExperiment<DPickedPeak<1,KernelTraits> > >::readPeakSupplementalData_ < DPickedPeak<1,KernelTraits> >( std::vector<void*>& data, DPickedPeak<1,KernelTraits>& peak, Size n);


	} // namespace Internal

} // namespace OpenMS




#endif
