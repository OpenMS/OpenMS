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
#include <OpenMS/KERNEL/PickedPeak1D.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

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
			: public SchemaHandler
		{
		 public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzDataHandler(MapType& exp, const String& filename, ProgressLogger& logger)
				: SchemaHandler(TAG_NUM,MAP_NUM,filename), // number of tags, number of maps
					exp_(&exp),
					cexp_(0),
					peak_count_(0),
					prec_(0),	acq_(0),
					meta_id_(),
					exp_sett_(),
					decoder_(),
					spec_write_counter_(1),
					logger_(logger)
	  	{
				fillMaps_(Schemes::MzData[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}

      /// Constructor for a read-only handler
      MzDataHandler(const MapType& exp, const String& filename, const ProgressLogger& logger)
				: SchemaHandler(TAG_NUM,MAP_NUM,filename), // number of tags, number of maps
					exp_(0),
					cexp_(&exp),
					peak_count_(0),
					prec_(0),	acq_(0),
					meta_id_(),
					exp_sett_(),
					decoder_(),
					spec_write_counter_(1),
					logger_(logger)
  		{
				fillMaps_(Schemes::MzData[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
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
			                 SUP_DATA_ARRAY_REF, ATT_PRECISION, ATT_ENDIAN, LENGTH, VERSION, ACCESSION, ATT_NUM};
			
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
				
			UInt peak_count_;
			SpectrumType spec_;
			Precursor* prec_;
			Acquisition* acq_;
			String meta_id_;
			/// encoded data which is read and has to be decoded
			std::vector<String> data_to_decode_;
			/// floating point numbers which have to be encoded and written
			std::vector<Real> data_to_encode_;
			std::vector<std::vector<Real> > decoded_list_;
			std::vector<std::vector<DoubleReal> > decoded_double_list_;
			std::vector<String> array_name_;
			std::vector<Precision> precisions_;
			std::vector<Endian> endians_;
			//@}

			/// stream to collect experimental settings
			std::stringstream exp_sett_;

			/// Decoder/Encoder for Base64-data in MzData
			Base64 decoder_;

			/// spectrum counter (spectra without peaks are not written)
			UInt spec_write_counter_;

			const ProgressLogger& logger_;

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

			/// write binary data to stream (first one)
			inline void writeBinary_(std::ostream& os, UInt size, const String& tag, const String& desc="", int id=-1)
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

				std::string str;				
				decoder_.encode(data_to_encode_, Base64::LITTLEENDIAN, str);
				data_to_encode_.clear();
				os << "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\""
					 << size << "\">"
					 << str
					 << "</data>\n\t\t\t</" << tag << ">\n";
			}

			inline double getDatum_(UInt member, UInt index)
			{
				if (precisions_[member]==DOUBLE)
				{
//					std::cout << "fetching double-precision index " << index << " of member" << member << std::endl;
					return decoded_double_list_[member][index];
				}
				else
				{
//					std::cout << "fetching single-precision index " << index << " of member" << member << std::endl;
					return (double) decoded_list_[member][index];
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
			void readPeakSupplementalData_( PeakType& /*peak*/, UInt /*n*/)
			{
			}
			
			private:
				MzDataHandler(); // not impelmented -> private
		};

		//--------------------------------------------------------------------------------

		template <typename MapType>
		void MzDataHandler<MapType>::characters(const XMLCh* const chars, unsigned int /*length*/)
		{
			char* transcoded_chars = sm_.convert(chars);
				
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << transcoded_chars;
				return;
			}

			// find the tag that the parser is in right now
 			for (UInt i=0; i<is_parser_in_tag_.size(); i++)
 			{
				if (is_parser_in_tag_[i])
				{
					switch(i) 
					{
	  				case COMMENTS:		// <comment> is child of more than one other tags
							if (is_parser_in_tag_[ACQDESC])
							{
								spec_.setComment( transcoded_chars );
							}
							else
							{
								warning(String("Unhandled tag \"comments\" with content:") + transcoded_chars);
							}
							break;
						case DATA:
							//chars may be split to several chunks => concatenate them
							data_to_decode_.back() += transcoded_chars;
							break;
					  case ARRAYNAME:
							array_name_.push_back(transcoded_chars);
							if (spec_.getMetaInfoDescriptions().find(meta_id_) != spec_.getMetaInfoDescriptions().end())
							{
								spec_.getMetaInfoDescriptions()[meta_id_].setName(transcoded_chars);
							}
							break;
						case NAMEOFFILE: 	// <nameOfFile> is child of more than one other tags
							if (is_parser_in_tag_[SUPSRCFILE])
							{
								spec_.getMetaInfoDescriptions()[meta_id_].getSourceFile().setNameOfFile( transcoded_chars );
							}
							else
							{
								warning(String("Unhandled tag \"nameOfFile\" with content: ") + transcoded_chars);
							}
							break;
						case PATHTOFILE: // <pathOfFile> is child of more than one other tags
							if (is_parser_in_tag_[SUPSRCFILE])
							{
								spec_.getMetaInfoDescriptions()[meta_id_].getSourceFile().setPathToFile( transcoded_chars );
							}
							else
							{
								warning(String("Unhandled tag \"pathToFile\" with content: ") + transcoded_chars);
							}
							break;
						case FILETYPE: // <fileType> is child of more than one other tags
							if (is_parser_in_tag_[SUPSRCFILE])
							{
								spec_.getMetaInfoDescriptions()[meta_id_].getSourceFile().setFileType( transcoded_chars );
							}
							else
							{
								warning(String("Unhandled tag \"fileType\" with content: ") + transcoded_chars);
							}
							break;	
					}
				}
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
		{
			//std::cout << "begin startelement" << std::endl;
			
// 			std::cout << "Start: '" << sm_.convert(qname) << "'" << std::endl;
			
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << '<' << sm_.convert(qname);
				UInt n=attributes.getLength();
				for (UInt i=0; i<n; ++i)
				{
					exp_sett_ << ' ' << sm_.convert(attributes.getQName(i)) << "=\""	<< sm_.convert(attributes.getValue(i)) << '\"';
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
					exp_sett_ << '<' << sm_.convert(qname) << '>'; 
					break;
				case CVPARAM:
				{
					String accession = getAttributeAsString_(ACCESSION, true, qname);
					String value = getAttributeAsString_(VALUE, false, qname);
					cvParam_(accession, value);
				}
				break;
				case USERPARAM:
				{
					String name = getAttributeAsString_(NAME, true, qname);
					String value = getAttributeAsString_(VALUE, false, qname);
					userParam_(name, value);
			  }
				break;
				case SUPARRAYBINARY:
					meta_id_ = getAttributeAsString_(ID, true, qname);
					break;
				case SPECTRUM:
					spec_ = SpectrumType();
					break;
			  case SPECTRUMLIST:
			  	{
				  	if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				  	//std::cout << Date::now() << " Reserving space for spectra" << std::endl;
				  	UInt count = asInt_(getAttributeAsString_(COUNT, true, qname));
				  	exp_->reserve(count);
				  	logger_.startProgress(0,count,"loading mzData file");
				  	//std::cout << Date::now() << " done" << std::endl;
			  	}
			  	break;
				case ACQSPEC:
					tmp_type = getAttributeAsString_(SPECTRUMTYPE, true, qname);
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
					
					spec_.getAcquisitionInfo().setMethodOfCombination(getAttributeAsString_(METHOD_OF_COMBINATION, true, qname));
					break;
				case ACQUISITION:
					{
						spec_.getAcquisitionInfo().insert(spec_.getAcquisitionInfo().end(), Acquisition());
						acq_ = &(spec_.getAcquisitionInfo().back());
						acq_->setNumber(asInt_(getAttributeAsString_(ACQNUMBER, true, qname)));
					}	
					break;
				case SPECTRUMINSTRUMENT:
				case ACQINSTRUMENT:
				{
					spec_.setMSLevel(asInt_(getAttributeAsString_(MSLEVEL, true, qname)));
					String start = getAttributeAsString_(MZRANGE_START, false, qname);
					String stop = getAttributeAsString_(MZRANGE_STOP, false, qname);
					
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
					meta_id_ = getAttributeAsString_(SUP_DATA_ARRAY_REF, true, qname);
					break;
				case DATA:
					// store precision for later
					precisions_.push_back((Precision)str2enum_(PRECISION, getAttributeAsString_(ATT_PRECISION, true, qname)));
					endians_.push_back((Endian)str2enum_(ENDIAN, getAttributeAsString_(ATT_ENDIAN, true, qname)));

					//reserve enough space in spectrum
					if (is_parser_in_tag_[MZARRAYBINARY])
					{
						peak_count_ = asInt_(getAttributeAsString_(LENGTH, true, qname));
						spec_.getContainer().reserve(peak_count_);
					}					
					break;
				case MZARRAYBINARY:
					array_name_.push_back("mz");
					data_to_decode_.resize(data_to_decode_.size()+1);
					break;
				case INTENARRAYBINARY:
					array_name_.push_back("intens");
					data_to_decode_.resize(data_to_decode_.size()+1);
					break;
				case ARRAYNAME:
					// Note: name is set in closing tag as it is CDATA
					data_to_decode_.resize(data_to_decode_.size()+1);
					break;
				case MZDATA:
					{
						String s = getAttributeAsString_(VERSION, true, qname);
						for (UInt index=0; index<Schemes::MzData_num; ++index)
						{
							if (s!=String(schema_) && s.hasSubstring(Schemes::MzData[index][0]))
							{
								schema_ = index;
								// refill maps with older schema
								for (UInt i=0; i<str2enum_array_.size(); i++)	str2enum_array_[i].clear();
								for (UInt i=0; i<enum2str_array_.size(); i++)	enum2str_array_[i].clear();
								fillMaps_(Schemes::MzData[schema_]);
							}
						}
					}
					break;
			}
			//std::cout << "end startelement" << std::endl;
		}


		template <typename MapType>
		void MzDataHandler<MapType>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static UInt scan_count = 0;
			//std::cout << "begin endelement" << std::endl;
			
// 			std::cout << "End: '" << sm_.convert(qname) << "'" << std::endl;
			
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << "</" << sm_.convert(qname) << ">\n";
				if (sm_.convert(qname) != enum2str_(TAGMAP,DESCRIPTION))	
				{
					return;
				}
			}
			
			bool skip = skip_tag_.top();

			int tag = leaveTag(qname);

			// Do something depending on the tag
			switch(tag) 
			{
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
					logger_.setProgress(++scan_count);
					decoded_list_.clear();
					decoded_double_list_.clear();
					data_to_decode_.clear();
					array_name_.clear();
					precisions_.clear();
					endians_.clear();
					break;
				case MZDATA:
					logger_.endProgress();
					scan_count = 0;
					break;
			}
			//std::cout << "end endelement" << std::endl;
			sm_.clear();
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
		void MzDataHandler<MapType>::cvParam_(const String& accession, const String& value)
		{
			//std::cout << "accession is '" << accession << "'." << std::endl;
			int ont = str2enum_(ONTOLOGYMAP, accession, "cvParam element"); // index of current ontology term

			std::string error = "";
			if(is_parser_in_tag_[SPECTRUMINSTRUMENT] || is_parser_in_tag_[ACQINSTRUMENT])
			{
				InstrumentSettings& sett = spec_.getInstrumentSettings();
				bool skip = false;
				
				switch (ont)
				{
					case SCANMODE:
						sett.setScanMode( (InstrumentSettings::ScanMode)str2enum_(SCANMODEMAP, value, (accession + " value").c_str()) );
						break;
					case TIMEMIN:
						spec_.setRT(asFloat_(value)*60); //Minutes to seconds
						if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRT())))
						{
							skip = true;
						}
						break;
					case TIMESEC:
						spec_.setRT(asFloat_(value));
						if (options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(spec_.getRT())))
						{
							skip = true;
						}
						break;
					case POLARITY:
						sett.setPolarity( (IonSource::Polarity)str2enum_(POLARITYMAP, value, (accession + " value").c_str()) );
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
						spec_.getPrecursorPeak().setCharge(asInt_(value));
						break;
					case INTENSITY:
						spec_.getPrecursorPeak().setIntensity(asFloat_(value));
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
						prec_->setActivationMethod((Precursor::ActivationMethod)str2enum_(ACTMETHODMAP, value, (accession + " value").c_str()));
						break;
					case ENERGY: 
						prec_->setActivationEnergy( asFloat_(value) );
						break;
					case EUNITS:
						prec_->setActivationEnergyUnit((Precursor::EnergyUnits)str2enum_(EUNITSMAP, value, (accession + " value").c_str()));
						break;
					default:
						error = "PrecursorList.Precursor.Activation.UserParam";
				}
			}
			else
			{
				warning(String("Invalid cvParam: accession=\"") + accession + ", value=\"" + value + "\"");
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
				if (precisions_[i]==DOUBLE)	// precision 64 Bit
				{
					if (endians_[i]==BIG)
					{
						//std::cout << "nr. " << i << ": decoding as high-precision big endian" << std::endl;
						decoder_.decode(data_to_decode_[i], Base64::BIGENDIAN, decoded_double);
					}
					else
					{
						//std::cout << "nr. " << i << ": decoding as high-precision little endian" << std::endl;
						decoder_.decode(data_to_decode_[i], Base64::LITTLEENDIAN, decoded_double);
					}
					// push_back the decoded double data - and an empty one into
					// the dingle-precision vector, so that we don't mess up the index
					//std::cout << "list size: " << decoded_double.size() << std::endl;
					decoded_double_list_.push_back(decoded_double);
					decoded_list_.push_back(std::vector<float>());
				}
				else
				{											// precision 32 Bit
					if (endians_[i]==BIG)
					{
						//std::cout << "nr. " << i << ": decoding as low-precision big endian" << std::endl;
						decoder_.decode(data_to_decode_[i], Base64::BIGENDIAN, decoded);
					}
					else
					{
						//std::cout << "nr. " << i << ": decoding as low-precision little endian" << std::endl;
						decoder_.decode(data_to_decode_[i], Base64::LITTLEENDIAN, decoded);
					}
					//std::cout << "list size: " << decoded.size() << std::endl;
					decoded_list_.push_back(decoded);
					decoded_double_list_.push_back(std::vector<double>());
				}
			}

			const int MZ = 0;
			const int INTENS = 1;

			// this works only if MapType::PeakType is at leat DRawDataPoint
			{
				//push_back the peaks into the container
				for (UInt n = 0 ; n < peak_count_ ; n++)
				{
					double mz = getDatum_(MZ, n);
					double intensity = getDatum_(INTENS, n);
					if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(mz)))
					 && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(intensity))))
					{
						spec_.insert(spec_.end(), PeakType());
						spec_.back().setIntensity(intensity);
						spec_.back().setPosition(mz);
						//read supplemental data for derived classes (do nothing for DPeak)
						readPeakSupplementalData_(spec_.back(), n);
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
	
					writeCVS2_(os, spec.getInstrumentSettings().getScanMode(), SCANMODEMAP,
										"1000036", "ScanMode",6);
					writeCVS2_(os, spec.getInstrumentSettings().getPolarity(), POLARITYMAP,
										"1000037", "Polarity",6);
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
							writeCVS2_(os, prec.getActivationMethod(), ACTMETHODMAP, "1000044", "Method",7);
							writeCVS_(os, prec.getActivationEnergy(), "1000045", "CollisionEnergy",7);
							writeCVS2_(os, prec.getActivationEnergyUnit(), EUNITSMAP,"1000046", "EnergyUnits",7);
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
	//				float* tmp = decoder_[0].getFloatBuffer(spec.size());
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
	
					// write the supplementary data for picked peaks (is a no-op otherwise)
					if (options_.getWriteSupplementalData())
					{
						this->writeDerivedPeakSupplementalData_(os, spec.getContainer());
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


		/**
			 @brief Partial specialization that writes supplemental data for picked peaks.

			 @note Partial specialization must be placed in .C file
		*/
		template <>
		template <>
		void MzDataHandler <MSExperiment<PickedPeak1D > >::writeDerivedPeakSupplementalData_ < DPeakArray<PickedPeak1D > >( std::ostream& os, DPeakArray<PickedPeak1D > const & container);

		/**
			 @brief Partial specialization that reads supplemental data for picked peaks.

			 @note Partial specialization must be placed in .C file
		*/
		template <>
		template <>
		void MzDataHandler <MSExperiment<PickedPeak1D > >::readPeakSupplementalData_ < PickedPeak1D >( PickedPeak1D& peak, UInt n);


	} // namespace Internal

} // namespace OpenMS




#endif
