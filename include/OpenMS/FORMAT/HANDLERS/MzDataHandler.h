// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: MzDataHandler.h,v 1.19 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataExpSettHandler.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DPeak.h>



namespace OpenMS
{
	namespace Internal
	{

		/**
			@brief XML handler for MzDataFile
			
			MapType has to be a MSExperiment or have the same interface.
			Do not use this class. It is only needed in MzDataFile.
			
			@todo implement and test SupDesc part (Jens)
		*/
		template <typename MapType>
		class MzDataHandler
			: public SchemaHandler
		{
		 public:
      /**@name Constructors and destructor */
      //@{
      ///
      MzDataHandler(MapType& exp)
				: SchemaHandler(TAG_NUM,MAP_NUM), // number of tags, number of maps
					exp_(&exp), 
					cexp_(0),
					peak_count_(0),
					exp_sett_str_(),
					exp_sett_(exp_sett_str_, IO_ReadWrite),
					decoder_(2),
					spec_write_counter_(1)
	  	{
				fillMaps_(Schemes::MzData[schema_]);	// fill maps with current schema
			}

      ///
      MzDataHandler(const MapType& exp)
				: SchemaHandler(TAG_NUM,MAP_NUM), // number of tags, number of maps
					exp_(0), 
					cexp_(&exp),
					peak_count_(0),
					exp_sett_str_(),
					exp_sett_(),
					decoder_(2),
					spec_write_counter_(1)
  		{
				fillMaps_(Schemes::MzData[schema_]);	// fill maps with current schema
			}

      ///
      virtual ~MzDataHandler(){}
      //@}

      /// This function is called for each closing tag in the XML file.
      virtual bool endElement( const QString & uri, const QString & local_name,
															 const QString & qname );

			/// This function is called for each opening XML tag in the file.
      virtual bool startElement(const QString & uri, const QString & local_name,
																const QString & qname, const QXmlAttributes & attributes );

			/// This function is called by the parser for each chunk of
		  /// characters between two tags.
      virtual bool characters( const QString & chars );

  		///Writes the contents to a stream
			void writeTo(std::ostream& os);

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
				SCANMODEMAP, POLARITYMAP, ACTMETHODMAP, ONTOLOGYMAP, TAGMAP, MAP_NUM};

			/// Possible precisions for Base64 data encoding
			enum Precision { UNKNOWN_PRECISION, REAL, DOUBLE};

			/// Possible endian-types for Base64 data encoding
			enum Endian { UNKNOWN_ENDIAN, LITTLE, BIG};

			typedef typename MapType::SpectrumType SpectrumType;
			typedef typename MapType::PeakType PeakType;

			/**@name temporary datastructures to hold parsed data */
			//@{
				
			Size peak_count_;
			SpectrumType* spec_;
			Precursor* prec_;
			Acquisition* acq_;
			MetaInfoDescription* meta_;
			QString meta_id_;
			std::vector<QString> data_;
			std::vector<QString> array_name_;
			std::vector<Precision> precisions_;
			std::vector<Endian> endians_;
			//@}

			/// string with the xml containing the experimental settings
			QString exp_sett_str_;
			/// stream to collect experimental settings
			QTextStream exp_sett_;

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
			void cvParam_(QString name, QString value);

			/** @brief read attributes of MzData's userParamType

			Example:
			&lt;userParam name="@p name" value="@p value"/&gt;
			@p name and @p value are stored as MetaValues
			*/
			void userParam_(QString name, QString value);

			/// write binary data to stream using the first decoder_ (previously filled)
			inline void writeBinary(std::ostream& os, Size size, const QString& tag,
															const QString& desc="", int id=-1)
			{
				os 	<< "\t\t\t<" << tag;
				if (id>=0)
					os << " id=\"" << id << "\"";
				os << ">\n";
				if (desc!=""){
					os << "\t\t\t\t<arrayName>" << desc << "</arrayName>\n";
				}
				os << "\t\t\t\t<data precision=\"32\" endian=\"little\" length=\""
					 << size << "\">"
					 << decoder_[0].encodeFloat()
					 << "</data>\n\t\t\t</" << tag << ">\n";
			}

			inline double getDatum(const std::vector<void*>& ptrs,
														 UnsignedInt member, UnsignedInt index)
			{
				if (precisions_[member]==DOUBLE)
					return static_cast<double*>(ptrs[member])[index];
				else
					return static_cast<float*>(ptrs[member])[index];
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
		bool MzDataHandler<MapType>::characters( const QString & chars )
		{
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << chars;
				return true;
			}

			// find the tag that the parser is in right now
 			for (Size i=0; i<is_parser_in_tag_.size(); i++)
				if (is_parser_in_tag_[i]){
					switch(i) {
					// Do something with the characters depending on the tag
  				case COMMENTS:		// <comment> is child of more than one other tags
						if (is_parser_in_tag_[ACQDESC])
							spec_->setComment( chars.ascii() );
						else if (useWarnings())
							warning(QXmlParseException(
																				 QString("Unhandled tag \"comments\" with content: %1\n").arg(chars))
										 );
						break;
					case DATA:
						data_.push_back(chars);		// store characters for later
						if (is_parser_in_tag_[MZARRAYBINARY]) array_name_.push_back("mz");
						if (is_parser_in_tag_[INTENARRAYBINARY]) array_name_.push_back("intens");
						break;
				  case ARRAYNAME:
						array_name_.push_back(chars);
						break;
					case NAMEOFFILE: 	// <nameOfFile> is child of more than one other tags
						if (is_parser_in_tag_[SUPSRCFILE])
							meta_->getSourceFile().setNameOfFile( chars.ascii() );
						else if (useWarnings())
							warning(QXmlParseException(
																				 QString("Unhandled tag \"nameOfFile\" with content: %1\n").arg(chars))
										 );
						break;
					case PATHTOFILE: // <pathOfFile> is child of more than one other tags
						if (is_parser_in_tag_[SUPSRCFILE])
							meta_->getSourceFile().setPathToFile( chars.ascii() );
						else if (useWarnings())
							warning(QXmlParseException(
																				 QString("Unhandled tag \"pathToFile\" with content: %1\n").arg(chars))
										 );
						break;
					case FILETYPE: // <fileType> is child of more than one other tags
						if (is_parser_in_tag_[SUPSRCFILE])
							meta_->getSourceFile().setFileType( chars.ascii() );
						else if (useWarnings())
							warning(QXmlParseException(
																				 QString("Unhandled tag \"fileType\" with content: %1\n").arg(chars))
										 );
						break;	
					}
				}
			return true;
		}

		template <typename MapType>
		bool MzDataHandler<MapType>::startElement(const QString & /*uri*/,
																							const QString & /*local_name*/,	const QString & qname, const QXmlAttributes & attributes )
		{
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << '<' << qname;
				Size n=attributes.count();
				for (Size i=0; i<n; ++i)
					exp_sett_ << ' ' << attributes.qName(i) << "=\""	<< attributes.value(i) << '\"';
				exp_sett_ << '>';
				return true;
			}

			int tag = str2enum_(TAGMAP,qname,"opening tag");	// index of current tag
			is_parser_in_tag_[tag] = true;

			// Do something depending on the tag
			switch(tag) {
			case DESCRIPTION: exp_sett_ << '<' << qname << '>'; break;
			case CVPARAM:	cvParam_(attributes.value("name"),attributes.value("value")); break;
		  case USERPARAM:	userParam_(attributes.value("name"),attributes.value("value")); break;
			case SPECTRUM: 
				exp_->insert(exp_->end(),SpectrumType());
				//std::cout << "Capacity: " << exp_->capacity() << std::endl;
				spec_ = &(exp_->back()); 
				break;
		  case SPECTRUMLIST: 
		  	//std::cout << Date::now() << " Reserving space for spectra" << std::endl;
		  	exp_->reserve( asSignedInt_(attributes.value("count")) ); 
		  	//std::cout << Date::now() << " done" << std::endl;
		  	break;
			case ACQSPEC:
				spec_->getAcquisitionInfo().setSpectrumType(attributes.value("spectrumType").ascii());
				spec_->getAcquisitionInfo().setMethodOfCombination(
																													 attributes.value("methodOfCombination").ascii()
																													);
				break;
			case ACQUISITION:
				{
					spec_->getAcquisitionInfo().insert(spec_->getAcquisitionInfo().end(), Acquisition());
					acq_ = &(spec_->getAcquisitionInfo().back());
					acq_->setNumber(asSignedInt_(attributes.value("acqNumber")));
				}	
				break;
			case SPECTRUMINSTRUMENT: case ACQINSTRUMENT:
				spec_->setMSLevel(asSignedInt_(attributes.value("msLevel")));
				if  (!attributes.value("mzRangeStart").isEmpty())
					spec_->getInstrumentSettings().setMzRangeStart(
																												 asDouble_(attributes.value("mzRangeStart"))
																												);
				if  (!attributes.value("mzRangeStop").isEmpty())
					spec_->getInstrumentSettings().setMzRangeStop(
																												asDouble_(attributes.value("mzRangeStop"))
																											 );
				break;
			case PRECURSOR:
				prec_ = &(spec_->getPrecursor());
				//UNHANDLED: attributes.value("spectrumRef");
				break;
			case SUPDATADESC:
				meta_ = new MetaInfoDescription();
				meta_id_ = attributes.value("supDataArrayRef");
				break;
			case DATA:
				// store precision for later
				precisions_.push_back((Precision)str2enum_(PRECISION,attributes.value("precision")));
				endians_.push_back((Endian)str2enum_(ENDIAN,attributes.value("endian")));
				if (is_parser_in_tag_[MZARRAYBINARY])
				{
					peak_count_ = asSignedInt_(attributes.value("length"));
					//std::cout << Date::now() << " Reserving space for peaks" << std::endl;
					spec_->getContainer().reserve(peak_count_);
				}
				break;
			case MZDATA:
				{
					QString s = attributes.value("version");
					for (UnsignedInt index=0; index<Schemes::MzData_num; ++index)
					{
						if (s!=schema_ && s.contains(Schemes::MzData[index][0].c_str()))
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
			return no_error_;
		}


		template <typename MapType>
		bool MzDataHandler<MapType>::endElement
		( const QString & /*uri*/, const QString & /*local_name*/, const QString & qname )
		{
			if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
			{
				exp_sett_ << "</" << qname << ">\n";
				if (qname != enum2str_(TAGMAP,DESCRIPTION))	return true;
			}

	  	int tag = str2enum_(TAGMAP,qname,"closing tag");  // index of current tag
			is_parser_in_tag_[tag] = false;

			// Do something depending on the tag
			switch(tag) {
			case DESCRIPTION:
				// delegate control to ExperimentalSettings handler
				{
					QXmlSimpleReader parser;
					srand(static_cast<unsigned>(time(0)));
					parser.setFeature("http://xml.org/sax/features/namespaces",false);
					parser.setFeature("http://xml.org/sax/features/namespace-prefixes", false);

					MzDataExpSettHandler handler( *((ExperimentalSettings*)exp_));
					parser.setContentHandler(&handler);
					parser.setErrorHandler(&handler);
					QXmlInputSource source;
					source.setData(exp_sett_str_);
					parser.parse(source);
				}
				break;
			case SPECTRUM:
				fillData_();
				data_.clear();
				array_name_.clear();
				precisions_.clear();
				endians_.clear();
				break;
			case SUPDESC:
				spec_->getMetaInfoDescriptions()[meta_id_.ascii()] = *meta_;
				delete meta_;
				break;
			}  
			return true;
		}

		template <typename MapType>
		void MzDataHandler<MapType>::userParam_(QString name, QString value)
		{
			if(is_parser_in_tag_[SPECTRUMINSTRUMENT] || is_parser_in_tag_[ACQINSTRUMENT])
				setAddInfo(spec_->getInstrumentSettings(),
									 name,value,"SpectrumSettings.SpectrumInstrument.UserParam");
			else if(is_parser_in_tag_[ACQUISITION])
				setAddInfo(*acq_,name,value,"SpectrumSettings.AcqSpecification.Acquisition.UserParam");
			else if (is_parser_in_tag_[IONSELECTION])
				setAddInfo(spec_->getPrecursorPeak(),name,value,
									 "PrecursorList.Precursor.IonSelection.UserParam");
			else if (is_parser_in_tag_[ACTIVATION])
				setAddInfo(*prec_,name,value,"PrecursorList.Precursor.Activation.UserParam");
			else if (useWarnings())
				warning(QXmlParseException(
																	 QString("Invalid userParam: name=\"%1\", value=\"%2\"\n").arg(name).arg(value)));
		}

		template <typename MapType>
		void MzDataHandler<MapType>::cvParam_(QString name, QString value)
		{
			int ont = str2enum_(ONTOLOGYMAP,name,"cvParam elment"); // index of current ontology term

			std::string error = "";
			if(is_parser_in_tag_[SPECTRUMINSTRUMENT] || is_parser_in_tag_[ACQINSTRUMENT])
			{
				InstrumentSettings& sett = spec_->getInstrumentSettings();
				switch (ont){
				case SCANMODE:
					sett.setScanMode( (InstrumentSettings::ScanMode)str2enum_(SCANMODEMAP,value) );
					break;
				case TIMEMIN:
					spec_->setRetentionTime(asFloat_(value)*60); //Minutes to seconds
					break;
				case TIMESEC:
					spec_->setRetentionTime(asFloat_(value));
					break;
				case POLARITY:
					sett.setPolarity( (IonSource::Polarity)str2enum_(POLARITYMAP,value) );
					break;
			  default:       error = "SpectrumDescription.SpectrumSettings.SpectrumInstrument";
				}
			}
			else if (is_parser_in_tag_[IONSELECTION]) {
				switch (ont){
				case MZ_ONT:
					spec_->getPrecursorPeak().getPosition()[0] = asFloat_(value);
					break;
				case CHARGESTATE:
					spec_->getPrecursorPeak().setCharge(asSignedInt_(value));
					break;
				case INTENSITY:
					spec_->getPrecursorPeak().getIntensity() = asFloat_(value);
					break;
				case IUNITS:
					setAddInfo(spec_->getPrecursorPeak(),"#IntensityUnits",value,
										 "Precursor.IonSelection.IntensityUnits");
					break;
				default:          error = "PrecursorList.Precursor.IonSelection.UserParam";
				}
			}
			else if (is_parser_in_tag_[ACTIVATION]) {
				switch (ont){
				case METHOD:
					prec_->setActivationMethod(
																		 (Precursor::ActivationMethod)str2enum_(ACTMETHODMAP,value));
					break;
				case ENERGY: prec_->setActivationEnergy( asFloat_(value) ); break;
				case EUNITS:
					prec_->setActivationEnergyUnit((Precursor::EnergyUnits)str2enum_(EUNITSMAP,value));
					break;
				default:     error = "PrecursorList.Precursor.Activation.UserParam";
				}
			}
			else if (useWarnings())
				warning(QXmlParseException(
																	 QString("Invalid cvParam: name=\"%1\", value=\"%2\"\n").arg(name).arg(value))
							 );

			if (useWarnings() && error != "")
				warning(QXmlParseException(
																	 QString("Invalid cvParam: name=\"%1\", value=\"%2\" in %3\n")
																	 .arg(name).arg(value).arg(error.c_str()))
							 );
		}

		template <typename MapType>
		void MzDataHandler<MapType>::fillData_()
		{
			if (data_.size() > decoder_.size()) // not enough decoder
				decoder_.resize(data_.size());

			std::vector<void*> ptrs(data_.size(),0);		//pointers to data of each decoder
			for (Size i=0; i<data_.size(); i++)
			{
				if (precisions_[i]==DOUBLE)	// precision 64 Bit
				{
					if (endians_[i]==BIG){
						ptrs[i] = static_cast<void*>(
																				 decoder_[i].decodeDoubleCorrected(data_[i].ascii(), data_[i].length()));
					}else{
						ptrs[i] = static_cast<void*>(
																				 decoder_[i].decodeDouble(data_[i].ascii(), data_[i].length()));
					}
				}else{											// precision 32 Bit
					if (endians_[i]==BIG){
						ptrs[i] = static_cast<void*>(
																				 decoder_[i].decodeFloatCorrected(data_[i].ascii(), data_[i].length()));
					}else{
						ptrs[i] = static_cast<void*>(
																				 decoder_[i].decodeFloat(data_[i].ascii(), data_[i].length()));
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
					spec_->insert(spec_->end(), PeakType());
					spec_->back().getIntensity() = getDatum(ptrs,INTENS,n);
					spec_->back().getPosition()[0] = getDatum(ptrs,MZ,n);
					//read supplemental data for derived classes (do nothing for DPeak)
					readPeakSupplementalData_(ptrs,spec_->back(),n);
				}
			}
		}

		template <typename MapType>
		void MzDataHandler<MapType>::writeTo(std::ostream& os)
		{
			os << "<!-- -*- Mode: XML; tab-width: 2; -*- -->\n<mzData version=\"1.05\" accessionNumber=\"OpenMS:\">\n";

			// delegate control to ExperimentalSettings handler
			Internal::MzDataExpSettHandler handler(*((const ExperimentalSettings*)cexp_));
			handler.writeTo(os);

			//determine how many spectra there are (count only those with peaks)
			UnsignedInt count_tmp_  = 0;
			for (UnsignedInt s=0; s<cexp_->size(); s++)
			{
				const SpectrumType& spec = (*cexp_)[s];
				if (spec.size()!=0) ++count_tmp_;
			}

			os << "\t<spectrumList count=\"" << count_tmp_ << "\">\n";

			for (UnsignedInt s=0; s<cexp_->size(); s++)
			{
				const SpectrumType& spec = (*cexp_)[s];

				//do not write empty spectra
				if (spec.size()==0) continue;

				os << "\t\t<spectrum id=\"" << spec_write_counter_++ << "\">\n"
					 << "\t\t\t<spectrumDesc>\n"
					 << "\t\t\t\t<spectrumSettings>\n";

				if (!spec.getAcquisitionInfo().empty())
				{
					os << "\t\t\t\t\t<acqSpecification spectrumType=\""
						 << spec.getAcquisitionInfo().getSpectrumType() << "\" methodOfCombination=\""
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
							<< "\t\t\t\t\t<precursor msLevel=\"1\" spectrumRef=\"0\">\n";
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
						writeUserParam_(os, peak);
						os << "\t\t\t\t\t\t</ionSelection>\n";
						os << "\t\t\t\t\t\t<activation>\n";
					}
					if (spec.getPrecursor() != Precursor())
					{
						const Precursor& prec = spec.getPrecursor();
						writeCVS_(os, prec.getActivationMethod(), ACTMETHODMAP, "1000044", "Method",6);
						writeCVS_(os, prec.getActivationEnergy(), "1000045", "CollisionEnergy",6);
						writeCVS_(os, prec.getActivationEnergyUnit(), EUNITSMAP,
											"1000046", "EnergyUnits",6);
						writeUserParam_(os, prec);
					}
						os << "\t\t\t\t\t\t</activation>\n";
					os << "\t\t\t\t\t</precursor>\n"
						 << "\t\t\t\t</precursorList>\n";
				}
				os << "\t\t\t</spectrumDesc>\n";

				// m/z
				float* tmp = decoder_[0].getFloatBuffer(spec.size());
				for (UnsignedInt i=0; i<spec.size(); i++)
					tmp[i] = spec.getContainer()[i].getPosition()[0];
				writeBinary(os,spec.size(),"mzArrayBinary");

				// intensity
				for (UnsignedInt i=0; i<spec.size(); i++)
					tmp[i] = spec.getContainer()[i].getIntensity();
				writeBinary(os,spec.size(),"intenArrayBinary");

				// write the supplementary data for picked peaks (is a no-op otherwise)
				this->writeDerivedPeakSupplementalData_(os, spec.getContainer());

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
		void MzDataHandler <MSExperiment<DPickedPeak<1,KernelTraits> > >::writeDerivedPeakSupplementalData_ < DPeakArrayNonPolymorphic<1, DPickedPeak<1,KernelTraits> > >( std::ostream& os, DPeakArrayNonPolymorphic<1, DPickedPeak<1,KernelTraits> > const & container);

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
