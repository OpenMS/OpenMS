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
// $Id: MzXMLHandler.h,v 1.7 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZXMLHANDLER_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <qtextstream.h>

namespace OpenMS
{
	namespace Internal
	{
  /**
  	@brief XML handlers for MzXMLFile

		MapType has to be a MSExperiment or have the same interface.
  		Do not use this class. It is only needed in MzXMLFile.
  */
	template <typename MapType>
  class MzXMLHandler
		: public SchemaHandler
  {
    public:
      /**@name Constructors and destructor */
      //@{
      ///
      MzXMLHandler(MapType& exp)
			: SchemaHandler(TAG_NUM,MAP_NUM), // number of tags, number of maps
				exp_(&exp),	cexp_(0),
				peak_(),
				spec_(0),
				analyzer_(0),
				decoder_(),
				peak_count_(0),
				spec_write_counter_(1)
  		{
				fillMaps_(Schemes::MzXML[schema_]);	// fill maps with current schema
			}

      ///
      MzXMLHandler(const MapType& exp)
			: SchemaHandler(TAG_NUM,MAP_NUM), // number of tags, number of maps
				exp_(0), cexp_(&exp),
				peak_(),
				decoder_(),
				peak_count_(0),
				spec_write_counter_(1)
  		{
				fillMaps_(Schemes::MzXML[schema_]);	// fill maps with current schema
			}

  		///
      virtual ~MzXMLHandler(){}
      //@}

      ///
      virtual bool endElement( const QString & uri, const QString & local_name,
															 const QString & qname ); 
      ///
      virtual bool startElement(const QString & uri, const QString & local_name,
																const QString & qname, const QXmlAttributes & attributes );
      ///
      virtual bool characters( const QString & chars );

  		///Write the contents to a stream
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
		enum Tags { TAGNULL=0, MSRUN, INDEX, OFFSET, SHA1, PARENTFILE, INSTRUMENT, DATAPROCESSING,
								SEPARATION, SPOTTING, SCAN, SCANORIGIN, PRECURSORMZ, MALDI, PEAKS, NAMEVALUE,
								COMMENT, SOFTWARE, INDEXOFFSET, OPERATOR, MANUFACTURER, MODEL, IONISATION, 									ANALYZER, DETECTOR, RESOLUTION, MZXML, PROCESSING, SEPARATIONTECH, TAG_NUM};


		/** @brief indices for attributes used by MzXML

			Used to access enum2str_() with ATTMAP.
			If you add terms, also add them to XMLSchemes.h
		*/
		enum Attributes {ATTNULL, POLARITY, SCANTYPE, CENTROIDED, DEISOTOPED,
										 DECONVOLUTED, RETTIME,IONENERGY, COLLENERGY, PRESSURE,
										 STARTMZ, ENDMZ, LOWMZ, HIGHMZ, BASEPEAKMZ, BASEPEAKINT,
										 TOTIONCURRENT, PEAKSCOUNT, NUM, MSLEVEL, SCANCOUNT,
										 FILENAME, FILETYPE, SOFTWAREVERSION, NAME, TYPE,
										 COMPLETION_TIME, PRECURSOR_INTENSITY, PRECURSOR_CHARGE,
										 FIRST_NAME, LAST_NAME, EMAIL, PHONE, URI, VALUE, CATEGORY,
										 PRECISION, BYTEORDER, PAIRORDER, SCHEMA, SPOT_INTEGRATION,
										 INTENSITY_CUTOFF};

		/** @brief indices for enum2str-maps used by mzXML

			Used to access enum2str_().
			If you add maps, also add them to XMLSchemes.h.
			Add no elements to the enum after MAP_NUM.
			Each map corresponds to a string in XMLSchemes.h.
		*/
		enum MapTypes {	POLARITYMAP=0, IONTYPEMAP, TYPEMAP, ANALYZERTYPEMAP, SCANMODEMAP,
										ATTMAP, TAGMAP, RESMETHODMAP, PEAKPROCMAP, PRECISIONMAP,MAP_NUM};

		/// Possible precisions for Base64 data encoding
		enum Precision { UNKNOWN_PRECISION, REAL, DOUBLE};

		typedef typename MapType::SpectrumType Spectrum;
		typedef typename MapType::SpectrumType::PrecursorPeakType PrecursorPeakType;

		/**@name temporary datastructures to hold parsed data */
    //@{
		typename MapType::SpectrumType::PeakType peak_;
		Spectrum* spec_;
		MassAnalyzer* analyzer_;
		MetaInfoDescription* meta_;
		QString meta_id_;
		Base64 decoder_;
		Size peak_count_;
		Precision precision_;
		const QXmlAttributes* atts_;
		//@}

		/// spectrum counter (spectra without peaks are not written)
		UnsignedInt spec_write_counter_;

		/// Add name, value and description to a given MetaInfo object
		void setAddInfo(MetaInfoInterface& info, const QString& name,
										const QString& value, const String& description)
		{
			info.metaRegistry().registerName(name.ascii(), description);
			info.setMetaValue(name.ascii(),value.ascii());
		}

		/// write metaInfo to xml (usually in nameValue-tag)
		inline void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta,
																int indent=4, String tag="nameValue")
		{
			std::vector<std::string> keys;  // Vector to hold keys to meta info
			meta.getKeys(keys);

			for (std::vector<std::string>::const_iterator it = keys.begin(); it!=keys.end(); ++it)
				if ( (*it)[0] != '#')  // internally used meta info start with '#'
			{
				String name = *it;
				os << String(indent,'\t') << "<" << tag << " name=\"";
				if (tag=="processingOperation")
				{
					std::vector<String> parts;
					name.split('#',parts);
					os << parts[0] << "\" type=\"" << parts[1];
				}
				else
				{
					os << name;
				}
				os << "\" value=\""
					 << meta.getMetaValue(*it) << "\"/>\n";
			}
		}

		/// check if value of attribute equals the required value, otherwise throw error
		inline void checkAttribute_(UnsignedInt attribute, const QString& required,
			const QString& required_alt="")
		{
			QString value = atts_->value(enum2str_(ATTMAP,attribute).c_str());
			if (value=="") return;
			if (value!=required && value!=required_alt)
			{
				no_error_ = false;
				error_message_ = QString("Unable to process data with %3 \"%1\" parsed by %2")
												 .arg(value).arg(file_).arg(enum2str_(ATTMAP,attribute));
			}
		}

		/// return value of attribute as QString
		inline QString getQAttribute(UnsignedInt attribute)
		{
			return atts_->value(enum2str_(ATTMAP,attribute).c_str());
		}

		/// return value of attribute as OpenMS::String
		inline String getAttribute(UnsignedInt attribute)
		{
			return atts_->value(enum2str_(ATTMAP,attribute).c_str()).ascii();
		}
  };







	//--------------------------------------------------------------------------------



	template <typename MapType>
  bool MzXMLHandler<MapType>::characters( const QString & chars )
  {
		if(is_parser_in_tag_[PEAKS])
		{

			if (precision_==DOUBLE)		//precision 64
			{
				double* data = decoder_.decodeDoubleCorrected(chars.ascii(), chars.length());
				//push_back the peaks into the container
				for (Size n = 0 ; n < (2 * peak_count_) ; n += 2)
				{
					peak_.getPosition()[0] = data[n];
					peak_.getIntensity() = data[n+1];
     				spec_->getContainer().push_back(peak_);
				}
			}else											//precision 32
			{
				float* data = decoder_.decodeFloatCorrected(chars.ascii(), chars.length());
				//push_back the peaks into the container
				for (Size n = 0 ; n < (2 * peak_count_) ; n += 2)
				{
					peak_.getPosition()[0] = data[n];
					peak_.getIntensity() = data[n+1];
     				spec_->getContainer().push_back(peak_);
				}
			}
		}
		else if (	is_parser_in_tag_[OFFSET] ||
							is_parser_in_tag_[INDEXOFFSET] ||
							is_parser_in_tag_[SHA1])
			;// do nothing
		else if (	is_parser_in_tag_[PRECURSORMZ])
		{
			spec_->getPrecursorPeak().getPosition()[0] = asFloat_(chars);
		}
		else if (	is_parser_in_tag_[COMMENT])
		{
			if (is_parser_in_tag_[INSTRUMENT])
			{
				 setAddInfo(exp_->getInstrument(),"#Comment" , chars, "Instrument.Comment");
			}else if (is_parser_in_tag_[DATAPROCESSING])
			{
				setAddInfo(exp_->getProcessingMethod(),"#Comment", chars,"DataProcessing.Comment");
			}else if (is_parser_in_tag_[SCAN])
			{
				spec_->setComment( chars.ascii() );
			}else if (useWarnings() && chars.stripWhiteSpace()!="")
				warning(QXmlParseException(QString("Unhandled characters: \"%1\"\n").arg(chars)));
		}
		else if (useWarnings() && chars.stripWhiteSpace()!="")
			warning(QXmlParseException(QString("Unhandled characters: \"%1\"\n").arg(chars)));

		return true;
  }

	template <typename MapType>
  	bool MzXMLHandler<MapType>::startElement(const QString & /*uri*/,
	const QString & /*local_name*/,	const QString & qname, const QXmlAttributes & attributes )
  {
		int tag = str2enum_(TAGMAP,qname,"opening tag");	// index of current tag
		is_parser_in_tag_[tag] = true;
		atts_ = &attributes;

		switch(tag)
			{
			case MSRUN: case MZXML:
				if (tag==MSRUN && !getQAttribute(SCANCOUNT).isEmpty())  // optional attribute
					exp_->reserve( asUnsignedInt_(getQAttribute(SCANCOUNT)) );

				// look for schema information
				if (!getQAttribute(SCHEMA).isEmpty())
				{
					QString s = getQAttribute(SCHEMA);
					for (UnsignedInt index=0; index<Schemes::MzXML_num; ++index)
					{
						if (s!=schema_ && s.contains(Schemes::MzXML[index][0].c_str()))
						{
							schema_ = index;
							// refill maps with older schema
							for (Size i=0; i<str2enum_array_.size(); i++)	str2enum_array_[i].clear();
							for (Size i=0; i<enum2str_array_.size(); i++)	enum2str_array_[i].clear();
							fillMaps_(Schemes::MzXML[schema_]);
						}
					}
				}
				// Additional attributes: startTime, endTime
				break;
			case PARENTFILE:
				exp_->getSourceFile().setNameOfFile( getAttribute(FILENAME).c_str() );
				exp_->getSourceFile().setFileType( getAttribute(FILETYPE) );
				// Additional attributes: fileSha1
				break;
			case INSTRUMENT:
				{
					if (attributes.length()==0) break;  // attributes only in mzXML 1.0

					exp_->getInstrument().setModel( attributes.value(enum2str_(TAGMAP,MODEL).c_str()).ascii());
					exp_->getInstrument().setVendor(attributes.value(
						enum2str_(TAGMAP,MANUFACTURER).c_str()).ascii()
					);

					MassAnalyzer analyzer;
					QString str = enum2str_(TAGMAP,ANALYZER);
					analyzer.setType(
						 (MassAnalyzer::AnalyzerType)str2enum_(ANALYZERTYPEMAP,atts_->value(str),str)
					);
					exp_->getInstrument().getMassAnalyzers().push_back(analyzer);
					str = enum2str_(TAGMAP,IONISATION);
					exp_->getInstrument().getIonSource().setIonizationMethod(
					 (IonSource::IonizationMethod)str2enum_(IONTYPEMAP,atts_->value(str),str));
				}
				break;
			case SOFTWARE:
				if (is_parser_in_tag_[DATAPROCESSING]){
					exp_->getSoftware().setVersion( getAttribute(SOFTWAREVERSION) );
					exp_->getSoftware().setName( getAttribute(NAME) );
					exp_->getSoftware().setComment( getAttribute(TYPE) );
					if  (!getQAttribute(COMPLETION_TIME).isEmpty())
						exp_->getSoftware().setCompletionTime(asFloat_(getQAttribute(COMPLETION_TIME)));
				}else if (is_parser_in_tag_[INSTRUMENT]){
					// not part of METADATA -> putting it into MetaInfo
					std::string swn = "#InstSoftware", swn_d = "Instrument software name",
						swv = "#InstSoftwareVersion", swv_d = "Instrument software version",
						swt = "#InstSoftwareType", swt_d = "Instrument software type",
						cmpl = "#InstSoftwareTime", cmpl_d = "Instrument software completion time";
					MetaInfoRegistry& registry =	MetaInfo().registry();
					registry.registerName(swn,swn_d);
					registry.registerName(swv,swv_d);
					registry.registerName(swt,swt_d);
					registry.registerName(cmpl,cmpl_d);
					exp_->getInstrument().setMetaValue(swn,getAttribute(NAME));
					exp_->getInstrument().setMetaValue(swv,getAttribute(SOFTWAREVERSION));
					exp_->getInstrument().setMetaValue(swt,getAttribute(TYPE));
					if  (!getQAttribute(COMPLETION_TIME).isEmpty())
						exp_->getInstrument().setMetaValue(cmpl,asFloat_(getQAttribute(COMPLETION_TIME)));
				}
				break;
			case PEAKS:
				{
					checkAttribute_(PRECISION,enum2str_(PRECISIONMAP,REAL),
																		enum2str_(PRECISIONMAP,DOUBLE));
					const QString str = enum2str_(ATTMAP,PRECISION);
					precision_ = (Precision) str2enum_(PRECISIONMAP,atts_->value(str),str);
					checkAttribute_(BYTEORDER,"network");
					checkAttribute_(PAIRORDER,"m/z-int");
				}
				break;
			case PRECURSORMZ:
				{
					typename MapType::SpectrumType::PrecursorPeakType& peak
						= spec_->getPrecursorPeak();
					peak.setIntensity( asFloat_(getQAttribute(PRECURSOR_INTENSITY)) );
					// optional attribute
					if (!getQAttribute(PRECURSOR_CHARGE).isEmpty())
						peak.setCharge(asSignedInt_(getQAttribute(PRECURSOR_CHARGE)));
					// Unhandled: windowWideness, precursorScanNum (optinal)
				}
				break;
			case SCAN:
				{
					exp_->push_back( Spectrum() );
					spec_ = &(*exp_)[exp_->size()-1];

					// required attributes
					peak_count_ = asSignedInt_( getQAttribute(PEAKSCOUNT) );
					spec_->setMSLevel( asSignedInt_(getQAttribute(MSLEVEL)) );

					//optinal attributes
					for (int i=0; i<attributes.length(); i++)
					{
						int att = str2enum_(ATTMAP,attributes.qName(i),"scan attribute");
						QString value = attributes.value(i);
						InstrumentSettings& sett = spec_->getInstrumentSettings();
						switch (att)
							{
							case POLARITY:
								sett.setPolarity( (IonSource::Polarity)
																		str2enum_(POLARITYMAP,value,"polarity") );
								break;
							case SCANTYPE:
								sett.setScanMode( (InstrumentSettings::ScanMode)
																		str2enum_(SCANMODEMAP,value,"scan mode") );
								break;
							case RETTIME:
								spec_->setRetentionTime(
																	asFloat_(value.remove(0,2).remove('S')));
								break;
							case STARTMZ:      sett.setMzRangeStart( asDouble_(value)); break;
							case ENDMZ:				 sett.setMzRangeStop( asDouble_(value)); break;
							case DEISOTOPED:
							  exp_->getProcessingMethod().setDeisotoping(asBool_(value));	break;
							case DECONVOLUTED:
								exp_->getProcessingMethod().setChargeDeconvolution(asBool_(value));	break;
							//	case CENTROIDED: case IONENERGY:	case COLLENERGY: 	case PRESSURE: case LOWMZ:
							//	case HIGHMZ: case BASEPEAKMZ:	case BASEPEAKINT:	case TOTIONCURRENT
							}
					}
				}
				break;
			case OPERATOR:
				{
					ContactPerson contact;
					contact.setName( QString("%2,%1").arg(getQAttribute(FIRST_NAME))
																					 .arg(getQAttribute(LAST_NAME)).ascii());
					if (!getQAttribute(EMAIL).isEmpty())
						contact.setEmail(getAttribute(EMAIL));
					contact.setContactInfo( QString("%1,%2").arg(getQAttribute(PHONE))
																									.arg(getQAttribute(URI)).ascii());
					exp_->getContacts().push_back(contact);
				}
				break;
			case MANUFACTURER:
				if (getAttribute(CATEGORY)==enum2str_(TAGMAP,MANUFACTURER))
					exp_->getInstrument().setVendor( getAttribute(VALUE) );
				break;
			case MODEL:
				if (getAttribute(CATEGORY)==enum2str_(TAGMAP,MODEL))
					exp_->getInstrument().setModel( getAttribute(VALUE) );
				break;
			case IONISATION:
				if (getAttribute(CATEGORY)==enum2str_(TAGMAP,IONISATION))
					exp_->getInstrument().getIonSource().setIonizationMethod(
						(IonSource::IonizationMethod)
						str2enum_(IONTYPEMAP,getQAttribute(VALUE),"ionization type")
					);
				break;
			case ANALYZER:
				if (getAttribute(CATEGORY)==enum2str_(TAGMAP,ANALYZER))
				{
					analyzer_ = new MassAnalyzer();
					analyzer_->setType( (MassAnalyzer::AnalyzerType)
						str2enum_(ANALYZERTYPEMAP,getQAttribute(VALUE),"analyzer type")
					);
				}
				break;
			case DETECTOR:
				if (getAttribute(CATEGORY)==enum2str_(TAGMAP,DETECTOR))
				{
					IonDetector& ion_d = exp_->getInstrument().getIonDetector();
					ion_d.setType( (IonDetector::Type)
						str2enum_(TYPEMAP,getQAttribute(VALUE),"detector type") );
				}
				break;
			case RESOLUTION:
				if (getAttribute(CATEGORY)==enum2str_(TAGMAP,RESOLUTION))
				{
					if (analyzer_ == 0) break;
					analyzer_->setResolutionMethod(
						(MassAnalyzer::ResolutionMethod)
						str2enum_(RESMETHODMAP,getQAttribute(VALUE),"resolution method"));
				}
				break;
			case DATAPROCESSING:
					//optinal attributes
					for (int i=0; i<attributes.length(); i++)
					{
						int att = str2enum_(ATTMAP,attributes.qName(i),"dataprocessing attribute");
						QString value = attributes.value(i);
						switch (att)
							{
								case DEISOTOPED:
								  exp_->getProcessingMethod().setDeisotoping(asBool_(value));
									break;
								case DECONVOLUTED:
									exp_->getProcessingMethod().setChargeDeconvolution(asBool_(value));
									break;
								case CENTROIDED:
									exp_->getProcessingMethod().setSpectrumType((SpectrumSettings::SpectrumType)
										str2enum_(PEAKPROCMAP,value,"peak processing"));
									break;
								//UNHANDLED: spotIntegration, intensityCutoff
							}
					}
					break;
			case NAMEVALUE:
				if (is_parser_in_tag_[INSTRUMENT])
				{
					setAddInfo(	exp_->getInstrument(), getQAttribute(NAME),
											getQAttribute(VALUE), "Instrument.Comment");
				}else if (is_parser_in_tag_[SCAN])
				{
					setAddInfo(	*spec_, getQAttribute(NAME),
											getQAttribute(VALUE), "Instrument.Comment");
				}else if (useWarnings())
					warning(QXmlParseException(QString("Unhandled tag %2.\n")
																			.arg(enum2str_(TAGMAP,NAMEVALUE).c_str())));
				break;
			case PROCESSING:
				setAddInfo(exp_->getProcessingMethod(),
						QString("%1#%2").arg(getQAttribute(NAME)).arg( getQAttribute(TYPE)),
						getQAttribute(VALUE), "Processing.Comment");

		}

		return no_error_;
	}


	template <typename MapType>
	bool MzXMLHandler<MapType>::endElement( const QString & /*uri*/,
	const QString & /*local_name*/,	const QString & qname )
  {
		int tag = str2enum_(TAGMAP,qname,"closing tag"); // index of current tag
		is_parser_in_tag_[tag] = false;

		if (tag==INSTRUMENT && analyzer_)
		{
			exp_->getInstrument().getMassAnalyzers().push_back(*analyzer_);
			delete analyzer_;
			analyzer_ = 0;
		}
		return true;
  }

	template <typename MapType>
	void MzXMLHandler<MapType>::writeTo(std::ostream& os)
	{
		//determine how many spectra there are (count only those with peaks)
		UnsignedInt count_tmp_  = 0;
		for (UnsignedInt s=0; s<cexp_->size(); s++)
		{
			const Spectrum& spec = (*cexp_)[s];
			if (spec.size()!=0) ++count_tmp_;
		}

		os << "<!-- -*- Mode: XML; tab-width: 2; -*- -->\n"
			 << "<mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_2.0\" "
		   << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
			 << "xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_2.0 "
			 << "http://sashimi.sourceforge.net/schema_revision/mzXML_2.0/mzXML_idx_2.0.xsd\">\n"
			 << "\t<msRun scanCount=\"" << count_tmp_ << "\">\n"
			 << "\t\t<parentFile fileName=\"" << cexp_->getSourceFile().getNameOfFile()
			 << "\" fileType=\"" << cexp_->getSourceFile().getFileType()
			 << "\" fileSha1=\"0000000000000000000000000000000000000000\"/>\n";

		if (cexp_->getInstrument() != Instrument())
		{
			const Instrument& inst = cexp_->getInstrument();
			os << "\t\t<msInstrument>\n"
				 << "\t\t\t<msManufacturer category=\"msManufacturer\" value=\""
				 <<	inst.getVendor() << "\"/>\n"
				 << "\t\t\t<msModel category=\"msModel\" value=\""
				 << inst.getModel() << "\"/>\n"
				 << "\t\t\t<msIonisation category=\"msIonisation\" value=\""
				 << enum2str_(IONTYPEMAP,inst.getIonSource().getIonizationMethod())
				 << "\"/>\n";

			const std::vector<MassAnalyzer>& analyzers = inst.getMassAnalyzers();
			if ( analyzers.size()>0 )
			{
				os << "\t\t\t<msMassAnalyzer category=\"msMassAnalyzer\" value=\""
					 << enum2str_(ANALYZERTYPEMAP,analyzers[0].getType())  << "\"/>\n";
			}
			else if (useWarnings())
				std::cerr << "Warning: mzXML supports only one analyzer! "
									<< "Skipping the other " << analyzers.size()
									<< " mass analyzers.\n";

			os << "\t\t\t<msDetector category=\"msDetector\" value=\""
				 << enum2str_(TYPEMAP,inst.getIonDetector().getType()) << "\"/>\n";
			try{
				std::string type = inst.getMetaValue("#InstSoftwareType").toString(),
				 	name = inst.getMetaValue("#InstSoftware"),
					version = inst.getMetaValue("#InstSoftwareVersion");
				float time = inst.getMetaValue("#InstSoftwareTime");
				os << "\t\t\t<software type=\"" << type
					 << "\" name=\"" << name
					 << "\" version=\"" << version
				 	 << "\" completionTime=\"" << time << "\"/>\n";
			}
			catch(Exception::InvalidValue exception)
			{

			}

			if ( analyzers.size()>0 )
			{
				if (analyzers[0].getResolutionMethod())
					os << "\t\t\t<msResolution category=\"msResolution\" value=\""
				 		 << enum2str_(RESMETHODMAP,analyzers[0].getResolutionMethod()) << "\"/>\n";
			}
			else if (useWarnings())
				std::cerr << "Warning: mzXML supports only one analyzer! "
									<< "Skipping the other " << analyzers.size()
									<< " mass analyzers.\n";

			if ( cexp_->getContacts().size()>0 )
			{
				const ContactPerson& cont = cexp_->getContacts()[0];
				std::vector<String> name;
				cont.getName().split(',',name);
				os << "\t\t\t<operator first=\"" << name[1]
					 << "\" last=\"" << name[0];

				std::vector<String> info;
				cont.getContactInfo().split(',',info);
				os << "\" phone=\"" << info[0]
					 << "\" email=\"" << cont.getEmail()
					 << "\" URI=\"" << info[1] << "\"/>\n";
			}
			writeUserParam_(os,inst,3);
			DataValue com = inst.getMetaValue("#Comment");
			if (!com.isEmpty())
				os << "\t\t\t<comment>" << com << "</comment>\n";
			os << "\t\t</msInstrument>\n";
		}

		const Software& software = cexp_->getSoftware();
		os << "\t\t<dataProcessing deisotoped=\""
			 << cexp_->getProcessingMethod().getDeisotoping()
			 << "\" chargeDeconvoluted=\""
			 << cexp_->getProcessingMethod().getChargeDeconvolution()
			 << "\" centroided=\""
			 << enum2str_(PEAKPROCMAP,cexp_->getProcessingMethod().getSpectrumType())
			 << "\">\n"
			 << "\t\t\t<software type=\"" << software.getComment()
			 << "\" name=\"" << software.getName()
			 << "\" version=\"" << software.getVersion();

		if (software.getCompletionTime() != float())
			os << "\" completionTime=\"" << software.getCompletionTime();
		os << "\"/>\n";
		writeUserParam_(os,cexp_->getProcessingMethod(),3,"processingOperation");

		DataValue com = cexp_->getProcessingMethod().getMetaValue("#Comment");
		if (!com.isEmpty())
			os << "\t\t\t<comment>" << com << "</comment>\n";
		os << "\t\t</dataProcessing>\n";

		// write scans
		for (UnsignedInt s=0; s<cexp_->size(); s++)
		{
			const Spectrum& spec = (*cexp_)[s];

			//do not write empty spectra
			if (spec.size()==0) continue;
			
			int MSLevel = spec.getMSLevel();

			if (MSLevel==1 && s!=0)
				os << String(MSLevel+1,'\t') << "</scan>\n";

			os << String(MSLevel+1,'\t')
				 << "<scan num=\"" << spec_write_counter_++ << "\" msLevel=\""
				 << spec.getMSLevel() << "\" peaksCount=\""
				 << spec.size() << "\" polarity=\""
				 << enum2str_(POLARITYMAP,spec.getInstrumentSettings().getPolarity());

			if (spec.getInstrumentSettings().getScanMode())
			{
				os << "\" scanType=\""
					 << enum2str_(SCANMODEMAP,spec.getInstrumentSettings().getScanMode());
			}
			os << "\" retentionTime=\"PT"
				 << spec.getRetentionTime() << "S\"";
			if (spec.getInstrumentSettings().getMzRangeStart()!=0)
				os << " startMz=\"" << spec.getInstrumentSettings().getMzRangeStart() << "\"";
			if (spec.getInstrumentSettings().getMzRangeStop()!=0)
				os << " endMz=\"" << spec.getInstrumentSettings().getMzRangeStop() << "\"";
			os << ">\n";

			const PrecursorPeakType& peak = spec.getPrecursorPeak();
			if (peak!= PrecursorPeakType())
			{
				os << String(MSLevel+2,'\t') << "<precursorMz precursorIntensity=\""
					 << peak.getIntensity();
				if (peak.getCharge()!=0)
					os << "\" precursorCharge=\"" << peak.getCharge();
				os << "\">"
				 	 << peak.getPosition()[0] << "</precursorMz>\n";
			}

			os << String(MSLevel+2,'\t') << "<peaks precision=\"32\""
				 << " byteOrder=\"network\" pairOrder=\"m/z-int\">";

			float* tmp = decoder_.getFloatBuffer(spec.size()*2);
			for (UnsignedInt i=0; i<spec.size(); i++)
			{
				tmp[2*i]   = spec.getContainer()[i].getPosition()[0];
				tmp[2*i+1] = spec.getContainer()[i].getIntensity();
			}
			os << decoder_.encodeFloatCorrected() << "</peaks>\n";

			writeUserParam_(os,spec,MSLevel+2);
			if (spec.getComment() != "")
				os << String(MSLevel+2,'\t') << "<comment>" << spec.getComment() << "</comment>\n";

			if (MSLevel==2)
				os << String(MSLevel+1,'\t') << "</scan>\n";
		}

		if (cexp_->size()) os << "\t\t</scan>\n";
		os << "\t</msRun>\n"
			 << "\t<indexOffset>0</indexOffset>\n"
			 << "</mzXML>\n";
	}

	} // namespace Internal

} // namespace OpenMS

#endif
