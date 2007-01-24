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

#include <OpenMS/FORMAT/HANDLERS/MzDataExpSettHandler.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/ProcessingMethod.h>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
	namespace Internal
	{

	MzDataExpSettHandler::MzDataExpSettHandler(ExperimentalSettings& exp, const String& filename)
		: SchemaHandler(TAG_NUM,MAP_NUM,filename), // number of tags, number of maps
  		exp_(&exp), 
  		cexp_(0)
	{
		fillMaps_(Schemes::MzDataExpSett[schema_]);
		// fill maps with current schema
	}

   MzDataExpSettHandler::MzDataExpSettHandler(const ExperimentalSettings& exp, const String& filename)
		: SchemaHandler(TAG_NUM,MAP_NUM,filename), // number of tags, number of maps
			exp_(0), 
			cexp_(&exp)
  {
		fillMaps_(Schemes::MzDataExpSett[schema_]);
		// fill maps with current schema
	}

  MzDataExpSettHandler::~MzDataExpSettHandler()
  {	
  }

  void MzDataExpSettHandler::characters(const XMLCh* const chars, const unsigned int /*length*/)
  {
		// find the tag that the parser is in right now
 		for (Size i=0; i<is_parser_in_tag_.size(); i++)
			if (is_parser_in_tag_[i]){
				switch(i) {
					// Do something with the characters depending on the tag
					case SAMPLENAME_TAG:   exp_->getSample().setName( XMLString::transcode(chars) ); break;
					case INSTNAME:         exp_->getInstrument().setName(XMLString::transcode(chars)); break;
					case SWVERSION:        exp_->getSoftware().setVersion( XMLString::transcode(chars) ); break;
					case CONTACTINST:      contact_->setInstitution( XMLString::transcode(chars) ); break;
					case CONTACTINFO:      contact_->setContactInfo( XMLString::transcode(chars) ); break;

					case NAMEOFFILE:	exp_->getSourceFile().setNameOfFile( XMLString::transcode(chars) );	break;
					case PATHTOFILE:	exp_->getSourceFile().setPathToFile( XMLString::transcode(chars) );	break;
					case FILETYPE:		exp_->getSourceFile().setFileType( XMLString::transcode(chars) );	break;
					case COMMENTS:		// <comment> is child of more than one other tags
						if (is_parser_in_tag_[SOFTWARE])
						{
							exp_->getSoftware().setComment( XMLString::transcode(chars) );
						}
						else
						{
							const Locator* loc = 0;
							setDocumentLocator(loc);
							String tmp = String("Unhandled tag \"comments\" with content: ") + XMLString::transcode(chars);
							warning(SAXParseException(XMLString::transcode(tmp.c_str()), *loc )); 
						}
						break;
					case NAME: 	// <name> is child of more than one other tags
						if (is_parser_in_tag_[CONTACT])
						{
							std::vector<String> tmp;
							if (String(XMLString::transcode(chars)).split(',',tmp))
							{
								contact_->setFirstName(tmp[1]);
								contact_->setLastName(tmp[0]);
							}
							else
							{
								if (String(XMLString::transcode(chars)).split(' ',tmp))
								{
									contact_->setFirstName(tmp[0]);
									contact_->setLastName(tmp[1]);
								}
								else
								{
									contact_->setLastName(XMLString::transcode(chars));
								}
							}
						}
						else if (is_parser_in_tag_[SOFTWARE])
						{
							exp_->getSoftware().setName( XMLString::transcode(chars) );
						}
						else
						{
							const Locator* loc = 0;
							setDocumentLocator(loc);
							String tmp = String("Unhandled tag \"name\" with content: ") + XMLString::transcode(chars);
							warning(SAXParseException(XMLString::transcode(tmp.c_str()), *loc )); 
						}
						break;
				}
			}
  }
	
  void MzDataExpSettHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
  {
  	
  	//cout << "Exp - Start: '" << XMLString::transcode(qname) << "'" << endl;
		
		int tag = str2enum_(TAGMAP,XMLString::transcode(qname),"opening tag");	// index of current tag
		is_parser_in_tag_[tag] = true;

		// Do something depending on the tag
		switch(tag) 
		{
			case CVPARAM:	cvParam_(attributes.getValue(XMLString::transcode("name")),attributes.getValue(XMLString::transcode("value"))); break;
		  case USERPARAM:	userParam_(attributes.getValue(XMLString::transcode("name")),attributes.getValue(XMLString::transcode("value"))); break;
			case CONTACT:  contact_ = new ContactPerson(); break;
			case ANALYZER: analyzer_ = new MassAnalyzer(); break;
 			case SOFTWARE:
 				if (attributes.getIndex(XMLString::transcode("completionTime"))!=-1)
 				{
					exp_->getSoftware().setCompletionTime( asDateTime_(XMLString::transcode(attributes.getValue(XMLString::transcode("completionTime")))) );
				}
				break;
		}
	}



	void MzDataExpSettHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
  	
  	//cout << "Exp - End: '" << XMLString::transcode(qname) << "'" << endl;
  		
		int tag = str2enum_(TAGMAP,XMLString::transcode(qname),"closing tag");  // index of current tag
		is_parser_in_tag_[tag] = false;

		// Do something depending on the tag
		switch(tag) {
			case CONTACT:
				exp_->getContacts().push_back(*contact_);
				delete contact_;
				break;
			case ANALYZER:
				exp_->getInstrument().getMassAnalyzers().push_back(*analyzer_);
				delete analyzer_;
				break;
		}
  }


	void MzDataExpSettHandler::userParam_(const XMLCh* name, const XMLCh* value)
	{
		if (is_parser_in_tag_[DETECTOR])
			setAddInfo_(exp_->getInstrument().getIonDetector(), XMLString::transcode(name),XMLString::transcode(value),"Descr.Instrument.Detector.UserParam");
		else if (is_parser_in_tag_[INSTSRC])
			setAddInfo_(exp_->getInstrument().getIonSource(), XMLString::transcode(name),XMLString::transcode(value),"Descr.Instrument.Source.UserParam");
		else if (is_parser_in_tag_[SAMPLEDESCRIPTION])
			setAddInfo_(exp_->getSample(),XMLString::transcode(name),XMLString::transcode(value),"Descr.Admin.SampleDescription.UserParam");
		else if (is_parser_in_tag_[ANALYZER])
			setAddInfo_(*analyzer_,XMLString::transcode(name),XMLString::transcode(value),"AnalyzerList.Analyzer.UserParam");
		else if (is_parser_in_tag_[INSTADDITIONAL])
			setAddInfo_(exp_->getInstrument(),XMLString::transcode(name),XMLString::transcode(value),"Description.Instrument.Additional");
		else if (is_parser_in_tag_[PROCMETHOD])
			setAddInfo_(exp_->getProcessingMethod(), XMLString::transcode(name),XMLString::transcode(value),"DataProcessing.ProcessingMethod.UserParam");
		else
		{
			const Locator* loc = 0;
			setDocumentLocator(loc);
			String tmp = String("Invalid userParam: name=\"") + XMLString::transcode(name) + "\", value=\"" + XMLString::transcode(value) + "\"";
			warning(SAXParseException(XMLString::transcode(tmp.c_str()), *loc )); 
		}
	}


	void MzDataExpSettHandler::cvParam_(const XMLCh* name, const XMLCh* value)
	{
		int ont = str2enum_(ONTOLOGYMAP,XMLString::transcode(name),"cvParam elment"); // index of current ontology term

		std::string error = "";
		if (is_parser_in_tag_[DETECTOR]){
			IonDetector& ion_d = exp_->getInstrument().getIonDetector();
			switch (ont){
			case DETECTTYPE: ion_d.setType( (IonDetector::Type)str2enum_(TYPEMAP,XMLString::transcode(value)) ); break;
			case DETECTRES:  ion_d.setResolution( asFloat_(XMLString::transcode(value)) ); break;
			case ADCFREQ:    ion_d.setADCSamplingFrequency( asFloat_(XMLString::transcode(value)) ); break;
			case ACQMODE:
				ion_d.setAcquisitionMode((IonDetector::AcquisitionMode)str2enum_(ACQMODEMAP,XMLString::transcode(value)) );
				break;
			default:         error = "Description.Instrument.Detector.UserParam";
			}
		} else if (is_parser_in_tag_[INSTSRC]) {
			IonSource& ion_s = exp_->getInstrument().getIonSource();
			switch (ont) {
			case IONTYPE:   ion_s.setIonizationMethod( (IonSource::IonizationMethod)str2enum_(IONTYPEMAP,XMLString::transcode(value)) ); break;
			case INLETTYPE:	ion_s.setInletType( (IonSource::InletType)str2enum_(INLETTYPEMAP,XMLString::transcode(value)) ); break;
			case IONMODE:   ion_s.setPolarity( (IonSource::Polarity)str2enum_(IONMODEMAP,XMLString::transcode(value)) ); break;
			default:        error = "Description.Instrument.Source.UserParam";
			}
		}
		else if (is_parser_in_tag_[SAMPLEDESCRIPTION]) {
			Sample& sample = exp_->getSample();
			switch (ont){
			case SAMPLENAME_ONT: sample.setName( XMLString::transcode(value) ); break;
			case SAMPLESTATE:
				sample.setState( (Sample::SampleState)str2enum_(SAMPLESTATEMAP,XMLString::transcode(value)) );
				break;
			case SAMPLEMASS:     sample.setMass( asFloat_(XMLString::transcode(value)) ); break;
			case SAMPLEVOLUME:   sample.setVolume( asFloat_(XMLString::transcode(value)) ); break;
			case SAMPLECONC: 	   sample.setConcentration( asFloat_(XMLString::transcode(value)) ); break;
			case SAMPLENUMBER:   sample.setNumber( XMLString::transcode(value) ); break;
		  default:             error = "Description.Admin.SampleDescription.UserParam";
			}
		}
		else if (is_parser_in_tag_[ANALYZER]) {
			typedef MassAnalyzer MA;
			switch (ont){
			case ANALYZTYPE:
				analyzer_->setType( (MA::AnalyzerType)str2enum_(ANALYZERTYPEMAP,XMLString::transcode(value)));
				break;
			case RESOLUTION:	analyzer_->setResolution( asFloat_(XMLString::transcode(value)) ); break;
			case ACCURACY:	  analyzer_->setAccuracy( asFloat_(XMLString::transcode(value)) ); break;
			case SCANRATE:   	analyzer_->setScanRate( asFloat_(XMLString::transcode(value)) ); break;
			case SCANTIME:	  analyzer_->setScanTime( asFloat_(XMLString::transcode(value)) ); break;
			case TOFLENGTH:   analyzer_->setTOFTotalPathLength( asFloat_(XMLString::transcode(value)) ); break;
			case ISOWIDTH:    analyzer_->setIsolationWidth( asFloat_(XMLString::transcode(value)) ); break;
			case MAGSTRENGTH: analyzer_->setMagneticFieldStrength( asFloat_(XMLString::transcode(value)) ); break;
			case FINALMSEXP:	analyzer_->setFinalMSExponent( asSignedInt_(XMLString::transcode(value)) ); break;
			case RESMETHOD:
				analyzer_->setResolutionMethod( (MA::ResolutionMethod)str2enum_(RESMETHODMAP,XMLString::transcode(value)));
				break;
			case RESTYPE:
				analyzer_->setResolutionType( (MA::ResolutionType)str2enum_(RESTYPEMAP,XMLString::transcode(value)));
				break;
			case SCANFCT:
				analyzer_->setScanFunction( (MA::ScanFunction)str2enum_(SCANFUNCTIONMAP,XMLString::transcode(value)));
				break;
			case SCANDIR:
				analyzer_->setScanDirection( (MA::ScanDirection)str2enum_(SCANDIRECTIONMAP,XMLString::transcode(value)));
				break;
			case SCANLAW:
				analyzer_->setScanLaw( (MA::ScanLaw)str2enum_(SCANLAWMAP,XMLString::transcode(value)));
				break;
			case TANDEM:
				analyzer_->setTandemScanMethod((MA::TandemScanningMethod)str2enum_(TANDEMMAP,XMLString::transcode(value)));
				break;
			case REFLECTRON:
				analyzer_->setReflectronState( (MA::ReflectronState)str2enum_(REFLECTRONMAP,XMLString::transcode(value)));
				break;
			default:          error = "AnalyzerList.Analyzer.UserParam";
			}
		}
		else if (is_parser_in_tag_[INSTADDITIONAL]) {
			switch (ont){
			case VENDOR: exp_->getInstrument().setVendor(XMLString::transcode(value)); break;
			case MODEL:	 exp_->getInstrument().setModel(XMLString::transcode(value)); break;
			case CUSTOM: exp_->getInstrument().setCustomizations(XMLString::transcode(value)); break;
			default:     error = "Description.Instrument.Additional";
			}
		}
		else if (is_parser_in_tag_[PROCMETHOD]) {
			ProcessingMethod& meth = exp_->getProcessingMethod();
			switch (ont)	{
			case DEISOTOPED:  meth.setDeisotoping(asBool_(XMLString::transcode(value))); break;
			case DECONVOLVED: meth.setChargeDeconvolution(asBool_(XMLString::transcode(value))); break;
			case PEAKPROC:
				meth.setSpectrumType( (SpectrumSettings::SpectrumType)str2enum_(PEAKPROCMAP,XMLString::transcode(value)));
				break;
			default:          error = "DataProcessing.ProcessingMethod.UserParam";
			}
		}
		else
		{
			const Locator* loc = 0;
			setDocumentLocator(loc);
			String tmp = String("Invalid cvParam: name=\"") + XMLString::transcode(name) + "\", value=\"" + XMLString::transcode(value) + "\"";
			warning(SAXParseException(XMLString::transcode(tmp.c_str()), *loc )); 
		}
		
		if (error != "")
		{
			const Locator* loc = 0;
			setDocumentLocator(loc);
			String tmp = String("Invalid cvParam: name=\"") + XMLString::transcode(name) +"\", value=\"" + XMLString::transcode(value) +"\" in " + error;
			warning(SAXParseException(XMLString::transcode(tmp.c_str()), *loc )); 
		}
	}

	void MzDataExpSettHandler::writeTo(std::ostream& os)
	{
			const Sample& sm = cexp_->getSample();
			os << "\t<description>\n"
				 << "\t\t<admin>\n"
				 << "\t\t\t<sampleName>"
				 << sm.getName()
				 << "</sampleName>\n";

			if( sm.getNumber()!="" || sm.getState() || sm.getMass()
					|| sm.getVolume() || sm.getConcentration()	|| !sm.isMetaEmpty())
			{
				os << "\t\t\t<sampleDescription>\n";
				writeCVS_(os, sm.getNumber(), "1000001", "SampleNumber");
				//writeCVS_(os, sm.getName(), "1000002", "SampleName");  Already stored in <sampleName></sampleName>
				writeCVS_(os, sm.getState(), SAMPLESTATEMAP, "1000003", "SampleState");
				writeCVS_(os, sm.getMass(), "1000004", "SampleMass");
				writeCVS_(os, sm.getVolume(), "1000005", "SampleVolume");
				writeCVS_(os, sm.getConcentration(), "1000006", "SampleConcentration");
				writeUserParam_(os, cexp_->getSample());
				os << "\t\t\t</sampleDescription>\n";
			}

			if (cexp_->getSourceFile().getNameOfFile()!="")
			{
				os << "\t\t\t<sourceFile>\n"
					 << "\t\t\t\t<nameOfFile>" << cexp_->getSourceFile().getNameOfFile() << "</nameOfFile>\n"
					 << "\t\t\t\t<pathToFile>" << cexp_->getSourceFile().getPathToFile() << "</pathToFile>\n";
				if (cexp_->getSourceFile().getFileType()!="")
					os << "\t\t\t\t<fileType>" << cexp_->getSourceFile().getFileType() << "</fileType>\n";
				os << "\t\t\t</sourceFile>\n";
			}

			for (UnsignedInt i=0; i < cexp_->getContacts().size(); ++i)
			{
				os << "\t\t\t<contact>\n"
					 << "\t\t\t\t<name>" << cexp_->getContacts()[i].getFirstName() << " " << cexp_->getContacts()[i].getLastName() << "</name>\n"
					 << "\t\t\t\t<institution>" << cexp_->getContacts()[i].getInstitution() << "</institution>\n";
				if (cexp_->getContacts()[i].getContactInfo()!="")
					os << "\t\t\t\t<contactInfo>" << cexp_->getContacts()[i].getContactInfo() << "</contactInfo>\n";
				os << "\t\t\t</contact>\n";
			}
			os << "\t\t</admin>\n";
			const Instrument& inst = cexp_->getInstrument();
			os << "\t\t<instrument>\n"
				 << "\t\t\t<instrumentName>" << inst.getName() << "</instrumentName>\n"
				 << "\t\t\t<source>\n";
			writeCVS_(os, inst.getIonSource().getInletType(), INLETTYPEMAP, "1000007", "InletType");
			writeCVS_(os,inst.getIonSource().getIonizationMethod(),IONTYPEMAP,
								"1000008","IonizationType");
			writeCVS_(os, inst.getIonSource().getPolarity(), IONMODEMAP, "1000009", "IonizationMode");
			writeUserParam_(os, inst.getIonSource());

			os << "\t\t\t</source>\n"
				 << "\t\t\t<analyzerList count=\"" << inst.getMassAnalyzers().size() << "\">\n";
			for (UnsignedInt i=0; i<inst.getMassAnalyzers().size(); ++i)
			{
				os << "\t\t\t\t<analyzer>\n";
				const MassAnalyzer& ana = inst.getMassAnalyzers()[i];
				writeCVS_(os, ana.getType(), ANALYZERTYPEMAP, "1000010", "AnalyzerType",5);
				writeCVS_(os, ana.getResolution(), "1000011", "Resolution",5);
				writeCVS_(os, ana.getResolutionMethod(), RESMETHODMAP,"1000012", "ResolutionMethod",5);
				writeCVS_(os, ana.getResolutionType(), RESTYPEMAP, "1000013", "ResolutionType",5);
				writeCVS_(os, ana.getAccuracy(), "1000014", "Accuracy",5);
				writeCVS_(os, ana.getScanRate(), "1000015", "ScanRate",5);
				writeCVS_(os, ana.getScanTime(), "1000016", "ScanTime",5);
				writeCVS_(os, ana.getScanFunction(), SCANFUNCTIONMAP,	"1000017", "ScanFunction",5);
				writeCVS_(os, ana.getScanDirection(), SCANDIRECTIONMAP,	"1000018", "ScanDirection",5);
				writeCVS_(os, ana.getScanLaw(), SCANLAWMAP, "1000019", "ScanLaw",5);
				writeCVS_(os, ana.getTandemScanMethod(), TANDEMMAP,"1000020", "TandemScanningMethod",5);
				writeCVS_(os, ana.getReflectronState(), REFLECTRONMAP, "1000021", "ReflectronState",5);
				writeCVS_(os, ana.getTOFTotalPathLength(), "1000022", "TOFTotalPathLength",5);
				writeCVS_(os, ana.getIsolationWidth(), "1000023", "IsolationWidth",5);
				writeCVS_(os, ana.getFinalMSExponent(), "1000024", "FinalMSExponent",5);
				writeCVS_(os, ana.getMagneticFieldStrength(), "1000025", "MagneticFieldStrength",5);
				writeUserParam_(os, ana, 5);
				os << "\t\t\t\t</analyzer>\n";
			}
			os << "\t\t\t</analyzerList>\n"
				 << "\t\t\t<detector>\n";
			writeCVS_(os, inst.getIonDetector().getType(), TYPEMAP, "1000026", "DetectorType");
			writeCVS_(os, inst.getIonDetector().getAcquisitionMode(), ACQMODEMAP,
								"1000027", "DetectorAcquisitionMode");
			writeCVS_(os, inst.getIonDetector().getResolution(), "1000028", "DetectorResolution");
			writeCVS_(os, inst.getIonDetector().getADCSamplingFrequency(),
								"1000029", "ADCSamplingFrequency");
			writeUserParam_(os, inst.getIonDetector());
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

			os << "\t\t</instrument>\n"
				 << "\t\t<dataProcessing>\n"
				 << "\t\t\t<software";
			if (cexp_->getSoftware().getCompletionTime()!=DateTime())
			{
				String tmp;
				cexp_->getSoftware().getCompletionTime().get(tmp);
				String time(tmp);
				time.substitute(' ','T');
				os << " completionTime=\"" << time << "\"";
			}
			os << ">\n"
				 << "\t\t\t\t<name>" << cexp_->getSoftware().getName() << "</name>\n"
				 << "\t\t\t\t<version>" << cexp_->getSoftware().getVersion() << "</version>\n";
			if (cexp_->getSoftware().getComment()!="")
				os << "\t\t\t\t<comments>none</comments>\n";
			os << "\t\t\t</software>\n";

			if (cexp_->getProcessingMethod().getSpectrumType()!=0)
			{
				os << "\t\t\t<processingMethod>\n"
					 << "\t\t\t\t<cvParam cvLabel=\"psi\" name=\"Deisotoping\" accession=\"PSI:1000033\" value=\""
					 << ((cexp_->getProcessingMethod().getDeisotoping())? "true" : "false")
					 << "\"/>\n"
					 << "\t\t\t\t<cvParam cvLabel=\"psi\" name=\"ChargeDeconvolution\" accession=\"PSI:1000034\" value=\""
					 << ((cexp_->getProcessingMethod().getChargeDeconvolution())? "true" : "false")
					 << "\"/>\n";
				writeCVS_(os, cexp_->getProcessingMethod().getSpectrumType(), PEAKPROCMAP,
									"1000035", "PeakProcessing");
				writeUserParam_(os, cexp_->getProcessingMethod());
				os << "\t\t\t</processingMethod>\n";
			}
			os << "\t\t</dataProcessing>\n"
				 << "\t</description>\n";
	}


	} // namespace Internal

} // namespace OpenMS
