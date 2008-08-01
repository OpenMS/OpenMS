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
		: XMLHandler(filename,""),
  		exp_(&exp), 
  		cexp_(0)
	{
		cv_terms_.resize(15);
		// SampleState
		String(";Solid;Liquid;Gas;Solution;Emulsion;Suspension").split(';',cv_terms_[0]);
		// IonizationMode
		String(";PositiveIonMode;NegativeIonMode").split(';',cv_terms_[1]);
		// ResolutionMethod
		String(";FWHM;TenPercentValley;Baseline").split(';',cv_terms_[2]);
		// ResolutionType
		String(";Constant;Proportional").split(';',cv_terms_[3]);
		// ScanFunction
		String(";SelectedIonDetection;MassScan").split(';',cv_terms_[4]);
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
		String(";ProductIonScan;PrecursorIonScan;ConstantNeutralLoss;SingleReactionMonitoring;MultipleReactionMonitoring;SingleIonMonitoring;MultipleIonMonitoring").split(';',cv_terms_[12]);
		// DetectorType
		String(";EM;Photomultiplier;FocalPlaneArray;FaradayCup;ConversionDynodeElectronMultiplier;ConversionDynodePhotomultiplier;Multi-Collector;ChannelElectronMultiplier").split(';',cv_terms_[13]);
		// AnalyzerType
		String(";Quadrupole;PaulIonTrap;RadialEjectionLinearIonTrap;AxialEjectionLinearIonTrap;TOF;Sector;FourierTransform;IonStorage").split(';',cv_terms_[14]);
	}

   MzDataExpSettHandler::MzDataExpSettHandler(const ExperimentalSettings& exp, const String& filename)
		: XMLHandler(filename,""),
			exp_(0), 
			cexp_(&exp)
  {
		cv_terms_.resize(15);
		// SampleState
		String(";Solid;Liquid;Gas;Solution;Emulsion;Suspension").split(';',cv_terms_[0]);
		// IonizationMode
		String(";PositiveIonMode;NegativeIonMode").split(';',cv_terms_[1]);
		// ResolutionMethod
		String(";FWHM;TenPercentValley;Baseline").split(';',cv_terms_[2]);
		// ResolutionType
		String(";Constant;Proportional").split(';',cv_terms_[3]);
		// ScanFunction
		String(";SelectedIonDetection;MassScan").split(';',cv_terms_[4]);
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
		String(";ProductIonScan;PrecursorIonScan;ConstantNeutralLoss;SingleReactionMonitoring;MultipleReactionMonitoring;SingleIonMonitoring;MultipleIonMonitoring").split(';',cv_terms_[12]);
		// DetectorType
		String(";EM;Photomultiplier;FocalPlaneArray;FaradayCup;ConversionDynodeElectronMultiplier;ConversionDynodePhotomultiplier;Multi-Collector;ChannelElectronMultiplier").split(';',cv_terms_[13]);
		// AnalyzerType
		String(";Quadrupole;PaulIonTrap;RadialEjectionLinearIonTrap;AxialEjectionLinearIonTrap;TOF;Sector;FourierTransform;IonStorage").split(';',cv_terms_[14]);
	}

  MzDataExpSettHandler::~MzDataExpSettHandler()
  {	
  }

  void MzDataExpSettHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
  {
  	if (open_tags_.back()=="sampleName")
  	{
  		exp_->getSample().setName( sm_.convert(chars) );
  	}
  	else if (open_tags_.back()=="instrumentName")
  	{
  		exp_->getInstrument().setName(sm_.convert(chars));
  	}
  	else if (open_tags_.back()=="version")
  	{
  		exp_->getSoftware().setVersion( sm_.convert(chars) );
  	}
  	else if (open_tags_.back()=="institution")
  	{
  		exp_->getContacts().back().setInstitution( sm_.convert(chars) );
  	}
  	else if (open_tags_.back()=="contactInfo")
  	{
  		exp_->getContacts().back().setContactInfo( sm_.convert(chars) );
  	}
  	else if (open_tags_.back()=="nameOfFile")
  	{
  		exp_->getSourceFile().setNameOfFile( sm_.convert(chars) );
  	}
  	else if (open_tags_.back()=="pathToFile")
  	{
  		exp_->getSourceFile().setPathToFile( sm_.convert(chars) );
  	}
  	else if (open_tags_.back()=="fileType")
  	{
  		exp_->getSourceFile().setFileType( sm_.convert(chars) );
  	}
  	else if (open_tags_.back()=="comments")// <comment> is child of more than one other tags
  	{
  		String parent_tag = *(open_tags_.end()-2);

			if (parent_tag=="software")
			{
				exp_->getSoftware().setComment( sm_.convert(chars) );
			}
			else
			{
				warning(String("Unhandled tag \"comments\" with content: ") + sm_.convert(chars), 0, 0); 
			}
  	}
  	else if (open_tags_.back()=="name")
  	{
  		String parent_tag = *(open_tags_.end()-2);
  		
			if (parent_tag=="contact")
			{
				std::vector<String> tmp;
				if (String(sm_.convert(chars)).split(',',tmp))
				{
					exp_->getContacts().back().setFirstName(tmp[1]);
					exp_->getContacts().back().setLastName(tmp[0]);
				}
				else
				{
					if (String(sm_.convert(chars)).split(' ',tmp))
					{
						exp_->getContacts().back().setFirstName(tmp[0]);
						exp_->getContacts().back().setLastName(tmp[1]);
					}
					else
					{
						exp_->getContacts().back().setLastName(sm_.convert(chars));
					}
				}
			}
			else if (parent_tag=="software")
			{
				exp_->getSoftware().setName( sm_.convert(chars) );
			}
			else
			{
				warning(String("Unhandled tag \"name\" with content: ") + sm_.convert(chars), 0, 0); 
			}
  	}
  }
	
  void MzDataExpSettHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
  {
  	String tag = sm_.convert(qname);
  	//cout << "Start: '" << tag << "'" << endl;
  	
  	open_tags_.push_back(tag);

		// Do something depending on the tag
		if (tag=="cvParam")
		{
			cvParam_(attributes.getValue(sm_.convert("accession")),attributes.getValue(sm_.convert("value")));
		}
		else if (tag=="userParam")
		{
			//decode name and value
			String name_transcoded = attributeAsString_(attributes, sm_.convert("name"));
			String value_transcoded = "";
			optionalAttributeAsString_(value_transcoded, attributes, sm_.convert("value"));
	
			String& parent_tag = *(open_tags_.end()-2);
			if (parent_tag=="detector")
			{
				exp_->getInstrument().getIonDetector().setMetaValue(name_transcoded,value_transcoded);
			}
			else if (parent_tag=="source")
			{
				exp_->getInstrument().getIonSource().setMetaValue(name_transcoded,value_transcoded);
			}
			else if (parent_tag=="sampleDescription")
			{
				exp_->getSample().setMetaValue(name_transcoded,value_transcoded);
			}
			else if (parent_tag=="analyzer")
			{
				exp_->getInstrument().getMassAnalyzers().back().setMetaValue(name_transcoded,value_transcoded);
			}
			else if (parent_tag=="additional")
			{
				exp_->getInstrument().setMetaValue(name_transcoded,value_transcoded);
			}
			else if (parent_tag=="processingMethod")
			{
				exp_->getProcessingMethod().setMetaValue(name_transcoded,value_transcoded);			
			}
			else
			{
				warning(String("Invalid userParam: name=\"") + name_transcoded + "\", value=\"" + value_transcoded + "\"", 0, 0); 	
			}

		}
		else if (tag=="contact")
		{
			exp_->getContacts().resize(exp_->getContacts().size()+1);
		}
		else if (tag=="analyzer")
		{
			exp_->getInstrument().getMassAnalyzers().resize(exp_->getInstrument().getMassAnalyzers().size()+1);
		}
		else if (tag=="software")
		{
			if (attributes.getIndex(sm_.convert("completionTime"))!=-1)
			{
				exp_->getSoftware().setCompletionTime( asDateTime_(sm_.convert(attributes.getValue(sm_.convert("completionTime")))) );
			}
		}
  	//cout << "Exp - Start - End: '" << tag << "'" << endl;
	}

	void MzDataExpSettHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const /*qname*/)
  {	
  	//cout << "Exp: '" << sm_.convert(qname) << "'" << endl;
  	open_tags_.pop_back();
  }

	void MzDataExpSettHandler::cvParam_(const XMLCh* accession, const XMLCh* value)
	{
		//decode name and value
		String accession_transcoded = sm_.convert(accession);
		String value_transcoded = "";
		if (value != NULL)
		{
			value_transcoded = sm_.convert(value);
		}
		
		String error = "";
		String& parent_tag = *(open_tags_.end()-2);
		
		if (parent_tag=="detector")
		{
			if (accession_transcoded=="PSI:1000026")
			{
				exp_->getInstrument().getIonDetector().setType( (IonDetector::Type)cvStringToEnum_(13,value_transcoded, "detector type") );
			}
			else if (accession_transcoded=="PSI:1000028")
			{
				exp_->getInstrument().getIonDetector().setResolution( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000029")
			{
				exp_->getInstrument().getIonDetector().setADCSamplingFrequency( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000027")
			{
				exp_->getInstrument().getIonDetector().setAcquisitionMode((IonDetector::AcquisitionMode)cvStringToEnum_(9, value_transcoded, "acquisition mode") );
			}
			else
			{
				error = "Description.Instrument.Detector.UserParam";
			}
		}
		else if (parent_tag=="source")
		{
			if (accession_transcoded=="PSI:1000008")
			{
				exp_->getInstrument().getIonSource().setIonizationMethod( (IonSource::IonizationMethod)cvStringToEnum_(10, value_transcoded, "ion source") );
			}
			else if (accession_transcoded=="PSI:1000007")
			{
				exp_->getInstrument().getIonSource().setInletType( (IonSource::InletType)cvStringToEnum_(11, value_transcoded,"inlet type") );
			}
			else if (accession_transcoded=="PSI:1000009")
			{
				exp_->getInstrument().getIonSource().setPolarity( (IonSource::Polarity)cvStringToEnum_(1, value_transcoded,"polarity") );
			}
			else 
			{
				error = "Description.Instrument.Source.UserParam";
			}
		}
		else if (parent_tag=="sampleDescription")
		{
			if (accession_transcoded=="PSI:1000001")
			{
				exp_->getSample().setNumber( value_transcoded );
			}
			else if (accession_transcoded=="PSI:1000003")
			{
				exp_->getSample().setState( (Sample::SampleState)cvStringToEnum_(0, value_transcoded, "sample state") );
			}
			else if (accession_transcoded=="PSI:1000004")
			{
				exp_->getSample().setMass( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000005")
			{
				exp_->getSample().setVolume( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000006")
			{
				exp_->getSample().setConcentration( asFloat_(value_transcoded) );
			}
			else 
			{
				error = "Description.Admin.SampleDescription.UserParam";
			}
		}
		else if (parent_tag=="analyzer")
		{
			if (accession_transcoded=="PSI:1000010")
			{
				exp_->getInstrument().getMassAnalyzers().back().setType( (MassAnalyzer::AnalyzerType)cvStringToEnum_(14, value_transcoded,"analyzer type"));
			}
			else if (accession_transcoded=="PSI:1000011")
			{
				exp_->getInstrument().getMassAnalyzers().back().setResolution( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000012")
			{
				exp_->getInstrument().getMassAnalyzers().back().setResolutionMethod( (MassAnalyzer::ResolutionMethod)cvStringToEnum_(2, value_transcoded,"resolution method"));
			}
			else if (accession_transcoded=="PSI:1000013")
			{
				exp_->getInstrument().getMassAnalyzers().back().setResolutionType( (MassAnalyzer::ResolutionType)cvStringToEnum_(3, value_transcoded, "resolution type"));
			}
			else if (accession_transcoded=="PSI:1000014")
			{
				exp_->getInstrument().getMassAnalyzers().back().setAccuracy( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000015")
			{
				exp_->getInstrument().getMassAnalyzers().back().setScanRate( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000016")
			{
				exp_->getInstrument().getMassAnalyzers().back().setScanTime( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000017")
			{
				exp_->getInstrument().getMassAnalyzers().back().setScanFunction( (MassAnalyzer::ScanFunction)cvStringToEnum_(4, value_transcoded, "scan function"));
			}
			else if (accession_transcoded=="PSI:1000018")
			{
				exp_->getInstrument().getMassAnalyzers().back().setScanDirection( (MassAnalyzer::ScanDirection)cvStringToEnum_(5, value_transcoded, "scan direction"));
			}
			else if (accession_transcoded=="PSI:1000019")
			{
				exp_->getInstrument().getMassAnalyzers().back().setScanLaw( (MassAnalyzer::ScanLaw)cvStringToEnum_(6, value_transcoded, "scan law"));

			}
			else if (accession_transcoded=="PSI:1000020")
			{
				exp_->getInstrument().getMassAnalyzers().back().setTandemScanMethod((MassAnalyzer::TandemScanningMethod)cvStringToEnum_(12, value_transcoded, "tandem scanning mode"));
			}
			else if (accession_transcoded=="PSI:1000021")
			{
				exp_->getInstrument().getMassAnalyzers().back().setReflectronState( (MassAnalyzer::ReflectronState)cvStringToEnum_(8, value_transcoded, "reflectron state"));
			}
			else if (accession_transcoded=="PSI:1000022")
			{
				exp_->getInstrument().getMassAnalyzers().back().setTOFTotalPathLength( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000023")
			{
				exp_->getInstrument().getMassAnalyzers().back().setIsolationWidth( asFloat_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000024")
			{
				exp_->getInstrument().getMassAnalyzers().back().setFinalMSExponent( asInt_(value_transcoded) );
			}
			else if (accession_transcoded=="PSI:1000025")
			{
				exp_->getInstrument().getMassAnalyzers().back().setMagneticFieldStrength( asFloat_(value_transcoded) );
			}
			else 
			{
				error = "AnalyzerList.Analyzer.UserParam";
			}
		}
		else if (parent_tag=="additional")
		{
			if (accession_transcoded=="PSI:1000030")
			{
				exp_->getInstrument().setVendor(value_transcoded);
			}
			else if (accession_transcoded=="PSI:1000031")
			{
				exp_->getInstrument().setModel(value_transcoded);
			}
			else if (accession_transcoded=="PSI:1000032")
			{
				exp_->getInstrument().setCustomizations(value_transcoded);
			}
			else 
			{
				error = "Description.Instrument.Additional";
			}
		}
		else if (parent_tag=="processingMethod")
		{
			if (accession_transcoded=="PSI:1000033")
			{
				exp_->getProcessingMethod().setDeisotoping(asBool_(value_transcoded));
			}
			else if (accession_transcoded=="PSI:1000034")
			{
				exp_->getProcessingMethod().setChargeDeconvolution(asBool_(value_transcoded));
			}
			else if (accession_transcoded=="PSI:1000035")
			{
				exp_->getProcessingMethod().setSpectrumType( (SpectrumSettings::SpectrumType)cvStringToEnum_(7, value_transcoded, "spectrum type"));
			}
			else 
			{
				error = "DataProcessing.ProcessingMethod.UserParam";
			}
		}
		else
		{
			warning(String("Unexpected cvParam: accession=\"") + accession_transcoded + "\", value=\"" + value_transcoded + "\" in tag " + parent_tag, 0, 0); 
		}
		
		if (error != "")
		{
			warning(String("Invalid cvParam: accession=\"") + accession_transcoded +"\", value=\"" + value_transcoded +"\" in " + error, 0, 0); 
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

			if (cexp_->getSourceFile().getNameOfFile()!="")
			{
				os << "\t\t\t<sourceFile>\n"
					 << "\t\t\t\t<nameOfFile>" << cexp_->getSourceFile().getNameOfFile() << "</nameOfFile>\n"
					 << "\t\t\t\t<pathToFile>" << cexp_->getSourceFile().getPathToFile() << "</pathToFile>\n";
				if (cexp_->getSourceFile().getFileType()!="")
					os << "\t\t\t\t<fileType>" << cexp_->getSourceFile().getFileType() << "</fileType>\n";
				os << "\t\t\t</sourceFile>\n";
			}

			for (UInt i=0; i < cexp_->getContacts().size(); ++i)
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
			writeCVS_(os, inst.getIonSource().getInletType(),11, "1000007", "InletType");
			writeCVS_(os, inst.getIonSource().getIonizationMethod(), 10, "1000008","IonizationType");
			writeCVS_(os, inst.getIonSource().getPolarity(), 1, "1000009", "IonizationMode");
			writeUserParam_(os, inst.getIonSource());
			os << "\t\t\t</source>\n";
						
			//no analyzer given => add empty entry as there must be one entry
			UInt num_analyzers = inst.getMassAnalyzers().size();
			if (num_analyzers == 0)
			{
				os << "\t\t\t<analyzerList count=\"1\">\n"
				   << "\t\t\t\t<analyzer>\n"
				   << "\t\t\t\t</analyzer>\n";
			}
			else
			{
				os << "\t\t\t<analyzerList count=\"" << inst.getMassAnalyzers().size() << "\">\n";
				for (UInt i=0; i<inst.getMassAnalyzers().size(); ++i)
				{
					os << "\t\t\t\t<analyzer>\n";
					const MassAnalyzer& ana = inst.getMassAnalyzers()[i];
					writeCVS_(os, ana.getType(), 14, "1000010", "AnalyzerType",5);
					writeCVS_(os, ana.getResolution(), "1000011", "Resolution",5);
					writeCVS_(os, ana.getResolutionMethod(), 2,"1000012", "ResolutionMethod",5);
					writeCVS_(os, ana.getResolutionType(), 3, "1000013", "ResolutionType",5);
					writeCVS_(os, ana.getAccuracy(), "1000014", "Accuracy",5);
					writeCVS_(os, ana.getScanRate(), "1000015", "ScanRate",5);
					writeCVS_(os, ana.getScanTime(), "1000016", "ScanTime",5);
					writeCVS_(os, ana.getScanFunction(), 4,	"1000017", "ScanFunction",5);
					writeCVS_(os, ana.getScanDirection(), 5,	"1000018", "ScanDirection",5);
					writeCVS_(os, ana.getScanLaw(), 6, "1000019", "ScanLaw",5);
					writeCVS_(os, ana.getTandemScanMethod(), 12,"1000020", "TandemScanningMethod",5);
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
			writeCVS_(os, inst.getIonDetector().getType(), 13, "1000026", "DetectorType");
			writeCVS_(os, inst.getIonDetector().getAcquisitionMode(), 9, "1000027", "DetectorAcquisitionMode");
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
				writeCVS_(os, cexp_->getProcessingMethod().getSpectrumType(), 7, "1000035", "PeakProcessing");
				writeUserParam_(os, cexp_->getProcessingMethod());
				os << "\t\t\t</processingMethod>\n";
			}
			os << "\t\t</dataProcessing>\n"
				 << "\t</description>\n";
	}


	} // namespace Internal

} // namespace OpenMS

