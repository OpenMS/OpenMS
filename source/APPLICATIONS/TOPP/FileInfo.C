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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>

#include "TOPPBase.h"


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FileInfo FileInfo
	
	@brief Shows basic information about the data in an MS file.
	
	With this tool information about the data range of a file is displayed. It prints that m/z, intensity
	and retention time range that data lies in and some statistics about the number of spectra 
	for each MS level is displayed.
	
	Additionally an overview of the metadata of the experiment can be displayed.
	
	@todo determine if raw or picked data is contained (Marc)
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileInfo
	: public TOPPBase
{
	public:
		TOPPFileInfo()
			: TOPPBase("FileInfo")
		{
			
		}
	
	protected:
		void printToolUsage_()
		{
			cerr << endl
		       << tool_name_ << " -- shows basic information about the file e.g. data ranges and file type." << endl
		       << "Version: " << VersionInfo::getVersion() << endl
		       << endl
		       << "Usage:" << endl
					 << "  " << tool_name_ << " [options]" << endl
					 << endl
					 << "Options are:" << endl
					 << "  -in <file>        input file" << endl
					 << "  -in_type <type>   input file type (default: determined from input file extension)" << endl
					 << "  -m                show meta information about the whole experiment" << endl
					 << endl
					 << "Valid input types are: 'mzData', 'mzXML', 'DTA2D', 'ANDIMS' (cdf) , 'FeatureFile'" << endl;
		}
	
		void printToolHelpOpt_()
		{
			cerr << endl
		       << tool_name_ << endl
		       << endl
		       << "INI options:" << endl
					 << "  in        input file name" << endl
					 << "  in_type   input file type (default: determined from input file name extension)" << endl
					 << "  m         show meta information about the whole experiment" << endl
					 << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"in\" value=\"example.mzData\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"in_type\" value=\"MZDATA\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"m\" value=\"\" type=\"string\"/>" << endl;
		}
	
		void setOptionsAndFlags_()
		{
			options_["-in"] = "in";
			options_["-in_type"] = "in_type";
		
			flags_["-m"] = "m";
		}
	
		ExitCodes main_(int , char**)
		{
	
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//file names
			String in = getParamAsString_("in");
			writeDebug_(String("Input file: ") + in, 1);
			
			//file type
			FileHandler fh;
			FileHandler::Type in_type = fh.nameToType(getParamAsString_("in_type",""));
			
			writeDebug_(String("Input file type (from command line): ") + fh.typeToName(in_type), 1);
			
			if (in_type==FileHandler::UNKNOWN)
			{
				in_type = fh.getTypeByFileName(in);
				writeDebug_(String("Input file type (from file extention): ") + fh.typeToName(in_type), 1);
			}	

			if (in_type==FileHandler::UNKNOWN)
			{
				in_type = fh.getTypeByContent(in);
				writeDebug_(String("Input file type (from file content): ") + fh.typeToName(in_type), 1);
			}

			cout << endl
					 << "-- General information --" << endl
				   << endl
					 << "file name: " << in << endl
					 << "file type: " <<  fh.typeToName(in_type) << endl
					 << endl;
			
			MSExperiment< DPeak<1> > exp;
			DFeatureMap<2> feat;
			ExperimentalSettings* exp_set;
			//-------------------------------------------------------------
			// MSExperiment
			//-------------------------------------------------------------
			if (in_type!=FileHandler::FEATURE)
			{
			
				if (! fh.loadExperiment(in,exp,in_type) )
				{
					writeLog_("Unsupported or corrupt input file. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;			
				}
		
				//basic info
				exp.updateRanges();
				vector<UnsignedInt> levels = exp.getMSLevels();
										
				cout << "Number of peaks: " << exp.getSize() << endl 
						 << endl
						 << "retention time range: " << exp.getMinRT() << " / " << exp.getMaxRT() << endl
						 << "m/z range: " << exp.getMinMZ() << " / " << exp.getMaxMZ() << endl
						 << "intensity range: " << exp.getMinInt() << " / " << exp.getMaxInt() << endl
						 << "MS levels: ";
			  if (levels.size()!=0)
			  {
			  	cout  << *(levels.begin());
					for (vector<UnsignedInt>::iterator it = ++levels.begin(); it != levels.end(); ++it)
					{
						cout << ", " << *it;
					}
				}	 
				cout << endl << endl; 	
		
				//count how many spectra per MS level there are
				vector<UnsignedInt> counts(5);
				for (MSExperiment< DPeak<1> >::iterator it = exp.begin(); it!=exp.end(); ++it)
				{
					counts[it->getMSLevel()]++;	
				}
				//output
				for (UnsignedInt i = 0; i!=5; ++i)
				{
					if (counts[i]!=0)
					{
						cout << "Spectra of MS Level " << i << ": " << counts[i] << endl;
					}
				}
				cout << endl;
				
				exp_set = &exp;
			}
			//-------------------------------------------------------------
			// Feature
			//-------------------------------------------------------------
			else
			{
				DFeatureMapFile().load(in,feat);
				feat.updateRanges();
				
				UnsignedInt mz_dim = DimensionDescription< DimensionDescriptionTagLCMS >::MZ;
				UnsignedInt rt_dim = DimensionDescription< DimensionDescriptionTagLCMS >::RT;	
				
				cout 
						 << "Number of features: " << feat.size() << endl
						 << endl
						 << "retention time range: " << feat.getMin()[rt_dim] << " / " << feat.getMax()[rt_dim] << endl
						 << "m/z range: " << feat.getMin()[mz_dim] << " / " << feat.getMax()[mz_dim] << endl
						 << "intensity range: " << feat.getMinInt() << " / " << feat.getMaxInt() << endl
						 << endl; 
		 
				exp_set = &feat;
			}
			
			
			// '-m' show meta info
			if (getParamAsBool_("m"))
			{
				String date;
				exp_set->getDate().get(date);
				//basic info
				cout << endl
				     << "-- Meta information --" << endl
				     << endl
				     << "Experiment Type  : " << ExperimentalSettings::NamesOfExperimentType[exp_set->getType()] << endl
				     << "Date             : " <<  date << endl;		
				     
				//basic info
				cout << endl
				     << "Sample" << endl
				     << "  Name             : " << exp_set->getSample().getName() << endl
				     << "  Organism         : " << exp_set->getSample().getOrganism()  << endl
				     << "  Comment          : " << exp_set->getSample().getComment()  << endl;		

				//instrument info
				cout << endl
				     << "Instument" << endl
				     << "  Name             : " << exp_set->getInstrument().getName() << endl
				     << "  Model            : " << exp_set->getInstrument().getModel()  << endl	
				     << "  Vendor           : " << exp_set->getInstrument().getVendor()  << endl
				     << "  Ion source       : " << IonSource::NamesOfIonizationMethod[exp_set->getInstrument().getIonSource().getIonizationMethod()]  << endl
				     << "  Detector         : " << IonDetector::NamesOfType[exp_set->getInstrument().getIonDetector().getType()]  << endl
						 << "  Mass Analyzer(s) : ";
				
				for (UnsignedInt i=0; i< exp_set->getInstrument().getMassAnalyzers().size(); ++i)
				{
					cout  << MassAnalyzer::NamesOfAnalyzerType[exp_set->getInstrument().getMassAnalyzers()[i].getType()] << ", ";
				}
				cout << endl << endl;
				
				//contact persons
				for (UnsignedInt i=0; i< exp_set->getContacts().size(); ++i)
				{
					cout << "Contact Person" << endl
					     << "  Name             : " << exp_set->getContacts()[i].getName() << endl
					     << "  Email            : " << exp_set->getContacts()[i].getEmail() << endl
					     << endl;
				}
			}
			
			cout << endl << endl;			
			
			return OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPFileInfo tool;
	return tool.main(argc,argv);
}

/// @endcond
