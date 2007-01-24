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

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>


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
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileInfo
	: public TOPPBase
{
	public:
		TOPPFileInfo()
			: TOPPBase("FileInfo","Shows basic information about the file e.g. data ranges and file type")
		{
			
		}
	
	protected:

		virtual void registerOptionsAndFlags_()
		{
			registerStringOption_("in","<file>","","input file");
			registerStringOption_("in_type","<type>","","input file type (default: determined from file extension or content)\n"
			                                            "Valid types are: 'mzData', 'mzXML', 'DTA2D', 'ANDIMS' (cdf) , 'FeatureFile'", false);
			registerFlag_("m","show meta information about the whole experiment");
		}
		
		ExitCodes main_(int , char**)
		{
	
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//file names
			String in = getStringOption_("in");
			inputFileReadable_(in);
			//file type
			FileHandler fh;
			FileHandler::Type in_type = fh.nameToType(getStringOption_("in_type"));
			
			if (in_type==FileHandler::UNKNOWN)
			{
				in_type = fh.getTypeByFileName(in);
				writeDebug_(String("Input file type (from file extention): ") + fh.typeToName(in_type), 2);
			}

			if (in_type==FileHandler::UNKNOWN)
			{
				in_type = fh.getTypeByContent(in);
				writeDebug_(String("Input file type (from content): ") + fh.typeToName(in_type), 2);
			}

			if (in_type==FileHandler::UNKNOWN)
			{
				writeLog_("Error: Could not determine input file type!");
				return PARSE_ERROR;
			}

			cout << endl
					 << "-- General information --" << endl
				   << endl
					 << "file name: " << in << endl
					 << "file type: " <<  fh.typeToName(in_type) << endl;
			
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

				cout << endl
						 << "peak type (metadata) : " << SpectrumSettings::NamesOfSpectrumType[exp.getProcessingMethod().getSpectrumType()] << endl
						 << "peak type (estimated): " << SpectrumSettings::NamesOfSpectrumType[PeakTypeEstimator().estimateType(exp[0].begin(),exp[0].end())] << endl
						 << endl;
		
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
				
				UnsignedInt mz_dim = DimensionDescription< LCMS_Tag >::MZ;
				UnsignedInt rt_dim = DimensionDescription< LCMS_Tag >::RT;	
				
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
			if (getFlag_("m"))
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
					     << "  First Name       : " << exp_set->getContacts()[i].getFirstName() << endl
					     << "  Last Name        : " << exp_set->getContacts()[i].getLastName() << endl
					     << "  Email            : " << exp_set->getContacts()[i].getEmail() << endl
					     << endl;
				}
			}
			
			cout << endl << endl;			
			
			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPFileInfo tool;
	return tool.main(argc,argv);
}

/// @endcond
