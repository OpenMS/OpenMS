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
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FeaturePairsXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QtCore/QString>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

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
	
	Optionally an overview of the metadata of the map and a statistical summary of intensities can be displayed.

	The tool can also be used to validate several XML formats against their XML schema.
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
			                                            "Valid types are: 'dta', 'mzData', 'mzXML', 'DTA2D', 'ANDIMS' (cdf) , 'FeatureXML'", false);
			registerFlag_("m","Show meta information about the whole experiment");
			registerFlag_("s","Computes a five-number statistics of intensities (and feature qualities)");
			registerFlag_("d","Show detailed listing of all scans (for mzData only)");
			registerFlag_("v","Validate the file only.\nThis feature works for mzData, featureXML, IdXML, featurePairsXML and ConsensusXML.");
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
			
			MSExperiment<RawDataPoint1D> exp;
			FeatureMap<> feat;
			ExperimentalSettings* exp_set;
			
			//validation
			if (getFlag_("v"))
			{
				bool valid = true;
				cout << endl << "Validating " << fh.typeToName(in_type) << " file" << endl;
				switch(in_type)
				{
					case FileHandler::MZDATA :
						valid = MzDataFile().isValid(in);
						break;
					case FileHandler::FEATURE :
						valid = FeatureXMLFile().isValid(in);
						break;
					case FileHandler::FEATURE_PAIRS :
						valid = FeaturePairsXMLFile().isValid(in);
						break;
					case FileHandler::IDXML :
						valid = IdXMLFile().isValid(in);
						break;
					case FileHandler::CONSENSUSXML :
						valid = ConsensusXMLFile().isValid(in);
						break;
					default:
						cout << "Aborted: Validation of this file type is not supported!" << endl;
						return EXECUTION_OK;
				};
				
				if (valid)
				{
					cout << "Success: fhe file is valid!" << endl;
				}
				else
				{
					cout << "Failed: errors are listed above!" << endl;
				}
				
				return EXECUTION_OK;
			}
			
			//-------------------------------------------------------------
			// MSExperiment
			//-------------------------------------------------------------
			if (in_type!=FileHandler::FEATURE)
			{
			
				if (! fh.loadExperiment(in,exp,in_type,log_type_) )
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
				vector<UInt> levels = exp.getMSLevels();
										
				cout << "Number of peaks: " << exp.getSize() << endl 
						 << endl
						 << "retention time range: " << exp.getMinRT() << " / " << exp.getMaxRT() << endl
						 << "m/z range: " << exp.getMinMZ() << " / " << exp.getMaxMZ() << endl
						 << "intensity range: " << exp.getMinInt() << " / " << exp.getMaxInt() << endl
						 << "MS levels: ";
			  if (levels.size()!=0)
			  {
			  	cout  << *(levels.begin());
					for (vector<UInt>::iterator it = ++levels.begin(); it != levels.end(); ++it)
					{
						cout << ", " << *it;
					}
				}	 
				cout << endl << endl; 	
		
				//count how many spectra per MS level there are
				vector<UInt> counts(5);
				for (MSExperiment<RawDataPoint1D>::iterator it = exp.begin(); it!=exp.end(); ++it)
				{
					counts[it->getMSLevel()]++;	
				}
				//output
				for (UInt i = 0; i!=5; ++i)
				{
					if (counts[i]!=0)
					{
						cout << "Spectra of MS Level " << i << ": " << counts[i] << endl;
					}
				}
				cout << endl;
				
				exp_set = &exp;
				
				
				// Detailed listing of scans
				if (getFlag_("d"))
				{
					UInt count=0;
					for (MSExperiment<RawDataPoint1D>::iterator it = exp.begin(); it!=exp.end(); ++it)
					{
						++count;
						cout << "scan " << count << ": peaks: " << it->size() << " -- RT " << it->getRT() << " -- m/z ";
						if (it->size()!=0)
						{
							cout << it->begin()->getMZ() << "-" << (it->end()-1)->getMZ();
						}
						cout << endl;
					}
				}
			}
			//-------------------------------------------------------------
			// Feature
			//-------------------------------------------------------------
			else
			{
				FeatureXMLFile().load(in,feat);
				feat.updateRanges();
				
				cout 
						 << "Number of features: " << feat.size() << endl
						 << endl
						 << "retention time range: " << feat.getMin()[RawDataPoint2D::RT] << " / " << feat.getMax()[RawDataPoint2D::RT] << endl
						 << "m/z range: " << feat.getMin()[RawDataPoint2D::MZ] << " / " << feat.getMax()[RawDataPoint2D::MZ] << endl
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
				
				for (UInt i=0; i< exp_set->getInstrument().getMassAnalyzers().size(); ++i)
				{
					cout  << MassAnalyzer::NamesOfAnalyzerType[exp_set->getInstrument().getMassAnalyzers()[i].getType()] << ", ";
				}
				cout << endl << endl;
				
				//contact persons
				for (UInt i=0; i< exp_set->getContacts().size(); ++i)
				{
					cout << "Contact Person" << endl
					     << "  First Name       : " << exp_set->getContacts()[i].getFirstName() << endl
					     << "  Last Name        : " << exp_set->getContacts()[i].getLastName() << endl
					     << "  Email            : " << exp_set->getContacts()[i].getEmail() << endl
					     << endl;
				}
			}

			// '-s' show statistics
			if (getFlag_("s"))
			{
				cout << endl
				     << "-- Statistics --" << endl
				     << endl;
				if (in_type!=FileHandler::FEATURE) //peaks
				{
					//copy intensities of  MS-level 1 peaks
					exp.updateRanges(1);
					UInt size = exp.getSize();
					DoubleReal* intensities = new  DoubleReal[ size ];
					UInt i = 0;
		      for (MSExperiment<RawDataPoint1D>::const_iterator spec = exp.begin(); spec != exp.end(); ++spec)
		      {
			      if (spec->getMSLevel()!=1)
			      {
			          continue;
			      }
			      for (MSExperiment<RawDataPoint1D>::SpectrumType::const_iterator it = spec->begin(); it!=spec->end(); ++it)
			      {
							intensities[i++] = it->getIntensity();
			      }
			      
		      }
		
					gsl_sort(intensities, 1, size);
		
					double mean_int, var_int, max_int, min_int;
					mean_int = gsl_stats_mean(intensities,1,size);
					var_int  = gsl_stats_variance(intensities,1,size);
					max_int  = gsl_stats_max(intensities,1,size);
					min_int  = gsl_stats_min(intensities,1,size);
		
					double median_int, upperq_int, lowerq_int;
					median_int = gsl_stats_median_from_sorted_data(intensities,1,size);
					upperq_int = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
					lowerq_int = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);
		
					delete [] intensities;
		
					cout << "Intensities:" << endl
							 << "  mean: " << QString::number(mean_int,'f',2).toStdString() << endl
							 << "  median: " << QString::number(median_int,'f',2).toStdString() << endl
							 << "  variance: " << QString::number(var_int,'f',2).toStdString() << endl
							 << "  min: " << QString::number(min_int,'f',2).toStdString() << endl
							 << "  max: " << QString::number(max_int,'f',2).toStdString() << endl
							 << "  lower_quartile: " << QString::number(lowerq_int,'f',2).toStdString() << endl
							 << "  upper_quartile: " << QString::number(upperq_int,'f',2).toStdString() << endl
							 << endl;
							 
				}
				else //features
				{
					UInt size = feat.size();
		
					DoubleReal* intensities = new DoubleReal[ size ];
					DoubleReal* qualities	 = new DoubleReal[ size ];
		
					for (unsigned int i = 0; i < size; 	++i)
					{
						intensities[i] = feat.at(i).getIntensity();
						qualities[i]   = feat.at(i).getOverallQuality();
					}
		
					gsl_sort(intensities, 1, size);
					gsl_sort(qualities, 1, size);
		
					double mean_int, var_int, max_int, min_int;
					mean_int = gsl_stats_mean(intensities,1,size);
					var_int  = gsl_stats_variance(intensities,1,size);
					max_int  = gsl_stats_max(intensities,1,size);
					min_int  = gsl_stats_min(intensities,1,size);
		
					double mean_q, var_q, max_q, min_q;
					mean_q = gsl_stats_mean(qualities,1,size);
					var_q  = gsl_stats_variance(qualities,1,size);
					max_q  = gsl_stats_max(qualities,1,size);
					min_q  = gsl_stats_min(qualities,1,size);
		
					double median_int, upperq_int, lowerq_int;
					median_int = gsl_stats_median_from_sorted_data(intensities,1,size);
					upperq_int = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
					lowerq_int = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);
		
					double median_q, upperq_q, lowerq_q;
					median_q = gsl_stats_median_from_sorted_data(qualities,1,size);
					upperq_q = gsl_stats_quantile_from_sorted_data(qualities,1,size,0.75);
					lowerq_q = gsl_stats_quantile_from_sorted_data (qualities,1,size,0.25);
		
					delete [] intensities;
					delete [] qualities;
					
					cout << "Intensities:" << endl
							 << "  mean: " << QString::number(mean_int,'f',2).toStdString() << endl
							 << "  median: " << QString::number(median_int,'f',2).toStdString() << endl
							 << "  variance: " << QString::number(var_int,'f',2).toStdString() << endl
							 << "  min: " << QString::number(min_int,'f',2).toStdString() << endl
							 << "  max: " << QString::number(max_int,'f',2).toStdString() << endl
							 << "  lower_quartile: " << QString::number(lowerq_int,'f',2).toStdString() << endl
							 << "  upper_quartile: " << QString::number(upperq_int,'f',2).toStdString() << endl
							 << endl
							 << "Qualities:" << endl
							 << "  mean: " << QString::number(mean_q,'f',4).toStdString() << endl
							 << "  median: " << QString::number(median_q,'f',4).toStdString()  << endl
							 << "  variance: " << QString::number(var_q,'f',4).toStdString()  << endl
							 << "  min: " << QString::number(min_q,'f',4).toStdString()  << endl
							 << "  max: " << QString::number(max_q,'f',4).toStdString()  << endl
							 << "  lower_quartile: " << QString::number(lowerq_q,'f',4).toStdString()  << endl
							 << "  upper_quartile: " << QString::number(upperq_q,'f',4).toStdString()  << endl
							 ;
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
