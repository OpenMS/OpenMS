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

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

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
	
	@brief Shows basic information about the data in a file.
	
	This tool can show basic information about the data in several peak, feature and consensus feature files. It can
	- show information about the data range of a file (m/z, RT, intensity)
	- show a statistical summary for intensities and qualities
	- show an overview of the metadata
  - validate several XML formats against their XML schema
	- check for corrupt data in a file (e.g. dupliacte spectra)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileInfo
	: public TOPPBase
{
	public:
		TOPPFileInfo()
			: TOPPBase("FileInfo","Shows basic information about the file, such as data ranges and file type.")
		{
			
		}
	
	protected:

		virtual void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file ");
			setValidFormats_("in",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,cdf,mgf,featureXML,consensusXML"));
			registerStringOption_("in_type","<type>","","input file type -- default: determined from file extension or content\n", false);
			setValidStrings_("in_type",StringList::create("mzData,mzXML,DTA,DTA2D,cdf,mgf,featureXML,consensusXML"));
			registerFlag_("m","Show meta information about the whole experiment");
			registerFlag_("s","Computes a five-number statistics of intensities and qualities");
			registerFlag_("d","Show detailed listing of all spectra (peak files only)");
			registerFlag_("c","Check for corrupt data in the file (peak files only)");
			registerFlag_("v","Validate the file only (for mzData, mzXML, featureXML, IdXML, featurePairsXML, ConsensusXML)");
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//file names
			String in = getStringOption_("in");
			//file type
			FileHandler fh;
			FileHandler::Type in_type = fh.nameToType(getStringOption_("in_type"));
			
			if (in_type==FileHandler::UNKNOWN)
			{
				in_type = fh.getType(in);
				writeDebug_(String("Input file type: ") + fh.typeToName(in_type), 2);
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
			
			MSExperiment<Peak1D> exp;
			FeatureMap<> feat;
			ConsensusMap cons;
			ExperimentalSettings* exp_set;
			
			//-------------------------------------------------------------
			// validation
			//-------------------------------------------------------------
			if (getFlag_("v"))
			{
				bool valid = true;
				cout << endl << "Validating " << fh.typeToName(in_type) << " file";
				switch(in_type)
				{
					case FileHandler::MZDATA :
						cout << " against schema version " << MzDataFile().getVersion() << endl;
						valid = MzDataFile().isValid(in);
						break;
					case FileHandler::MZML :
						cout << " against schema version " << MzMLFile().getVersion() << endl;
						valid = MzMLFile().isValid(in);
						break;
					case FileHandler::FEATUREXML :
						cout << " against schema version " << FeatureXMLFile().getVersion() << endl;
						valid = FeatureXMLFile().isValid(in);
						break;
					case FileHandler::IDXML :
						cout << " against schema version " << IdXMLFile().getVersion() << endl;
						valid = IdXMLFile().isValid(in);
						break;
					case FileHandler::CONSENSUSXML :
						cout << " against schema version " << ConsensusXMLFile().getVersion() << endl;
						valid = ConsensusXMLFile().isValid(in);
						break;
					case FileHandler::MZXML :
						cout << " against schema version " << MzXMLFile().getVersion() << endl;
						valid = MzXMLFile().isValid(in);
						break;
					default:
						cout << endl << "Aborted: Validation of this file type is not supported!" << endl;
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

			Map<String,int> meta_names;			
			//-------------------------------------------------------------
			// Features
			//-------------------------------------------------------------
			if (in_type==FileHandler::FEATUREXML)
			{
				FeatureXMLFile().load(in,feat);
				feat.updateRanges();
				
				cout << "Number of features: " << feat.size() << endl
						 << endl
						 << "retention time range: " << feat.getMin()[Peak2D::RT] << " / " << feat.getMax()[Peak2D::RT] << endl
						 << "m/z range: " << feat.getMin()[Peak2D::MZ] << " / " << feat.getMax()[Peak2D::MZ] << endl
						 << "intensity range: " << feat.getMinInt() << " / " << feat.getMaxInt() << endl
						 << endl; 
		 
				exp_set = &feat;
			}
			//-------------------------------------------------------------
			// Consensus features
			//-------------------------------------------------------------
			else if (in_type==FileHandler::CONSENSUSXML)
			{
				ConsensusXMLFile().load(in,cons);
				cons.updateRanges();
				
				cout << "Number of conensus features: " << cons.size() << endl
						 << endl
						 << "retention time range: " << cons.getMin()[Peak2D::RT] << " / " << cons.getMax()[Peak2D::RT] << endl
						 << "m/z range: " << cons.getMin()[Peak2D::MZ] << " / " << cons.getMax()[Peak2D::MZ] << endl
						 << "intensity range: " << cons.getMinInt() << " / " << cons.getMaxInt() << endl
						 << endl; 
		 		
		 		//file descriptions
		 		const ConsensusMap::FileDescriptions& descs = cons.getFileDescriptions();
		 		if (descs.size()!=0)
		 		{
		 			cout << "File descriptions" << endl;
			 		for (ConsensusMap::FileDescriptions::const_iterator it=descs.begin(); it!=descs.end(); ++it)
			 		{
						cout << " - " << it->second.filename << endl
								 << "   identifier: " << it->first << endl
								 << "   label     : " << it->second.label << endl
								 << "   size      : " << it->second.size << endl;
			 		}
		 		}

				exp_set = new ExperimentalSettings();
			}
			//-------------------------------------------------------------
			// Peaks
			//-------------------------------------------------------------
			else
			{
				if (! fh.loadExperiment(in,exp,in_type,log_type_) )
				{
					writeLog_("Unsupported or corrupt input file. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;			
				}
				
				//determine type (search for the first scan with at least 5 peaks)
				UInt type = SpectrumSettings::UNKNOWN;
				UInt i=0; 
				while(i<exp.size() && exp[i].size()<5) ++i;
				if (i!=exp.size())
				{
					type = PeakTypeEstimator().estimateType(exp[i].begin(),exp[i].end());
				}
				cout << endl
						 << "peak type (metadata) : " << SpectrumSettings::NamesOfSpectrumType[exp.getProcessingMethod().getSpectrumType()] << endl
						 << "peak type (estimated): " << SpectrumSettings::NamesOfSpectrumType[type] << endl;
				//if raw data, determine the spacing
				if (type==SpectrumSettings::RAWDATA)
				{
					vector<Real> spacing;
					for (UInt j=1; j<exp[i].size(); ++j)
					{
						spacing.push_back(exp[i][j].getMZ()-exp[i][j-1].getMZ());
					}
					std::sort(spacing.begin(),spacing.end());
					cout << "raw data spacing: " << spacing[spacing.size()/2] << endl;	
				}
				cout << endl;
		
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
				for (MSExperiment<Peak1D>::iterator it = exp.begin(); it!=exp.end(); ++it)
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
				
				//show meta data array names
				for (MSExperiment<Peak1D>::iterator it = exp.begin(); it!=exp.end(); ++it)
				{
					for (i=0; i<it->getMetaDataArrays().size();++i)
					{
						String name = it->getMetaDataArrays()[i].getName();
						if (meta_names.has(name))
						{
							meta_names[name]++;
						}
						else
						{
							meta_names[name] = 1;
						}
					}
				}
				if (!meta_names.empty())
				{
					for (Map<String,int>::ConstIterator it=meta_names.begin();it!=meta_names.end();++it)
					{
						cout << "Meta data array: " << it->first << " (for " << it->second << " spectra)" << endl;
					}		
					cout << endl;
				}
				
				// Detailed listing of scans
				if (getFlag_("d"))
				{
					cout << endl 
					     << "-- Detailed spectrum listing --" << endl 
					     << endl;
					UInt count=0;
					for (MSExperiment<Peak1D>::iterator it = exp.begin(); it!=exp.end(); ++it)
					{
						++count;
						cout << "spectrum " << count << " - mslevel:" << it->getMSLevel() << " scanMode:" << InstrumentSettings::NamesOfScanMode[it->getInstrumentSettings().getScanMode()] << " peaks:" << it->size() << " RT:" << it->getRT() << " m/z:";
						if (it->size()!=0)
						{
							cout << it->begin()->getMZ() << "-" << (it->end()-1)->getMZ();
						}
						cout << endl;
					}
				}

				//Check for corrupt data
				if (getFlag_("c"))
				{
					cout << endl 
					     << "-- Checking for corrupt data --" << endl 
					     << endl;
					std::vector<DoubleReal> rts;
					rts.reserve(exp.size());
					for (UInt s=0; s<exp.size();++s)
					{
						//ms level = 0
						if (exp[s].getMSLevel()==0)
						{
							cout << "Error: MS-level 0 in spectrum (RT: " << exp[s].getRT() << ")" << std::endl;
						}
						//duplicate scans
						if (exp[s].getMSLevel()==1)
						{
							if (find(rts.begin(),rts.end(),exp[s].getRT())!=rts.end())
							{
								cout << "Error: Duplicate spectrum retention time: " << exp[s].getRT() << std::endl;
							}
							rts.push_back(exp[s].getRT());
						}
						//scan size = 0
						if (exp[s].size()==0)
						{
							cout << "Warning: No peaks in spectrum (RT: " << exp[s].getRT() << ")" << std::endl;
						}
						//duplicate meta data array names
						Map<String,int> names;
						for (UInt m=0; m<exp[s].getMetaDataArrays().size(); ++m)
						{
							String name = exp[s].getMetaDataArrays()[m].getName();
							if (names.has(name))
							{
								cout << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")" << std::endl;
							}
							else
							{
								names[name] = 0;
							}
						}
						//check peaks
						std::vector<DoubleReal> mzs;
						mzs.reserve(exp[s].size());
						for (UInt p=0; p<exp[s].size();++p)
						{
							//negative intensity
							if (exp[s][p].getIntensity()<0.0)
							{
								cout << "Warning: Negative intensity of peak (RT: " << exp[s].getRT() << " RT: " << exp[s][p].getMZ() << ")" << std::endl;
							}
							//duplicate m/z
							if (find(mzs.begin(),mzs.end(),exp[s][p].getMZ())!=mzs.end())
							{
								cout << "Warning: Duplicate peak m/z " << exp[s][p].getMZ() << " in spectrum (RT: " << exp[s].getRT() << ")" << std::endl;
							}
							mzs.push_back(exp[s][p].getMZ());
						}
					}
				}
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
				if (in_type==FileHandler::FEATUREXML) //features
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
		
					double mean, var, max, min;
					mean = gsl_stats_mean(intensities,1,size);
					var  = gsl_stats_variance(intensities,1,size);
					max  = gsl_stats_max(intensities,1,size);
					min  = gsl_stats_min(intensities,1,size);
		
					double mean_q, var_q, max_q, min_q;
					mean_q = gsl_stats_mean(qualities,1,size);
					var_q  = gsl_stats_variance(qualities,1,size);
					max_q  = gsl_stats_max(qualities,1,size);
					min_q  = gsl_stats_min(qualities,1,size);
		
					double median, upperq, lowerq;
					median = gsl_stats_median_from_sorted_data(intensities,1,size);
					upperq = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
					lowerq = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);
		
					double median_q, upperq_q, lowerq_q;
					median_q = gsl_stats_median_from_sorted_data(qualities,1,size);
					upperq_q = gsl_stats_quantile_from_sorted_data(qualities,1,size,0.75);
					lowerq_q = gsl_stats_quantile_from_sorted_data (qualities,1,size,0.25);
		
					delete [] intensities;
					delete [] qualities;
					
					cout << "Intensities:" << endl
							 << "  mean: " << QString::number(mean,'f',2).toStdString() << endl
							 << "  median: " << QString::number(median,'f',2).toStdString() << endl
							 << "  variance: " << QString::number(var,'f',2).toStdString() << endl
							 << "  min: " << QString::number(min,'f',2).toStdString() << endl
							 << "  max: " << QString::number(max,'f',2).toStdString() << endl
							 << "  lower_quartile: " << QString::number(lowerq,'f',2).toStdString() << endl
							 << "  upper_quartile: " << QString::number(upperq,'f',2).toStdString() << endl
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
				else if (in_type==FileHandler::CONSENSUSXML) //consensus features
				{
					UInt size = cons.size();
		
					DoubleReal* intensities = new DoubleReal[ size ];
					DoubleReal* qualities	 = new DoubleReal[ size ];
		
					for (unsigned int i = 0; i < size; 	++i)
					{
						intensities[i] = cons.at(i).getIntensity();
						qualities[i]   = cons.at(i).getQuality();
					}
		
					gsl_sort(intensities, 1, size);
					gsl_sort(qualities, 1, size);
		
					double mean, var, max, min;
					mean = gsl_stats_mean(intensities,1,size);
					var  = gsl_stats_variance(intensities,1,size);
					max  = gsl_stats_max(intensities,1,size);
					min  = gsl_stats_min(intensities,1,size);
		
					double mean_q, var_q, max_q, min_q;
					mean_q = gsl_stats_mean(qualities,1,size);
					var_q  = gsl_stats_variance(qualities,1,size);
					max_q  = gsl_stats_max(qualities,1,size);
					min_q  = gsl_stats_min(qualities,1,size);
		
					double median, upperq, lowerq;
					median = gsl_stats_median_from_sorted_data(intensities,1,size);
					upperq = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
					lowerq = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);
		
					double median_q, upperq_q, lowerq_q;
					median_q = gsl_stats_median_from_sorted_data(qualities,1,size);
					upperq_q = gsl_stats_quantile_from_sorted_data(qualities,1,size,0.75);
					lowerq_q = gsl_stats_quantile_from_sorted_data (qualities,1,size,0.25);
		
					delete [] intensities;
					delete [] qualities;
					
					cout << "Intensities:" << endl
							 << "  mean: " << QString::number(mean,'f',2).toStdString() << endl
							 << "  median: " << QString::number(median,'f',2).toStdString() << endl
							 << "  variance: " << QString::number(var,'f',2).toStdString() << endl
							 << "  min: " << QString::number(min,'f',2).toStdString() << endl
							 << "  max: " << QString::number(max,'f',2).toStdString() << endl
							 << "  lower_quartile: " << QString::number(lowerq,'f',2).toStdString() << endl
							 << "  upper_quartile: " << QString::number(upperq,'f',2).toStdString() << endl
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
				else //peaks
				{
					//copy intensities of  MS-level 1 peaks
					exp.updateRanges(1);
					UInt size = exp.getSize();
					DoubleReal* intensities = new  DoubleReal[ size ];
					UInt i = 0;
		      for (MSExperiment<Peak1D>::const_iterator spec = exp.begin(); spec != exp.end(); ++spec)
		      {
			      if (spec->getMSLevel()!=1)
			      {
			          continue;
			      }
			      for (MSExperiment<Peak1D>::SpectrumType::const_iterator it = spec->begin(); it!=spec->end(); ++it)
			      {
							intensities[i++] = it->getIntensity();
			      }
		      }
		
					gsl_sort(intensities, 1, size);
		
					double mean, var, max, min;
					mean = gsl_stats_mean(intensities,1,size);
					var  = gsl_stats_variance(intensities,1,size);
					max  = gsl_stats_max(intensities,1,size);
					min  = gsl_stats_min(intensities,1,size);
		
					double median, upperq, lowerq;
					median = gsl_stats_median_from_sorted_data(intensities,1,size);
					upperq = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
					lowerq = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);
		
					delete [] intensities;
		
					cout << "Intensities:" << endl
							 << "  mean: " << QString::number(mean,'f',2).toStdString() << endl
							 << "  median: " << QString::number(median,'f',2).toStdString() << endl
							 << "  variance: " << QString::number(var,'f',2).toStdString() << endl
							 << "  min: " << QString::number(min,'f',2).toStdString() << endl
							 << "  max: " << QString::number(max,'f',2).toStdString() << endl
							 << "  lower_quartile: " << QString::number(lowerq,'f',2).toStdString() << endl
							 << "  upper_quartile: " << QString::number(upperq,'f',2).toStdString() << endl
							 << endl;
				
					//Statistics for meta information
					for (Map<String,int>::ConstIterator it=meta_names.begin();it!=meta_names.end();++it)
					{
						String name = it->first;
						cout << "Meta data: " << name << endl;
						vector<Real> m_values;
						DoubleReal sum = 0.0;
			      for (MSExperiment<Peak1D>::const_iterator spec = exp.begin(); spec != exp.end(); ++spec)
			      {
			      	for (UInt meta=0; meta<spec->getMetaDataArrays().size(); ++meta)
			      	{
			      		if (spec->getMetaDataArrays()[meta].getName()!=name) continue;
					      for (UInt peak=0; peak < spec->size(); ++peak)
					      {
									m_values.push_back(spec->getMetaDataArrays()[meta][peak]);
					      	sum += spec->getMetaDataArrays()[meta][peak];
					      }
					    }
			      }
						sort(m_values.begin(),m_values.end());
						cout << "  count: " << m_values.size() << endl
								 << "  min: " << QString::number(m_values.front(),'f',2).toStdString() << endl
								 << "  max: " << QString::number(m_values.back(),'f',2).toStdString() << endl
								 << "  mean: " << QString::number(sum/m_values.size(),'f',2).toStdString() << endl
								 << "  median: " << QString::number(m_values[m_values.size()/2],'f',2).toStdString() << endl
								 << endl;
					}
				}

			}

			cout << endl << endl;			

			return EXECUTION_OK;
		}
};

int main( int argc, const char** argv )
{
	TOPPFileInfo tool;
	return tool.main(argc,argv);
}

/// @endcond
