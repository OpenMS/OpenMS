// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FileFilter FileFilter
	
	@brief Extracts portions of the data from an mzData, featureXML or consensusXML file.
	
	With this tool it is possible to extract m/z, retention time and intensity ranges from an input file
	and to write all data that lies within the given ranges to an output file.
	
	Depending on the input file type, additional specific operations are possible:
	- mzData 
		- extract spectra of a certain MS level
		- filter by signal-to-noise estimation
		- filter by scan mode of the spectra
  - featureXML
    - filter by feature charge
    - filter by overall feature quality

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileFilter
	: public TOPPBase
{
	public:
		TOPPFileFilter()
			: TOPPBase("FileFilter","extracts/modifies portions of data from an mzData, featureXML or consensusXML file")
		{
		}
	
	protected:

		typedef MSExperiment<Peak1D> MapType;

		void registerOptionsAndFlags_()
		{
      registerInputFile_("in","<file>","","input file ");
   		setValidFormats_("in",StringList::create("mzData,featureXML,consensusXML"));

      registerOutputFile_("out","<file>","","output file");
	  	setValidFormats_("out",StringList::create("mzData,featureXML,consensusXML"));
      
			registerStringOption_("mz","[min]:[max]",":","m/z range to extract", false);
			registerStringOption_("rt","[min]:[max]",":","retention time range to extract", false);
			registerStringOption_("int","[min]:[max]",":","intensity range to extract", false);
      
			addText_("peak data options:");
      registerDoubleOption_("sn", "<s/n ratio>", 0, "write peaks with S/N > 'sn' values only", false);
			registerStringOption_("level","i[,j]...","1,2,3","MS levels to extract", false);
			registerStringOption_("remove_mode","<mode>","","Remove scans by scan mode\n",false);
			StringList mode_list;
			for (UInt i=0; i<InstrumentSettings::SIZE_OF_SCANMODE; ++i)
			{
				mode_list.push_back(InstrumentSettings::NamesOfScanMode[i]);
			}
			setValidStrings_("remove_mode",mode_list);
      registerFlag_("sort","sorts the output data according to RT and m/z."
      										 "\nNote: Spectrum meta data arrays are erased, as they would be invalid after sorting by m/z.");
      
      addText_("feature data options:");
      registerStringOption_("charge","[min]:[max]",":","charge range to extract", false);
      registerStringOption_("q","[min]:[max]",":","OverallQuality range to extract [0:1]", false);

			addEmptyLine_();
			addText_("Other options of the FileFilter only apply if S/N estimation is done.\n"
							 "They can be given only in the 'algorithm' section  of the INI file.");
	    
			registerSubsection_("algorithm","S/N algorithm section");

		}
	
		Param getSubsectionDefaults_(const String& /*section*/) const
		{
			SignalToNoiseEstimatorMedian<  MapType::SpectrumType > sn;
			Param tmp;
			tmp.insert("SignalToNoise:",sn.getParameters());
			return tmp;
		}

		ExitCodes main_(int , const char**)
		{

			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			String in = getStringOption_("in");
			String out = getStringOption_("out");
        
      //input file type
      FileHandler::Type in_type = FileHandler::getType(in);
      writeDebug_(String("Input file type: ") + FileHandler::typeToName(in_type), 2);
  
      if (in_type==FileHandler::UNKNOWN)
      {
        writeLog_("Error: Could not determine input file type!");
        return PARSE_ERROR;
      }
	
      FileHandler::Type out_type = in_type;            

			//ranges
			String mz, rt, it, level, charge, q, tmp;
			double mz_l, mz_u, rt_l, rt_u, it_l, it_u, sn, charge_l, charge_u, q_l, q_u;
			vector<UInt> levels;		
			//initialize ranges
			mz_l = rt_l = it_l = charge_l = q_l = -1 * numeric_limits<double>::max();
			mz_u = rt_u = it_u = charge_u = q_u = numeric_limits<double>::max();
			
			rt = getStringOption_("rt");
			mz = getStringOption_("mz");
			it = getStringOption_("int");
			level = getStringOption_("level");
			sn = getDoubleOption_("sn");
			charge = getStringOption_("charge");
      q = getStringOption_("q");
      
			//convert bounds to numbers
			try
			{
				//rt
				parseRange_(rt,rt_l,rt_u);
				writeDebug_("rt lower/upper bound: " + String(rt_l) + " / " + String(rt_u),1);	
				
				//mz
				parseRange_(mz,mz_l,mz_u);
				writeDebug_("mz lower/upper bound: " + String(mz_l) + " / " + String(mz_u),1);	
				
				//int
				parseRange_(it,it_l,it_u);
				writeDebug_("int lower/upper bound: " + String(it_l) + " / " + String(it_u),1);	
	
				//levels
				tmp = level;
				if (level.has(',')) //several levels given
				{
					vector<String> tmp2;
					level.split(',',tmp2);
					for (vector<String>::iterator it = tmp2.begin(); it != tmp2.end(); ++it)
					{
						levels.push_back(it->toInt());
					}
				}
				else //one level given
				{
					levels.push_back(level.toInt());
				}
				
				String tmp3("MS levels: ");
				tmp3 = tmp3 + *(levels.begin());
				for (vector<UInt>::iterator it = ++levels.begin(); it != levels.end(); ++it)
				{
					tmp3 = tmp3 + ", " + *it;
				}
				writeDebug_(tmp3,1);	


        //charge (features only)
        parseRange_(charge,charge_l,charge_u);
        writeDebug_("charge lower/upper bound: " + String(charge_l) + " / " + String(charge_u),1); 

        //charge (features only)
        parseRange_(q,q_l,q_u);
        writeDebug_("quality lower/upper bound: " + String(q_l) + " / " + String(q_u),1); 

			}
			catch(Exception::ConversionError&)
			{
				writeLog_(String("Invalid boundary '") + tmp + "' given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;			
			}
			
      
      if (in_type == FileHandler::MZDATA)
      {
  			//-------------------------------------------------------------
  			// loading input
  			//-------------------------------------------------------------
  			
  			MapType exp;
  			MzDataFile f;
  			f.setLogType(log_type_);
  			f.getOptions().setRTRange(DRange<1>(rt_l,rt_u));
  			f.getOptions().setMZRange(DRange<1>(mz_l,mz_u));
  			f.getOptions().setIntensityRange(DRange<1>(it_l,it_u));
  			f.load(in,exp);						
  		
  			//-------------------------------------------------------------
  			// calculations
  			//-------------------------------------------------------------
  			
  			//remove ms level first (might be a lot of spectra)
  			exp.erase(remove_if(exp.begin(), exp.end(), InMSLevelRange<MapType::SpectrumType>(levels, true)), exp.end());
  			
  			//remove by scan mode (might be a lot of spectra)
  			bool rem_mode = setByUser_("remove_mode");
  			writeDebug_(String("Remove by mode: ") + String(rem_mode),3);
  			if (rem_mode)
  			{
  				String mode = getStringOption_("remove_mode");
  				writeDebug_(String("Removing mode: ") + mode,3);
  				for (UInt i=0; i<InstrumentSettings::SIZE_OF_SCANMODE; ++i)
  				{
  					if (InstrumentSettings::NamesOfScanMode[i]==mode)
  					{
  						exp.erase(remove_if(exp.begin(), exp.end(), HasScanMode<MapType::SpectrumType>((InstrumentSettings::ScanMode)i)), exp.end());
  					}
  				}
  			}
  				
  			//remove empty scans
  			exp.erase(remove_if(exp.begin(), exp.end(), IsEmptySpectrum<MapType::SpectrumType>()), exp.end());
 				
 				//sort by RT and m/z
   			bool sort = getFlag_("sort");
  			writeDebug_(String("Sorting output data: ") + String(sort),3);
  			if (sort)
  			{
  				//if meta data arrays are present, remove them and warn
  				if (exp.clearMetaDataArrays())
  				{
  					writeLog_("Warning: Spectrum meta data arrays cannot be sorted. They are deleted.");
  				}
  				
  				//sort
  				exp.sortSpectra(true);
  			}

				// calculate S/N values and write them instead
				if (sn > 0)
				{
					SignalToNoiseEstimatorMedian < MapType::SpectrumType > snm;
					Param const& dc_param = getParam_().copy("algorithm:SignalToNoise:",true);
					snm.setParameters(dc_param);
					for(MapType::Iterator it = exp.begin(); it != exp.end(); ++it)
					{
						snm.init(it->begin(), it->end());
						for (MapType::SpectrumType::Iterator spec = it->begin(); spec != it->end(); ++spec)
						{
							std::cout << "sn is: " << snm.getSignalToNoise(spec) << "\n";
							if (snm.getSignalToNoise(spec) < sn) spec->setIntensity(0);
						}
						it->erase(remove_if(it->begin(), it->end(), InIntensityRange<MapType::PeakType>(1,numeric_limits<MapType::PeakType::IntensityType>::max(), true)) , it->end());
					}
				}
 
  			//-------------------------------------------------------------
  			// writing output
  			//-------------------------------------------------------------
  			
  			f.store(out,exp);
      }
      else if (out_type == FileHandler::FEATUREXML)
      {
        //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------

        typedef FeatureMap<> FeatureMapType;
        FeatureMapType feature_map;
        FeatureXMLFile f;
        //f.setLogType(log_type_);
        // this does not work yet implicitly - not supported by FeatureXMLFile
        f.getOptions().setRTRange(DRange<1>(rt_l,rt_u));
        f.getOptions().setMZRange(DRange<1>(mz_l,mz_u));
        f.getOptions().setIntensityRange(DRange<1>(it_l,it_u));
        f.load(in,feature_map);                 

      
        //-------------------------------------------------------------
        // calculations
        //-------------------------------------------------------------
        typedef FeatureMapType::FeatureType FeatureType;
        
        //copy all properties
        FeatureMapType map_sm = feature_map;
        //.. but delete feature information 
        map_sm.clear();
        
        bool rt_ok, mz_ok, int_ok, charge_ok, q_ok;
        
        // only keep charge ch_l:ch_u   (WARNING: featurefiles without charge information have charge=0, see Ctor of KERNEL/Feature.h)
        for (uint i = 0; i<feature_map.size(); ++i)
        {
          if (f.getOptions().getRTRange().encloses(DPosition<1>(feature_map[i].getRT()))) { rt_ok = true; } else {rt_ok = false;}
          if (f.getOptions().getMZRange().encloses(DPosition<1>(feature_map[i].getMZ()))) { mz_ok = true; } else {mz_ok = false;}
          if (f.getOptions().getIntensityRange().encloses(DPosition<1>(feature_map[i].getIntensity()))) { int_ok = true; } else {int_ok = false;}
          if ((charge_l <= feature_map[i].getCharge()) && (feature_map[i].getCharge() <= charge_u)) { charge_ok = true; } else {charge_ok = false;}
          if ((q_l <= feature_map[i].getOverallQuality()) && (feature_map[i].getOverallQuality() <= q_u)) { q_ok = true; } else {q_ok = false;}
          
          //std::cout << feature_map[i].getRT() << " " << feature_map[i].getMZ() << " " << feature_map[i].getIntensity() << " " << feature_map[i].getCharge() << " "<< feature_map[i].getOverallQuality() << " "; 
          if (rt_ok == true && mz_ok == true && int_ok == true && charge_ok == true && q_ok == true)
          {
            //std::cout << rt_ok << mz_ok << int_ok << charge_ok << "\n";
            map_sm.push_back (feature_map[i]);              
          }//else {std::cout << "\n";}
        }
        map_sm.updateRanges();
        
        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
                
        f.store(out,map_sm);
      }
      else if (out_type == FileHandler::CONSENSUSXML)
      {
        //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------

        ConsensusMap consensus_map;
        ConsensusXMLFile f;
        //f.setLogType(log_type_);
        // this does not work yet implicitly - not supported by FeatureXMLFile
        f.getOptions().setRTRange(DRange<1>(rt_l,rt_u));
        f.getOptions().setMZRange(DRange<1>(mz_l,mz_u));
        f.getOptions().setIntensityRange(DRange<1>(it_l,it_u));
        f.load(in,consensus_map);                 

      
        //-------------------------------------------------------------
        // calculations
        //-------------------------------------------------------------
        consensus_map.updateRanges();
        
        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
                
        f.store(out,consensus_map);
      }
      else
      {
        writeLog_("Unknown input file type given. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;          
      }
			
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPFileFilter tool;
	return tool.main(argc,argv);
}

/// @endcond
