// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm, Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FileFilter FileFilter

	@brief Extracts portions of the data from an mzML, featureXML or consensusXML file.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FileFilter \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool yielding output @n in mzML, featureXML @n or consensusXML format</td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool that profits on reduced input </td>
		</tr>

	</table>
</CENTER>

	With this tool it is possible to extract m/z, retention time and intensity ranges from an input file
	and to write all data that lies within the given ranges to an output file.

	Depending on the input file type, additional specific operations are possible:
	- mzML
		- extract spectra of a certain MS level
		- filter by signal-to-noise estimation
		- filter by scan mode of the spectra
	- featureXML
		- filter by feature charge
		- filter by feature size (number of subordinate features)
		- filter by overall feature quality
	- consensusXML
		- filter by size (number of elements in consensus features)
		- filter by consensus feature charge
		- filter by map (extracts specified maps and re-evaluates consensus centroid)@n e.g. FileFilter -map 2 3 5 -in file1.consensusXML -out file2.consensusXML@n If a single map is specified, the feature itself can be extracted.@n e.g. FileFilter -map 5 -in file1.consensusXML -out file2.featureXML


	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FileFilter.cli

	For the parameters of the S/N algorithm section see the class documentation there: @n
		@ref OpenMS::SignalToNoiseEstimatorMedian "sn"@n

	@todo add tests for selecting modes (port remove modes) (Andreas)
	@improvement MS2 and higher spectra should be filtered according to precursor m/z and RT. The MzMLFile, MzDataFile, MzXMLFile have to be changed for that (Hiwi)
               Currently when specifying mz or RT filters, they will also be applied to MS levels >=2 (not really what you usually want). To work around this,
               you need to extract the MS2 levels beforehand, do the filtering on MS1 and merge them back together.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileFilter
	: public TOPPBase
{
	public:

	TOPPFileFilter()
		: TOPPBase("FileFilter","Extracts or manipulates portions of data from peak, feature or consensus-feature files.")
	{
	}

	protected:

	typedef MSExperiment<Peak1D> MapType;

	void registerOptionsAndFlags_()
	{
		String formats("mzML,featureXML,consensusXML");

		registerInputFile_("in","<file>","","input file ");
		setValidFormats_("in",StringList::create(formats));

		registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
		setValidStrings_("in_type",StringList::create(formats));

		registerOutputFile_("out","<file>","","output file");
		setValidFormats_("out",StringList::create(formats));
		
		registerStringOption_("out_type", "<type>", "", "output file type -- default: determined from file extension or content\n", false);
		setValidStrings_("out_type",StringList::create(formats));

		registerStringOption_("mz","[min]:[max]",":","m/z range to extract", false);
		registerStringOption_("rt","[min]:[max]",":","retention time range to extract", false);
		registerStringOption_("int","[min]:[max]",":","intensity range to extract", false);

		registerFlag_("sort","sorts the output according to RT and m/z.");

		addText_("peak data options:");
		registerDoubleOption_("sn", "<s/n ratio>", 0, "write peaks with S/N > 'sn' values only", false);
		registerIntList_("level","\"i,j,...\"",IntList::create("1,2,3"),"MS levels to extract", false);
		registerIntList_("map","\"i,j,...\"",IntList::create(""),"maps to be extracted from a consensus", false);
		registerFlag_("sort_peaks","sorts the peaks according to m/z.");
		registerFlag_("no_chromatograms", "No conversion to space-saving real chromatograms, e.g. from SRM scans.");
		registerFlag_("remove_chromatograms", "Removes chromatograms stored in a file.");
		registerFlag_("map_and", "AND connective of map selection instead of OR.");

		addEmptyLine_();
		addText_("Remove spectra: ");
		registerFlag_("remove_zoom","Remove zoom (enhanced resolution) scans");

		registerStringOption_("remove_mode","<mode>","","Remove scans by scan mode\n",false);
		StringList mode_list;
		for (Size i=0; i<InstrumentSettings::SIZE_OF_SCANMODE; ++i)
		{
			mode_list.push_back(InstrumentSettings::NamesOfScanMode[i]);
		}
		setValidStrings_("remove_mode",mode_list);
		addEmptyLine_();
		registerStringOption_("remove_activation","<activation>","","Remove MSn scans where any of its precursors features a certain activation method\n",false);
		StringList activation_list;
		for (Size i=0; i<Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
		{
			activation_list.push_back(Precursor::NamesOfActivationMethod[i]);
		}
		setValidStrings_("remove_activation",activation_list);
		addEmptyLine_();

		addText_("Select spectra (remove all others):");
		registerFlag_("select_zoom", "Select zoom (enhanced resolution) scans");
		registerStringOption_("select_mode", "<mode>", "", "Selects scans by scan mode\n", false);
		setValidStrings_("select_mode", mode_list);
		registerStringOption_("select_activation", "<activation>", "", "Select MSn scans where any of its percursors features a certain activation method\n", false);
		setValidStrings_("select_activation", activation_list);
		addEmptyLine_();


		addText_("feature data options:");
		registerStringOption_("charge","[min]:[max]",":","charge range to extract", false);
		registerStringOption_("size","[min]:[max]",":","size range to extract", false);
		registerStringOption_("q","[min]:[max]",":","OverallQuality range to extract [0:1]", false);

		addText_("consensus feature data options:");
		registerStringOption_("size","[min]:[max]",":","size range to extract", false);
		registerStringOption_("charge","[min]:[max]",":","charge range to extract", false);

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

		//input file name and type
		String in = getStringOption_("in");
		FileHandler fh;

		FileTypes::Type in_type = fh.nameToType(getStringOption_("in_type"));
		
		if (in_type==FileTypes::UNKNOWN)
		{
			in_type = fh.getType(in);
			writeDebug_(String("Input file type: ") + fh.typeToName(in_type), 2);
		}
		
		if (in_type==FileTypes::UNKNOWN)
		{
			writeLog_("Error: Could not determine input file type!");
			return PARSE_ERROR;
		}
		
		
		//output file name and type
		String out = getStringOption_("out");

		FileTypes::Type out_type = fh.nameToType(getStringOption_("out_type"));
		
		if (out_type==FileTypes::UNKNOWN)
		{
			out_type = fh.getTypeByFileName(out);
			writeDebug_(String("Output file type: ") + fh.typeToName(out_type), 2);
		}
		
		if (out_type==FileTypes::UNKNOWN)
		{
			writeLog_("Error: Could not determine output file type!");
			return PARSE_ERROR;
		}
		
		
		bool no_chromatograms(getFlag_("no_chromatograms"));

		//ranges
		String mz, rt, it, charge, size, q;
		IntList levels, maps;
		double mz_l, mz_u, rt_l, rt_u, it_l, it_u, sn, charge_l, charge_u, size_l, size_u, q_l, q_u;
		//initialize ranges
		mz_l = rt_l = it_l = charge_l = size_l = q_l = -1 * numeric_limits<double>::max();
		mz_u = rt_u = it_u = charge_u = size_u = q_u = numeric_limits<double>::max();

		rt = getStringOption_("rt");
		mz = getStringOption_("mz");
		it = getStringOption_("int");
		levels = getIntList_("level");
		maps = getIntList_("map");
		sn = getDoubleOption_("sn");
		charge = getStringOption_("charge");
		size = getStringOption_("size");
		q = getStringOption_("q");

		//convert bounds to numbers
		try
		{
			//rt
			parseRange_(rt,rt_l,rt_u);
			//mz
			parseRange_(mz,mz_l,mz_u);
			//int
			parseRange_(it,it_l,it_u);
			//charge (features only)
			parseRange_(charge,charge_l,charge_u);
			//size (features and consensus features only)
			parseRange_(size,size_l,size_u);
			//overall quality (features only)
			parseRange_(q,q_l,q_u);
		}
		catch(Exception::ConversionError&)
		{
			String tmp;
			for(IntList::iterator it = levels.begin(); it != levels.end();++it)
			{
				tmp += *it;
			}

			writeLog_("Invalid boundary '" + tmp + "' given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}

		//sort by RT and m/z
 		bool sort = getFlag_("sort");
		writeDebug_("Sorting output data: " + String(sort),3);

		if (in_type == FileTypes::MZML)
		{
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------

  			MapType exp;
  			MzMLFile f;
  			f.setLogType(log_type_);
  			f.getOptions().setRTRange(DRange<1>(rt_l,rt_u));
  			f.getOptions().setMZRange(DRange<1>(mz_l,mz_u));
  			f.getOptions().setIntensityRange(DRange<1>(it_l,it_u));
  			f.load(in,exp);

			if (!no_chromatograms)
			{
				// convert the spectra chromatograms to real chromatograms
				ChromatogramTools chrom_tools;
				chrom_tools.convertSpectraToChromatograms(exp, true);
			}

			bool remove_chromatograms(getFlag_("remove_chromatograms"));
			if (remove_chromatograms)
			{
				exp.setChromatograms(vector<MSChromatogram<> >());
			}

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------

			//remove ms level first (might be a lot of spectra)
			exp.erase(remove_if(exp.begin(), exp.end(), InMSLevelRange<MapType::SpectrumType>(levels, true)), exp.end());

			//remove by scan mode (might be a lot of spectra)
			String remove_mode = getStringOption_("remove_mode");
			if (!remove_mode.empty())
			{
				writeDebug_("Removing mode: " + remove_mode,3);
				for (Size i=0; i<InstrumentSettings::SIZE_OF_SCANMODE; ++i)
				{
					if (InstrumentSettings::NamesOfScanMode[i]==remove_mode)
					{
						exp.erase(remove_if(exp.begin(), exp.end(), HasScanMode<MapType::SpectrumType>((InstrumentSettings::ScanMode)i)), exp.end());
					}
				}
			}

			//select by scan mode (might be a lot of spectra)
			String select_mode = getStringOption_("select_mode");
			if (!select_mode.empty())
			{
				writeDebug_("Selecting mode: " + select_mode,3);
				for (Size i=0; i<InstrumentSettings::SIZE_OF_SCANMODE; ++i)
				{
					if (InstrumentSettings::NamesOfScanMode[i]==select_mode)
					{
					exp.erase(remove_if(exp.begin(), exp.end(), HasScanMode<MapType::SpectrumType>((InstrumentSettings::ScanMode)i, true)), exp.end());
					}
				}
			}



			//remove by activation mode (might be a lot of spectra)
			String remove_activation = getStringOption_("remove_activation");
			if (!remove_activation.empty())
			{
				writeDebug_("Removing scans with activation mode: " + remove_activation,3);
				for (Size i=0; i<Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
				{
					if (Precursor::NamesOfActivationMethod[i]==remove_activation)
					{
						exp.erase(remove_if(exp.begin(), exp.end(), HasActivationMethod<MapType::SpectrumType>(StringList::create(remove_activation))), exp.end());
					}
				}
			}

			//select by activation mode
			String select_activation = getStringOption_("select_activation");
			if (!select_activation.empty())
			{
				writeDebug_("Selecting scans with activation mode: " + select_activation, 3);
				for (Size i = 0; i < Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
				{
					if (Precursor::NamesOfActivationMethod[i] == select_activation)
					{
						exp.erase(remove_if(exp.begin(), exp.end(), HasActivationMethod<MapType::SpectrumType>(StringList::create(select_activation),true)), exp.end());
					}
				}
			}

			//remove zoom scans (might be a lot of spectra)
			if (getFlag_("remove_zoom"))
			{
				writeDebug_("Removing zoom scans",3);
				exp.erase(remove_if(exp.begin(), exp.end(),IsZoomSpectrum<MapType::SpectrumType>()), exp.end());
			}

			if (getFlag_("select_zoom"))
			{
				writeDebug_("Selecting zoom scans", 3);
				exp.erase(remove_if(exp.begin(), exp.end(), IsZoomSpectrum<MapType::SpectrumType>(true)), exp.end());
			}

 			//remove empty scans
 			exp.erase(remove_if(exp.begin(), exp.end(), IsEmptySpectrum<MapType::SpectrumType>()), exp.end());

      //sort
      if (sort)
      {
        exp.sortSpectra(true);
      }
			if (getFlag_("sort_peaks"))
			{
				for (Size i=0; i<exp.size(); ++i)
				{
					exp[i].sortByPosition();
				}
			}

			// calculate S/N values and delete datapoints below S/N threshold
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
						if (snm.getSignalToNoise(spec) < sn) spec->setIntensity(0);
					}
					it->erase(remove_if(it->begin(), it->end(), InIntensityRange<MapType::PeakType>(1,numeric_limits<MapType::PeakType::IntensityType>::max(), true)) , it->end());
				}
			}

  			//-------------------------------------------------------------
  			// writing output
  			//-------------------------------------------------------------

			//annotate output with data processing info
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FILTERING));
			f.store(out,exp);
		}
		else if (in_type == FileTypes::FEATUREXML)
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
			map_sm.clear(false);

			bool rt_ok, mz_ok, int_ok, charge_ok, size_ok, q_ok;

			// only keep charge ch_l:ch_u   (WARNING: featurefiles without charge information have charge=0, see Ctor of KERNEL/Feature.h)
			for (uint i = 0; i<feature_map.size(); ++i)
			{
				if (f.getOptions().getRTRange().encloses(DPosition<1>(feature_map[i].getRT()))) { rt_ok = true; } else {rt_ok = false;}
				if (f.getOptions().getMZRange().encloses(DPosition<1>(feature_map[i].getMZ()))) { mz_ok = true; } else {mz_ok = false;}
				if (f.getOptions().getIntensityRange().encloses(DPosition<1>(feature_map[i].getIntensity()))) { int_ok = true; } else {int_ok = false;}
				if ((charge_l <= feature_map[i].getCharge()) && (feature_map[i].getCharge() <= charge_u)) { charge_ok = true; } else {charge_ok = false;}
				if ((size_l <= feature_map[i].getSubordinates().size()) && (feature_map[i].getSubordinates().size() <= size_u)) { size_ok = true; } else { size_ok = false;}
				if ((q_l <= feature_map[i].getOverallQuality()) && (feature_map[i].getOverallQuality() <= q_u)) { q_ok = true; } else {q_ok = false;}

				//std::cout << feature_map[i].getRT() << " " << feature_map[i].getMZ() << " " << feature_map[i].getIntensity() << " " << feature_map[i].getCharge() << " "<< feature_map[i].getOverallQuality() << " ";
				if (rt_ok == true && mz_ok == true && int_ok == true && charge_ok == true && size_ok == true && q_ok == true)
				{
					//std::cout << rt_ok << mz_ok << int_ok << charge_ok << size_ok << q_ok << "\n";
					map_sm.push_back (feature_map[i]);
				}//else {std::cout << "\n";}
			}
			map_sm.updateRanges();

			// sort if desired
      if (sort)
      {
        map_sm.sortByPosition();
      }

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------

			//annotate output with data processing info
			addDataProcessing_(map_sm, getProcessingInfo_(DataProcessing::FILTERING));

			f.store(out,map_sm);
		}
		else if (in_type == FileTypes::CONSENSUSXML)
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
			
			// copy all properties
			ConsensusMap consensus_map_filtered = consensus_map;
			//.. but delete feature information
			consensus_map_filtered.resize(0);
			
			bool charge_ok, size_ok;
			
			for ( ConsensusMap::const_iterator citer = consensus_map.begin(); citer != consensus_map.end(); ++citer)
			{
				if ((charge_l <= citer->getCharge()) && (citer->getCharge() <= charge_u)) { charge_ok = true; } else {charge_ok = false;}
				if ((citer->size() >= size_l) && (citer->size() <= size_u)) { size_ok = true; } else { size_ok = false;}
				
				if (charge_ok == true && size_ok == true)
				{
					consensus_map_filtered.push_back(*citer);
				}
			}
			
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			consensus_map_filtered.updateRanges();
			
			// sort if desired
			if (sort)
			{
				consensus_map_filtered.sortByPosition();
			}
			
			if (out_type == FileTypes::FEATUREXML)
			{				
				if (maps.size() == 1) // When extracting a feature map from a consensus map, only one map ID should be specified. Hence 'maps' should contain only one interger.
				{
					FeatureMap<> feature_map_filtered;
					FeatureXMLFile ff;
					
					for (ConsensusMap::Iterator cm_it=consensus_map_filtered.begin(); cm_it!=consensus_map_filtered.end(); ++cm_it)
					{
						
						for(ConsensusFeature::HandleSetType::const_iterator fh_iter = (*cm_it).getFeatures().begin();
							fh_iter != (*cm_it).getFeatures().end();
							++fh_iter)
						{
							if ((int)(*fh_iter).getMapIndex() == maps[0])
							{
								Feature feature;
								feature.setRT((*fh_iter).getRT());
								feature.setMZ((*fh_iter).getMZ());
								feature.setIntensity((*fh_iter).getIntensity());
								feature.setCharge((*fh_iter).getCharge());
								feature_map_filtered.push_back(feature);
							}
						}
					}
					
					//-------------------------------------------------------------
					// writing output
					//-------------------------------------------------------------
					
					//annotate output with data processing info
					addDataProcessing_(feature_map_filtered, getProcessingInfo_(DataProcessing::FILTERING));
					
					feature_map_filtered.applyMemberFunction(&UniqueIdInterface::setUniqueId);
					
					ff.store(out,feature_map_filtered);
				}
				else {
					writeLog_("When extracting a feature map from a consensus map, only one map ID should be specified. The 'map' parameter contains more than one. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}
			else if (out_type == FileTypes::CONSENSUSXML)
			{
				// generate new consensuses with features that appear in the 'maps' list        
				ConsensusMap cm_new; // new consensus map
				
				for (IntList::iterator map_it=maps.begin();map_it!=maps.end();++map_it)
				{
					cm_new.getFileDescriptions()[*map_it].filename = consensus_map_filtered.getFileDescriptions()[*map_it].filename;
					cm_new.getFileDescriptions()[*map_it].size = consensus_map_filtered.getFileDescriptions()[*map_it].size;
					cm_new.getFileDescriptions()[*map_it].unique_id = consensus_map_filtered.getFileDescriptions()[*map_it].unique_id;
				}

				for (ConsensusMap::Iterator cm_it = consensus_map_filtered.begin(); cm_it != consensus_map_filtered.end(); ++cm_it) // iterate over consensuses in the original consensus map
				{					
					ConsensusFeature consensus_feature_new(*cm_it); // new consensus feature
					consensus_feature_new.clear();

					ConsensusFeature::HandleSetType::const_iterator fh_it = cm_it->getFeatures().begin();
					ConsensusFeature::HandleSetType::const_iterator fh_it_end = cm_it->getFeatures().end();
					for(; fh_it != fh_it_end; ++fh_it) // iterate over features in consensus
					{
						if (maps.contains(fh_it->getMapIndex()))
						{
							consensus_feature_new.insert(*fh_it);
						}												
					}
					
					consensus_feature_new.computeConsensus(); // evaluate position of the consensus
					bool and_connective=getFlag_("map_and");

					if ((consensus_feature_new.size() != 0 && !and_connective) || (consensus_feature_new.size()==maps.size() && and_connective)) // add the consensus to the consensus map only if it is non-empty 
					{            
						cm_new.push_back(consensus_feature_new);
					}
				}				
								
				// assign unique ids
				cm_new.applyMemberFunction(&UniqueIdInterface::setUniqueId);
			
				//-------------------------------------------------------------
				// writing output
				//-------------------------------------------------------------
				
				if (maps.empty())
				{
					//annotate output with data processing info
					addDataProcessing_(consensus_map_filtered, getProcessingInfo_(DataProcessing::FILTERING));
					
					f.store(out, consensus_map_filtered);
				}
				else
				{
					//annotate output with data processing info
					addDataProcessing_(cm_new, getProcessingInfo_(DataProcessing::FILTERING));
					
					f.store(out, cm_new);
				}
			}
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
