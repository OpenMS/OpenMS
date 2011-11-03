// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Lars Nilse, Chris Bielow, Hendrik Brauer $
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
	- featureXML / consensusXML:
    - remove items with a certain meta value annotation. Allowing for >, < and = comparisons. List types are compared by length, not content. Integer, Double and String are compared using their build-in operators.
		- filter sequences, e.g. "LYSNLVER" or the modification "(Phospho)"@n e.g. FileFilter -id:sequences_whitelist Phospho -in file1.consensusXML -out file2.consensusXML
		- filter accessions, e.g. "sp|P02662|CASA1_BOVIN"
		- remove features with annotations
		- remove features without annotations
		- remove unassigned peptide identifications
		- filter id with best score of features with multiple peptide identifications@n e.g. FileFilter -id:remove_unannotated_features -id:remove_unassigned_ids -id:keep_best_score_id -in file1.featureXML -out file2.featureXML
		- remove features with id clashes (different sequences mapped to one feature)

	The Priority of the id-flags is (decreasing order): remove_annotated_features / remove_unannotated_features -> remove_clashes -> keep_best_score_id -> sequences_whitelist / accessions_whitelist

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

	private:
	static bool checkPeptideIdentification_(BaseFeature& feature, const bool remove_annotated_features, const bool remove_unannotated_features, const StringList& sequences, const StringList& accessions, const bool keep_best_score_id, const bool remove_clashes)
	{
		//flag: remove_annotated_features and non-empty peptideIdentifications
		if (remove_annotated_features && !feature.getPeptideIdentifications().empty())
		{
			return false;
		}
		//flag: remove_unannotated_features and no peptideIdentifications
		if (remove_unannotated_features && feature.getPeptideIdentifications().empty())
		{
			return false;
		}
		//flag: remove_clashes
		if (remove_clashes && !feature.getPeptideIdentifications().empty())
		{
			String temp = feature.getPeptideIdentifications().begin()->getHits().begin()->getSequence().toString();
			//loop over all peptideIdentifications
			for (vector<PeptideIdentification>::const_iterator pep_id_it = feature.getPeptideIdentifications().begin(); pep_id_it != feature.getPeptideIdentifications().end(); ++pep_id_it)
			{
				//loop over all peptideHits
				for (vector<PeptideHit>::const_iterator pep_hit_it = pep_id_it->getHits().begin(); pep_hit_it != pep_id_it->getHits().end(); ++pep_hit_it)
				{
					if (pep_hit_it->getSequence().toString()!=temp)
					{
						return false;
					}
				}
			}
		}
		//flag: keep_best_score_id
		if (keep_best_score_id && !feature.getPeptideIdentifications().empty())
		{
			PeptideIdentification temp = feature.getPeptideIdentifications().front();
			//loop over all peptideIdentifications
			for (vector<PeptideIdentification>::const_iterator pep_id_it = feature.getPeptideIdentifications().begin(); pep_id_it != feature.getPeptideIdentifications().end(); ++pep_id_it)
			{
				//loop over all peptideHits
				for (vector<PeptideHit>::const_iterator pep_hit_it = pep_id_it->getHits().begin(); pep_hit_it != pep_id_it->getHits().end(); ++pep_hit_it)
				{
					if ((pep_id_it->isHigherScoreBetter() && pep_hit_it->getScore() > temp.getHits().front().getScore()) ||
						(!pep_id_it->isHigherScoreBetter() && pep_hit_it->getScore() < temp.getHits().front().getScore()))
					{
						temp = *pep_id_it;
					}
				}
			}
			feature.setPeptideIdentifications(vector<PeptideIdentification> (1,temp));
			// not filtering sequences or accessions
      if (sequences.size() == 0 && accessions.size() == 0)
			{
				return true;
			}
		}
		//flag: sequences or accessions
		if (sequences.size() > 0 || accessions.size() > 0)
		{
			bool sequen = false;
			bool access = false;
			//loop over all peptideIdentifications
			for (vector<PeptideIdentification>::const_iterator pep_id_it = feature.getPeptideIdentifications().begin(); pep_id_it != feature.getPeptideIdentifications().end(); ++pep_id_it)
			{
				//loop over all peptideHits
				for (vector<PeptideHit>::const_iterator pep_hit_it = pep_id_it->getHits().begin(); pep_hit_it != pep_id_it->getHits().end(); ++pep_hit_it)
				{
					//loop over all sequence entries of the StringList
					for (StringList::ConstIterator seq_it = sequences.begin(); seq_it != sequences.end(); ++seq_it)
					{
						if (pep_hit_it->getSequence().toString().hasSubstring(*seq_it)
							|| pep_hit_it->getSequence().toUnmodifiedString().hasSubstring(*seq_it) )
						{
							sequen = true;
						}
					}
					//loop over all accessions of the peptideHits
					for (vector<String>::const_iterator p_acc_it = pep_hit_it->getProteinAccessions().begin(); p_acc_it != pep_hit_it->getProteinAccessions().end(); ++p_acc_it)
					{
						//loop over all accessions entries of the StringList
						for (StringList::ConstIterator acc_it = accessions.begin(); acc_it != accessions.end(); ++acc_it)
						{
							if (p_acc_it->hasSubstring(*acc_it))
							{
								access = true;
							}
						}
					}
				}
			}
			if (sequences.size() > 0 && accessions.size() > 0)
			{
				return (sequen && access);
			}
			if (sequences.size() > 0)
			{
				return sequen;
			}else
			{
				return access;
			}
		}
		return true;
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

		addEmptyLine_();
		addText_("Peak data options:");
		registerDoubleOption_("sn", "<s/n ratio>", 0, "write peaks with S/N > 'sn' values only", false);
    registerIntList_("rm_pc_charge","i j ...", IntList(), "Remove MS(2) spectra with these precursor charges. All spectra without precursor are kept!", false);
		registerIntList_("level","i j ...", IntList::create("1,2,3"),"MS levels to extract", false);
		registerFlag_("sort_peaks","sorts the peaks according to m/z.");
		registerFlag_("no_chromatograms", "No conversion to space-saving real chromatograms, e.g. from SRM scans.");
		registerFlag_("remove_chromatograms", "Removes chromatograms stored in a file.");

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
		registerStringOption_("select_activation", "<activation>", "", "Select MSn scans where any of its precursors features a certain activation method\n", false);
		setValidStrings_("select_activation", activation_list);
		addEmptyLine_();

    addText_("Feature&Consensus data options:");
    registerStringOption_("charge","[min]:[max]",":","charge range to extract", false);
    registerStringOption_("size","[min]:[max]",":","size range to extract", false);
    registerStringList_("remove_meta", "<name> 'lt|eq|gt' <value>", StringList(), "Expects a 3-tuple, with name of meta value, the comparison operator (equal, less or greater) and the value to compare to. All comparisons are done after converting the given value to the corresponding data value type of the meta value!", false);


		addText_("Feature data options:");
		registerStringOption_("q","[min]:[max]",":","Overall quality range to extract [0:1]", false);

		addText_("Consensus feature data options:");
    registerIntList_("map","i j ...",IntList::create(""),"maps to be extracted from a consensus", false);
    registerFlag_("map_and", "AND connective of map selection instead of OR.");

		addEmptyLine_();
		registerTOPPSubsection_("id","id section");
		addText_("The Priority of the id-flags is: remove_annotated_features / remove_unannotated_features -> remove_clashes -> keep_best_score_id -> sequences_whitelist / accessions_whitelist");
		registerFlag_("id:remove_clashes", "remove features with id clashes (different sequences mapped to one feature)", true);
		registerFlag_("id:keep_best_score_id", "in case of multiple peptide identifications, keep only the id with best score");
		registerStringList_("id:sequences_whitelist", "<sequence>", StringList(), "keep only features with white listed sequences, e.g. LYSNLVER or the modification (Oxidation)", false);
		registerStringList_("id:accessions_whitelist", "<accessions>", StringList(), "keep only features with white listed accessions, e.g. sp|P02662|CASA1_BOVIN", false);
		registerFlag_("id:remove_annotated_features", "remove features with annotations");
		registerFlag_("id:remove_unannotated_features", "remove features without annotations");
		registerFlag_("id:remove_unassigned_ids", "remove unassigned peptide identifications");

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

  bool checkMetaOk(const MetaInfoInterface& mi, const StringList& meta_info)
  {
    if (!mi.metaValueExists(meta_info[0])) return true; // not having the meta value means passing the test

    DataValue v_data = mi.getMetaValue(meta_info[0]);
    DataValue v_user;
    if      (v_data.valueType() == DataValue::STRING_VALUE) v_user = String (meta_info[2]);
    else if (v_data.valueType() == DataValue::INT_VALUE) v_user = String (meta_info[2]).toInt();
    else if (v_data.valueType() == DataValue::DOUBLE_VALUE) v_user = String (meta_info[2]).toDouble();
    else if (v_data.valueType() == DataValue::STRING_LIST) v_user = StringList::create(meta_info[2]);
    else if (v_data.valueType() == DataValue::INT_LIST) v_user = IntList::create(meta_info[2]);
    else if (v_data.valueType() == DataValue::DOUBLE_LIST) v_user = DoubleList::create(meta_info[2]);
    else if (v_data.valueType() == DataValue::EMPTY_VALUE) v_user = DataValue::EMPTY;

    if (meta_info[1]=="lt")
    {
      return v_data < v_user;
    }
    else if (meta_info[1]=="eq")
    {
      return v_data == v_user;
    }
    else if (meta_info[1]=="gt")
    {
      return v_data > v_user;
    }
    else
    {
      writeLog_("Internal Error. Meta value filtering got invalid comparison operator ('" + meta_info[1] + "'), which should have been catched before! Aborting!");
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Illegal meta value filtering operator!");
    }
  }

	ExitCodes main_(int , const char**)
	{

		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------

		//input file name and type
		String in = getStringOption_("in");
		FileHandler fh;

		FileTypes::Type in_type = fh.getType(in);
		//only use flag in_type, if the in_type cannot be determined by file
		if (in_type==FileTypes::UNKNOWN)
		{
			in_type = fh.nameToType(getStringOption_("in_type"));
			writeDebug_(String("Input file type: ") + fh.typeToName(in_type), 2);
		}
	
		//output file name and type
		String out = getStringOption_("out");

		FileTypes::Type out_type = fh.getTypeByFileName(out);

		//only use flag out_type, if the out_type cannot be determined by file
		if (out_type==FileTypes::UNKNOWN)
		{
			out_type = fh.nameToType(getStringOption_("out_type"));
			writeDebug_(String("Output file type: ") + fh.typeToName(out_type), 2);
		}
		//use in_type as out_type, if out_type cannot be determined by file or out_type flag
		if (out_type==FileTypes::UNKNOWN)
		{
			out_type = in_type;
			writeDebug_(String("Output file type: ") + fh.typeToName(out_type), 2);
		}

		bool no_chromatograms(getFlag_("no_chromatograms"));

		//ranges
		String mz, rt, it, charge, size, q;
		double mz_l, mz_u, rt_l, rt_u, it_l, it_u, sn, charge_l, charge_u, size_l, size_u, q_l, q_u;
		//initialize ranges
		mz_l = rt_l = it_l = charge_l = size_l = q_l = -1 * numeric_limits<double>::max();
		mz_u = rt_u = it_u = charge_u = size_u = q_u = numeric_limits<double>::max();

		rt = getStringOption_("rt");
		mz = getStringOption_("mz");
		it = getStringOption_("int");
		IntList levels = getIntList_("level");
		IntList maps = getIntList_("map");
		sn = getDoubleOption_("sn");
		charge = getStringOption_("charge");
		size = getStringOption_("size");
		q = getStringOption_("q");

		//id-filtering parameters
		bool remove_annotated_features = getFlag_("id:remove_annotated_features");
		bool remove_unannotated_features = getFlag_("id:remove_unannotated_features");
		bool remove_unassigned_ids = getFlag_("id:remove_unassigned_ids");
		StringList sequences = getStringList_("id:sequences_whitelist");
		StringList accessions = getStringList_("id:accessions_whitelist");
		bool keep_best_score_id = getFlag_("id:keep_best_score_id");
		bool remove_clashes = getFlag_("id:remove_clashes");

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

      //remove forbidden precursor charges
      IntList rm_pc_charge = getIntList_("rm_pc_charge");
      if (rm_pc_charge.size()>0) exp.erase(remove_if(exp.begin(), exp.end(), HasPrecursorCharge<MapType::SpectrumType>(rm_pc_charge, false)), exp.end());

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

			// calculate S/N values and delete data points below S/N threshold
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

		else if (in_type == FileTypes::FEATUREXML || in_type == FileTypes::CONSENSUSXML)
    {
      bool meta_ok = true; // assume true by default (as meta might not be checked below)
      StringList meta_info = getStringList_("remove_meta");
      bool check_meta = (meta_info.size()>0);
      if (check_meta && meta_info.size()!=3)
      {
        writeLog_("Param 'remove_meta' has invalid number of arguments. Expected 3, got " + String(meta_info.size()) + ". Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
      if (check_meta && !(meta_info[1]=="lt" || meta_info[1]=="eq" || meta_info[1]=="gt"))
      {
        writeLog_("Param 'remove_meta' has invalid second argument. Expected one of 'lt', 'eq' or 'gt'. Got '" + meta_info[1] + "'. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }

		  if (in_type == FileTypes::FEATUREXML)
      {
			  //-------------------------------------------------------------
			  // loading input
			  //-------------------------------------------------------------

			  FeatureMap<> feature_map;
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

			  //copy all properties
			  FeatureMap<> map_sm = feature_map;
			  //.. but delete feature information
			  map_sm.clear(false);

			  bool rt_ok, mz_ok, int_ok, charge_ok, size_ok, q_ok, annotation_ok;

			  // only keep charge ch_l:ch_u   (WARNING: feature files without charge information have charge=0, see Ctor of KERNEL/Feature.h)
			  for (FeatureMap<>::Iterator fm_it = feature_map.begin(); fm_it != feature_map.end(); ++fm_it)
			  {
				  rt_ok = f.getOptions().getRTRange().encloses(DPosition<1>(fm_it->getRT()));
				  mz_ok = f.getOptions().getMZRange().encloses(DPosition<1>(fm_it->getMZ()));
				  int_ok = f.getOptions().getIntensityRange().encloses(DPosition<1>(fm_it->getIntensity()));
				  charge_ok = ((charge_l <= fm_it->getCharge()) && (fm_it->getCharge() <= charge_u));
				  size_ok = ((size_l <= fm_it->getSubordinates().size()) && (fm_it->getSubordinates().size() <= size_u));
				  q_ok = ((q_l <= fm_it->getOverallQuality()) && (fm_it->getOverallQuality() <= q_u));
				  annotation_ok = checkPeptideIdentification_(*fm_it, remove_annotated_features, remove_unannotated_features, sequences, accessions, keep_best_score_id, remove_clashes);
          if (check_meta) meta_ok = checkMetaOk(*fm_it, meta_info);

				  if (rt_ok && mz_ok && int_ok && charge_ok && size_ok && q_ok && annotation_ok && meta_ok)
				  {
					  map_sm.push_back (*fm_it);
				  }
			  }
			  //delete unassignedPeptideIdentifications
			  if (remove_unassigned_ids)
			  {
				  map_sm.getUnassignedPeptideIdentifications().clear();
			  }
			  //update minimum and maximum position/intensity
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
			  f.getOptions().setRTRange(DRange<1>(rt_l,rt_u));
			  f.getOptions().setMZRange(DRange<1>(mz_l,mz_u));
			  f.getOptions().setIntensityRange(DRange<1>(it_l,it_u));
			  f.load(in,consensus_map);

			  //-------------------------------------------------------------
			  // calculations
			  //-------------------------------------------------------------

			  // copy all properties
			  ConsensusMap consensus_map_filtered = consensus_map;
			  //.. but delete feature information
			  consensus_map_filtered.resize(0);
			
			  bool charge_ok, size_ok, annotation_ok;

        for (ConsensusMap::Iterator cm_it = consensus_map.begin(); cm_it != consensus_map.end(); ++cm_it)
			  {
				  charge_ok = ((charge_l <= cm_it->getCharge()) && (cm_it->getCharge() <= charge_u));
				  size_ok = ((cm_it->size() >= size_l) && (cm_it->size() <= size_u));
				  annotation_ok = checkPeptideIdentification_(*cm_it, remove_annotated_features, remove_unannotated_features, sequences, accessions, keep_best_score_id, remove_clashes);
          if (check_meta) meta_ok = checkMetaOk(*cm_it, meta_info);

				  if (charge_ok && size_ok && annotation_ok && meta_ok)
				  {
					  consensus_map_filtered.push_back(*cm_it);
				  }
			  }
			  //delete unassignedPeptideIdentifications
			  if (remove_unassigned_ids)
			  {
				  consensus_map_filtered.getUnassignedPeptideIdentifications().clear();
			  }
			  //update minimum and maximum position/intensity
			  consensus_map_filtered.updateRanges();
			
			  // sort if desired
			  if (sort)
			  {
				  consensus_map_filtered.sortByPosition();
			  }
			
			  if (out_type == FileTypes::FEATUREXML)
			  {				
				  if (maps.size() == 1) // When extracting a feature map from a consensus map, only one map ID should be specified. Hence 'maps' should contain only one integer.
				  {
					  FeatureMap<> feature_map_filtered;
					  FeatureXMLFile ff;
					
					  for (ConsensusMap::Iterator cm_it=consensus_map_filtered.begin(); cm_it!=consensus_map_filtered.end(); ++cm_it)
					  {
						
						  for(ConsensusFeature::HandleSetType::const_iterator fh_iter = cm_it->getFeatures().begin(); fh_iter != cm_it->getFeatures().end(); ++fh_iter)
						  {
							  if ((int)fh_iter->getMapIndex() == maps[0])
							  {
								  Feature feature;
								  feature.setRT(fh_iter->getRT());
								  feature.setMZ(fh_iter->getMZ());
								  feature.setIntensity(fh_iter->getIntensity());
								  feature.setCharge(fh_iter->getCharge());
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
				  else
          {
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

				  cm_new.setProteinIdentifications(consensus_map_filtered.getProteinIdentifications());

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
				writeLog_("Error: Unknown output file type given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
		}
		else
		{
			writeLog_("Error: Unknown input file type given. Aborting!");
			printUsage_();
			return INCOMPATIBLE_INPUT_DATA;
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
