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
// $Maintainer: Oliver Kohlbacher $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
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
// Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FileInfo FileInfo
	@brief Shows basic information about the data in a OpenMS readable file.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FileInfo \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
		</tr>
	</table>
</CENTER>

	This tool can show basic information about the data in several peak, feature and consensus feature files. It can
	- show information about the data range of a file (m/z, RT, intensity)
	- show a statistical summary for intensities, qualities, feature widths
	- show an overview of the metadata
	- validate several XML formats against their XML schema
	- check for corrupt data in a file (e.g., duplicate spectra)

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FileInfo.cli

	In order to enrich the resulting data of your anaysis pipeline or to quickly compare different outcomes of your pipeline you can invoke the aforementioned information of your input data and (intermediary) results.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{
	//helper struct for identification data
	struct IdData
	{
		String identifier;
		vector<ProteinIdentification> proteins;
		vector<PeptideIdentification> peptides;
	};

	/// A little helper class to gather (and dump) some statistics from a vector<double>.  Uses statistical functions implemented in GSL.
	struct SomeStatistics
	{
		/**@brief Initialize SomeStatistics from data.

		@note: GSL statistics uses double and so we write double not DoubleReal here and where we use this.
		*/
		SomeStatistics& operator()(vector<double>& data)
		{
			// Sanity check: avoid core dump if no data points present.
			if (data.size() > 0)
			{
				sort(data.begin(), data.end());
				mean = gsl_stats_mean(&data.front(), 1, data.size());
				variance = gsl_stats_variance_m(&data.front(), 1, data.size(), mean);
				min = data.front();
				lowerq = gsl_stats_quantile_from_sorted_data (&data.front(), 1, data.size(), 0.25);
				median = gsl_stats_median_from_sorted_data(&data.front(), 1, data.size());
				upperq = gsl_stats_quantile_from_sorted_data(&data.front(), 1, data.size(), 0.75);
				max = data.back();
			}
			return *this;
		}
		double mean, variance, min, lowerq, median, upperq, max;
	};

	/// Write SomeStatistics to a stream.
	static ostream& operator << (ostream& os, const SomeStatistics& rhs)
	{
		return os <<
			"  mean: " << rhs.mean << endl <<
			"  min: " << rhs.min << endl <<
			"  lower quartile: " << rhs.lowerq << endl <<
			"  median: " << rhs.median << endl <<
			"  upper quartile: " << rhs.upperq << endl <<
			"  max: " << rhs.max << endl <<
			"  variance: " << rhs.variance << endl;
	}
}

class TOPPFileInfo
	: public TOPPBase
{
 public:
	TOPPFileInfo()
		: TOPPBase("FileInfo", "Shows basic information about the file, such as data ranges and file type.")
	{

	}

 protected:

	virtual void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file ");
#ifdef USE_ANDIMS
		setValidFormats_("in", StringList::create("mzData,mzXML,mzML,DTA,DTA2D,cdf,mgf,featureXML,consensusXML,idXML,pepXML,fid"));
#else
		setValidFormats_("in", StringList::create("mzData,mzXML,mzML,DTA,DTA2D,mgf,featureXML,consensusXML,idXML,pepXML,fid"));
#endif
		registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content", false);
#ifdef USE_ANDIMS
		setValidStrings_("in_type",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,cdf,mgf,featureXML,consensusXML,idXML,pepXML,fid"));
#else
		setValidStrings_("in_type",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,mgf,featureXML,consensusXML,idXML,pepXML,fid"));
#endif
		registerOutputFile_("out","<file>","","Optional output file. If '-' or left out, the output is written to the command line.", false);
		registerFlag_("m", "Show meta information about the whole experiment");
		registerFlag_("p", "Shows data processing information");
		registerFlag_("s", "Computes a five-number statistics of intensities, qualities, and widths");
		registerFlag_("d", "Show detailed listing of all spectra and chromatograms (peak files only)");
		registerFlag_("c", "Check for corrupt data in the file (peak files only)");
		registerFlag_("v", "Validate the file only (for mzML, mzData, mzXML, featureXML, idXML, consensusXML, pepXML)");
	}

	ExitCodes outputTo(ostream& os)
	{
		//-------------------------------------------------------------
		// Parameter handling
		//-------------------------------------------------------------

		// File names
		String in = getStringOption_("in");

		// File type
		FileHandler fh;
		FileTypes::Type in_type = fh.nameToType(getStringOption_("in_type"));

		if (in_type == FileTypes::UNKNOWN)
		{
			in_type = fh.getType(in);
			writeDebug_(String("Input file type: ") + fh.typeToName(in_type), 2);
		}

		if (in_type == FileTypes::UNKNOWN)
		{
			writeLog_("Error: Could not determine input file type!");
			return PARSE_ERROR;
		}

		os << endl
				 << "-- General information --" << endl
				 << endl
				 << "file name: " << in << endl
				 << "file type: " <<  fh.typeToName(in_type) << endl;

		MSExperiment<Peak1D> exp;
		FeatureMap<> feat;
		ConsensusMap cons;
		IdData id_data;

		//-------------------------------------------------------------
		// Validation
		//-------------------------------------------------------------
		if (getFlag_("v"))
		{
			bool valid = true;
			os << endl << "Validating " << fh.typeToName(in_type) << " file";
			switch (in_type)
			{
			case FileTypes::MZDATA :
				os << " against XML schema version " << MzDataFile().getVersion() << endl;
				valid = MzDataFile().isValid(in,os);
				break;
			case FileTypes::MZML :
				os << " against XML schema version " << MzMLFile().getVersion() << endl;
				valid = MzMLFile().isValid(in,os);
				break;
			case FileTypes::FEATUREXML :
				os << " against XML schema version " << FeatureXMLFile().getVersion() << endl;
				valid = FeatureXMLFile().isValid(in,os);
				break;
			case FileTypes::IDXML :
				os << " against XML schema version " << IdXMLFile().getVersion() << endl;
				valid = IdXMLFile().isValid(in,os);
				break;
			case FileTypes::CONSENSUSXML :
				os << " against XML schema version " << ConsensusXMLFile().getVersion() << endl;
				valid = ConsensusXMLFile().isValid(in,os);
				break;
			case FileTypes::MZXML :
				os << " against XML schema version " << MzXMLFile().getVersion() << endl;
				valid = MzXMLFile().isValid(in,os);
				break;
			case FileTypes::PEPXML:
				os << " against XML schema version " << PepXMLFile().getVersion() << endl;
				valid = PepXMLFile().isValid(in, os);
				break;
			default:
				os << endl << "Aborted: Validation of this file type is not supported!" << endl;
				return EXECUTION_OK;
			};

			if (valid)
			{
				os << "Success: the file is valid!" << endl;
			}
			else
			{
				os << "Failed: errors are listed above!" << endl;
			}

			if (in_type == FileTypes::MZML)
			{
				if (!valid)
				{
					os << endl << "Semantic validation is not perfomed due to previous errors! " << endl;
				}
				else
				{
					os << endl << "Semantically validating " << fh.typeToName(in_type) << " file:" << endl;
					StringList errors, warnings;
					valid = MzMLFile().isSemanticallyValid(in, errors, warnings);
					for (Size i=0; i<warnings.size(); ++i)
					{
						os << "Warning: " << warnings[i] << endl;
					}
					for (Size i=0; i<errors.size(); ++i)
					{
						os << "Error: " << errors[i] << endl;
					}
					if (valid)
					{
						os << "Success: the file is semantically valid!" << endl;
					}
					else
					{
						os << "Failed: errors are listed above!" << endl;
					}
				}
			}
			else if (in_type == FileTypes::MZDATA)
			{
				if (!valid)
				{
					os << endl << "Semantic validation is not perfomed due to previous errors! " << endl;
				}
				else
				{
					os << endl << "Semantically validating " << fh.typeToName(in_type) << " file (EXPERIMENTAL) :" << endl;
					StringList errors, warnings;
					valid = MzDataFile().isSemanticallyValid(in, errors, warnings);
					for (Size i=0; i<warnings.size(); ++i)
					{
						os << "Warning: " << warnings[i] << endl;
					}
					for (Size i=0; i<errors.size(); ++i)
					{
						os << "Error: " << errors[i] << endl;
					}
					if (valid)
					{
						os << "Success: the file is semantically valid!" << endl;
					}
					else
					{
						os << "Failed: errors are listed above!" << endl;
					}
				}
			}

			return EXECUTION_OK;
		}

		//-------------------------------------------------------------
		// Content statistics
		//-------------------------------------------------------------
		Map<String,int> meta_names;
		if (in_type == FileTypes::FEATUREXML) //features
		{
			FeatureXMLFile().load(in,feat);
			feat.updateRanges();

			os << "Number of features: " << feat.size() << endl
				 << endl
				 << "Ranges:" << endl
				 << "  retention time:  " << String::number(feat.getMin()[Peak2D::RT],2) << " : " << String::number(feat.getMax()[Peak2D::RT],2) << endl
				 << "  mass-to-charge:  " << String::number(feat.getMin()[Peak2D::MZ],2) << " : " << String::number(feat.getMax()[Peak2D::MZ],2) << endl
				 << "  intensity:       " << String::number(feat.getMinInt(),2) << " : " << String::number(feat.getMaxInt(),2) << endl
				 << endl;

			// Charge distribution and TIC
			Map<UInt,UInt> charges;
			DoubleReal tic = 0.0;
			for (Size i = 0; i < feat.size(); ++i)
			{
				charges[feat[i].getCharge()]++;
				tic += feat[i].getIntensity();
			}

			os << "Total ion current in features: " << tic << endl;
			os << "Charge distribution" << endl;
			for (Map<UInt, UInt>::const_iterator it = charges.begin(); it != charges.end(); ++it)
			{
				os << "charge " << it->first << ": " << it->second << endl;
			}
		}
		else if (in_type == FileTypes::CONSENSUSXML) //consensus features
		{
			ConsensusXMLFile().load(in,cons);
			cons.updateRanges();

			map<Size,UInt> num_consfeat_of_size;
			for (ConsensusMap::const_iterator cmit = cons.begin(); cmit != cons.end(); ++cmit )
			{
				++num_consfeat_of_size[cmit->size()];
			}

			os << endl << "Number of consensus features:" << endl;
			for (map<Size,UInt>::reverse_iterator i = num_consfeat_of_size.rbegin(); i != num_consfeat_of_size.rend(); ++i )
			{
				os << "  of size " << setw(2) << i->first << ": " << setw(6) << i->second << endl;
			}
			os << "  total:      " << setw(6) << cons.size() << endl << endl;

			os << "Ranges:" << endl
				 << "  retention time:  " << String::number(cons.getMin()[Peak2D::RT],2) << " : " << String::number(cons.getMax()[Peak2D::RT],2) << endl
				 << "  mass-to-charge:  " << String::number(cons.getMin()[Peak2D::MZ],2) << " : " << String::number(cons.getMax()[Peak2D::MZ],2) << endl
				 << "  intensity:       " << String::number(cons.getMinInt(),2) << " : " << String::number(cons.getMaxInt(),2) << endl;

			// file descriptions
			const ConsensusMap::FileDescriptions& descs = cons.getFileDescriptions();
			if (descs.size() != 0)
			{
				os << endl <<
					"File descriptions:" << endl;
				for (ConsensusMap::FileDescriptions::const_iterator it=descs.begin(); it != descs.end(); ++it)
				{
					os << " - " << it->second.filename << endl
							 << "   identifier: " << it->first << endl
							 << "   label     : " << it->second.label << endl
							 << "   size      : " << it->second.size << endl;
				}
			}
		}
		else if (in_type==FileTypes::IDXML) //identifications
		{
			UInt spectrum_count = 0;
			Size peptide_hit_count = 0;
			UInt runs_count = 0;
			Size protein_hit_count = 0;
			vector<String> peptides;
			vector<String> proteins;

			// reading input
			IdXMLFile().load(in, id_data.proteins, id_data.peptides, id_data.identifier);

			// calculations
			for (Size i = 0; i < id_data.peptides.size(); ++i)
			{
				if (!id_data.peptides[i].empty())
				{
					++spectrum_count;
					peptide_hit_count += id_data.peptides[i].getHits().size();
					const vector<PeptideHit>& temp_hits = id_data.peptides[i].getHits();
					for (Size j = 0; j < temp_hits.size(); ++j)
					{
						if (find(peptides.begin(), peptides.end(), temp_hits[j].getSequence().toString()) == peptides.end())
						{
							peptides.push_back(temp_hits[j].getSequence().toString());
						}
					}
				}
			}
			for (Size i = 0; i < id_data.proteins.size(); ++i)
			{
				++runs_count;
				protein_hit_count += id_data.proteins[i].getHits().size();
				const vector<ProteinHit>& temp_hits = id_data.proteins[i].getHits();
				for (Size j = 0; j < temp_hits.size(); ++j)
				{
					if (find(proteins.begin(), proteins.end(), temp_hits[j].getAccession()) == proteins.end())
					{
						proteins.push_back(temp_hits[j].getAccession());
					}
				}
			}

			os << "Number of runs: " << runs_count << endl;
			os << "Number of protein hits: " << protein_hit_count << endl;
			os << "Number of unique protein hits: " << proteins.size() << endl;
			os << endl;
			os << "Number of spectra: " << spectrum_count << endl;
			os << "Number of peptide hits: " << peptide_hit_count << endl;
			os << "Number of unique peptide hits: " << peptides.size() << endl;
		}

		else if (in_type == FileTypes::PEPXML)
		{
			os << "\nFor pepXML files, only validation against the XML schema is implemented at this point." << endl;
		}

		else //peaks
		{
			if (!fh.loadExperiment(in,exp,in_type,log_type_))
			{
				writeLog_("Unsupported or corrupt input file. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			//check if the meta data indicates that this is peak data
			UInt meta_type = SpectrumSettings::UNKNOWN;
			if (exp.size() > 0)
			{
				for (Size i = 0; i < exp[0].getDataProcessing().size(); ++i)
				{
					if (exp[0].getDataProcessing()[i].getProcessingActions().count(DataProcessing::PEAK_PICKING)==1)
					{
						meta_type = SpectrumSettings::PEAKS;
					}
				}
			}
			//determine type (search for the first scan with at least 5 peaks)
			UInt type = SpectrumSettings::UNKNOWN;
			UInt i = 0;
			while(i < exp.size() && exp[i].size() < 5)
			{
				++i;
			}
			if (i != exp.size())
			{
				type = PeakTypeEstimator().estimateType(exp[i].begin(),exp[i].end());
			}
			os << endl
					 << "peak type (metadata) : " << SpectrumSettings::NamesOfSpectrumType[meta_type] << endl
					 << "peak type (estimated): " << SpectrumSettings::NamesOfSpectrumType[type] << endl;
			//if raw data, determine the spacing
			if (type==SpectrumSettings::RAWDATA)
			{
				vector<Real> spacing;
				for (Size j = 1; j < exp[i].size(); ++j)
				{
					spacing.push_back(exp[i][j].getMZ() - exp[i][j-1].getMZ());
				}
				sort(spacing.begin(),spacing.end());
				os << "estimated raw data spacing: " << spacing[spacing.size()/2] << " (min: " << spacing[0] << " max: " << spacing.back() << ")" << endl;
			}
			os << endl;

			//basic info
			exp.updateRanges();
			vector<UInt> levels = exp.getMSLevels();

			os << "Number of spectra: "	<< exp.size() << endl;
			os << "Number of peaks:   " << exp.getSize() << endl
				 << endl
				 << "Ranges:" << endl
				 << "  retention time:  " << String::number(exp.getMinRT(),2) << " : " << String::number(exp.getMaxRT(),2) << endl
				 << "  mass-to-charge:  " << String::number(exp.getMinMZ(),2) << " : " << String::number(exp.getMaxMZ(),2) << endl
				 << "  intensity:       " << String::number(exp.getMinInt(),2) << " : " << String::number(exp.getMaxInt(),2) << endl
				 << endl;

			os << "MS levels: ";
			if (levels.size() != 0)
			{
				os << *(levels.begin());
				for (vector<UInt>::iterator it = ++levels.begin(); it != levels.end(); ++it)
				{
					os << ", " << *it;
				}
			}
			os << endl;

			//count how many spectra per MS level there are
			vector<UInt> counts(5);
			for (MSExperiment<Peak1D>::iterator it = exp.begin(); it != exp.end(); ++it)
			{
				counts[it->getMSLevel()]++;
			}
			//output
			for (Size i = 0; i != 5; ++i)
			{
				if (counts[i]!=0)
				{
					os << "Spectra of MS Level " << i << ": " << counts[i] << endl;
				}
			}
			os << endl;

			// show meta data array names
			for (MSExperiment<Peak1D>::iterator it = exp.begin(); it != exp.end(); ++it)
			{
				for (i = 0; i < it->getFloatDataArrays().size(); ++i)
				{
					String name = it->getFloatDataArrays()[i].getName();
					if (meta_names.has(name))
					{
						meta_names[name]++;
					}
					else
					{
						meta_names[name] = 1;
					}
				}
				for (i = 0; i < it->getIntegerDataArrays().size(); ++i)
				{
					String name = it->getIntegerDataArrays()[i].getName();
					if (meta_names.has(name))
					{
						meta_names[name]++;
					}
					else
					{
						meta_names[name] = 1;
					}
				}
				for (i = 0; i < it->getStringDataArrays().size(); ++i)
				{
					String name = it->getStringDataArrays()[i].getName();
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
					os << "Meta data array: " << it->first << " (for " << it->second << " spectra)" << endl;
				}
				os << endl;
			}

			// some chromatogram information
			if (exp.getChromatograms().size() != 0)
			{
				os << "Number of chromatograms: "	<< exp.getChromatograms().size() << endl;

				Size num_chrom_peaks(0);
				Map<ChromatogramSettings::ChromatogramType, Size> chrom_types;
				for (vector<MSChromatogram<> >::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
				{
					num_chrom_peaks += it->size();
					if (chrom_types.has(it->getChromatogramType()))
					{
						chrom_types[it->getChromatogramType()]++;
					}
					else
					{
						chrom_types[it->getChromatogramType()] = 1;
					}
				}
				os << "Number of chromatographic peaks: " << num_chrom_peaks << endl << endl;

				os << "#Chromatograms of types: " << endl;
				for (Map<ChromatogramSettings::ChromatogramType, Size>::const_iterator it = chrom_types.begin(); it != chrom_types.end(); ++it)
				{
					switch (it->first)
					{
						case ChromatogramSettings::MASS_CHROMATOGRAM:                         os << "   Mass chromatogram:                         " << it->second << endl; break;
						case ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM:            os << "   Total ion current chromatogram:            " << it->second << endl; break;
						case ChromatogramSettings::SELECTED_ION_CURRENT_CHROMATOGRAM:         os << "   Selected ion current chromatogram:         " << it->second << endl; break;
						case ChromatogramSettings::BASEPEAK_CHROMATOGRAM:                     os << "   Basepeak chromaogram:                      " << it->second << endl; break;
						case ChromatogramSettings::SELECTED_ION_MONITORING_CHROMATOGRAM:      os << "   Selected ion monitoring chromatogram:      " << it->second << endl; break;
						case ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM: os << "   Selected reaction monitoring chromatogram: " << it->second << endl; break;
						case ChromatogramSettings::ELECTROMAGNETIC_RADIATION_CHROMATOGRAM:    os << "   Electromagnetic radiation chromatogram:    " << it->second << endl; break;
						case ChromatogramSettings::ABSORPTION_CHROMATOGRAM:                   os << "   Absorption chromatogram:                   " << it->second << endl; break;
						case ChromatogramSettings::EMISSION_CHROMATOGRAM:                     os << "   Emission chromatogram:                     " << it->second << endl; break;
						default: 								                                              os << "   Unknown chromatogram:                      " << it->second << endl;
					}
				}
				if (getFlag_("d") && chrom_types.has(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM))
				{
					os << endl << " -- Detailed chromatogram listing -- " << endl;
					os << "\n#Selected Reaction Monitoring Transitions:" << endl;
					os << "#Q1 Q3 RT-begin RT-end name comment" << endl;
					for (vector<MSChromatogram<> >::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
					{
						if (it->getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
						{
							os << it->getPrecursor().getMZ() << " " << it->getProduct().getMZ() << " " << it->front().getRT() << " " << it->back().getRT() << " " << it->getName() << " " << it->getComment() << endl;
						}
					}
				}
			}

			// Detailed listing of scans
			if (getFlag_("d") && exp.size() > 0)
			{
				os << endl
						 << "-- Detailed spectrum listing --" << endl
						 << endl;
				UInt count=0;
				for (MSExperiment<Peak1D>::iterator it = exp.begin(); it!=exp.end(); ++it)
				{
					++count;
					os << "spectrum " << count << " - mslevel:" << it->getMSLevel() << " scanMode:" << InstrumentSettings::NamesOfScanMode[it->getInstrumentSettings().getScanMode()] << " peaks:" << it->size() << " RT:" << it->getRT() << " m/z:";
					if (it->size()!=0)
					{
						os << it->begin()->getMZ() << "-" << (it->end()-1)->getMZ();
					}
					os << endl;
				}
			}

			//Check for corrupt data
			if (getFlag_("c"))
			{
				os << endl
						 << "-- Checking for corrupt data --" << endl
						 << endl;
				// RTs sorted?
				if (!exp.isSorted(false))
				{
					os << "Error: Spectrum retention times are not sorted in ascending order" << endl;
				}
				vector<DoubleReal> ms1_rts;
				ms1_rts.reserve(exp.size());
				for (Size s = 0; s < exp.size(); ++s)
				{
					// ms level = 0
					if (exp[s].getMSLevel() == 0)
					{
						os << "Error: MS-level 0 in spectrum (RT: " << exp[s].getRT() << ")" << endl;
					}
					//scan size = 0
					if (exp[s].size() == 0)
					{
						os << "Warning: No peaks in spectrum (RT: " << exp[s].getRT() << ")" << endl;
					}
					//duplicate meta data array names
					Map<String,int> names;
					for (Size m = 0; m < exp[s].getFloatDataArrays().size(); ++m)
					{
						String name = exp[s].getFloatDataArrays()[m].getName();
						if (names.has(name))
						{
							os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")" << endl;
						}
						else
						{
							names[name] = 0;
						}
					}
					for (Size m = 0; m < exp[s].getIntegerDataArrays().size(); ++m)
					{
						String name = exp[s].getIntegerDataArrays()[m].getName();
						if (names.has(name))
						{
							os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")" << endl;
						}
						else
						{
							names[name] = 0;
						}
					}
					for (Size m = 0; m < exp[s].getStringDataArrays().size(); ++m)
					{
						String name = exp[s].getStringDataArrays()[m].getName();
						if (names.has(name))
						{
							os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")" << endl;
						}
						else
						{
							names[name] = 0;
						}
					}
					//duplicate scans (part 1)
					if (exp[s].getMSLevel() == 1) 
					{
						ms1_rts.push_back(exp[s].getRT());
					}
				}
				//duplicate scans (part 2)
				sort(ms1_rts.begin(), ms1_rts.end());
				for (Size i = 1; i < ms1_rts.size(); ++i)
				{
					if (ms1_rts[i-1]==ms1_rts[i]) os << "Error: Duplicate spectrum retention time: " << ms1_rts[i] << endl;
				}
				//check peaks
				for (Size s = 0; s < exp.size(); ++s)
				{
					//peaks sorted?
					if (!exp[s].isSorted())
					{
						os << "Error: Peak m/z positions are not sorted in ascending order in spectrum (RT: " << exp[s].getRT() << ")" << endl;
					}
					vector<DoubleReal> mzs;
					mzs.reserve(exp[s].size());
					for (Size p = 0; p < exp[s].size(); ++p)
					{
						//negative intensity
						if (exp[s][p].getIntensity() < 0.0)
						{
							os << "Warning: Negative peak intensity peak (RT: " << exp[s].getRT() << " MZ: " << exp[s][p].getMZ() << " intensity: " << exp[s][p].getIntensity() << ")" << endl;
						}
						//duplicate m/z (part 1)
						mzs.push_back(exp[s][p].getMZ());
					}
					//duplicate m/z (part 2)
					sort(mzs.begin(), mzs.end());
					for (Size i = 1; i < mzs.size(); ++i)
					{
						if (mzs[i-1]==mzs[i]) os << "Error: Duplicate peak m/z " << mzs[i] << " in spectrum (RT: " << exp[s].getRT() << ")" << endl;
					}
				}
			}
		}

		//-------------------------------------------------------------
		// meta information
		//-------------------------------------------------------------
		if (getFlag_("m"))
		{
			//basic info
			os << endl
					 << "-- Meta information --" << endl
					 << endl;

			if (in_type == FileTypes::FEATUREXML) //features
			{
				os << "Document id       : " << feat.getIdentifier() << endl << endl;
			}
			else if (in_type == FileTypes::CONSENSUSXML) //consensus features
			{
				os << "Document id       : " << cons.getIdentifier() << endl << endl;
			}
			else if (in_type == FileTypes::IDXML) //identifications
			{
				os << "Document id       : " << id_data.identifier << endl << endl;
			}
			else if (in_type == FileTypes::PEPXML)
			{
				// TODO
			}
			else //peaks
			{

				os << "Document id       : " << exp.getIdentifier() << endl
						 << "Date              : " << exp.getDateTime().get() << endl;

				//basic info
				os << endl
						 << "Sample" << endl
						 << "  Name             : " << exp.getSample().getName() << endl
						 << "  Organism         : " << exp.getSample().getOrganism()  << endl
						 << "  Comment          : " << exp.getSample().getComment()  << endl;

				//instrument info
				os << endl
						 << "Instrument" << endl
						 << "  Name             : " << exp.getInstrument().getName() << endl
						 << "  Model            : " << exp.getInstrument().getModel()  << endl
						 << "  Vendor           : " << exp.getInstrument().getVendor()  << endl
						 << "  Ion source(s)    : ";
				for (Size i = 0; i< exp.getInstrument().getIonSources().size(); ++i)
				{
					os << IonSource::NamesOfIonizationMethod[exp.getInstrument().getIonSources()[i].getIonizationMethod()];
					if (i != exp.getInstrument().getIonSources().size() - 1) 
					{
						os << ", ";
					}
				}
				os << endl << "  Mass Analyzer(s) : ";
				for (Size i=0; i< exp.getInstrument().getMassAnalyzers().size(); ++i)
				{
					os << MassAnalyzer::NamesOfAnalyzerType[exp.getInstrument().getMassAnalyzers()[i].getType()];
					if (i != exp.getInstrument().getMassAnalyzers().size()-1)
					{
						os << ", ";
					}
				}
				os << endl << "  Detector(s)      : ";
				for (Size i = 0; i < exp.getInstrument().getIonDetectors().size(); ++i)
				{
					os << IonDetector::NamesOfType[exp.getInstrument().getIonDetectors()[i].getType()];
					if (i!=exp.getInstrument().getIonDetectors().size()-1) os << ", ";
				}
				os << endl << endl;

				//contact persons
				for (Size i = 0; i < exp.getContacts().size(); ++i)
				{
					os << "Contact Person" << endl
							 << "  First Name       : " << exp.getContacts()[i].getFirstName() << endl
							 << "  Last Name        : " << exp.getContacts()[i].getLastName() << endl
							 << "  Email            : " << exp.getContacts()[i].getEmail() << endl
							 << endl;
				}
			}
		}


		//-------------------------------------------------------------
		// data processing
		//-------------------------------------------------------------
		if (getFlag_("p"))
		{
			//basic info
			os << endl
				 << "-- Data processing information --" << endl
				 << endl;

			//get data processing info
			vector<DataProcessing> dp;
			if (in_type == FileTypes::FEATUREXML) //features
			{
				dp = feat.getDataProcessing();
			}
			else if (in_type == FileTypes::CONSENSUSXML) //consensus features
			{
				dp = cons.getDataProcessing();
			}
			else if (in_type == FileTypes::IDXML) //identifications
			{
			}
			else if (in_type == FileTypes::PEPXML)
			{
				// TODO
			}
			else //peaks
			{
				if (exp.size() != 0)
				{
					os << "Note: The data is taken from the first spectrum!" << endl << endl;
					dp = exp[0].getDataProcessing();
				}
			}

			//print data
			if (dp.size() == 0)
			{
					os << "No information about data processing available!" << endl << endl;
			}
			else
			{
				for (Size i = 0; i < dp.size(); ++i)
				{
					os << "Processing " << (i+1) << ":" << endl;
					os << "  Software name    : " << dp[i].getSoftware().getName() << endl;
					os << "  Software version : " << dp[i].getSoftware().getVersion() << endl;
					os << "  Completion time  : " << dp[i].getCompletionTime().get() << endl;
					os << "  Actions          :";
					for (set<DataProcessing::ProcessingAction>::const_iterator it = dp[i].getProcessingActions().begin(); 
							 it != dp[i].getProcessingActions().end(); ++it)
					{
						if (it != dp[i].getProcessingActions().begin()) os << ",";
						os << " " << DataProcessing::NamesOfProcessingAction[*it];
					}
					os << endl << endl;
				}
			}
		}

		//-------------------------------------------------------------
		// statistics
		//-------------------------------------------------------------
		if (getFlag_("s"))
		{
			os << endl
					 << "-- Statistics --" << endl
					 << endl;
			OpenMS::SomeStatistics some_statistics;

			if (in_type == FileTypes::FEATUREXML) //features
			{
				Size size = feat.size();

				vector<double> intensities(size);
				vector<double> overall_qualities(size);
				vector<double> mz_qualities(size);
				vector<double> rt_qualities(size);
				vector<double> peak_widths(size);

				Size idx = 0;
				for (FeatureMap<>::const_iterator fm_iter = feat.begin();
						 fm_iter != feat.end(); ++fm_iter, ++idx)
				{
					intensities[idx] = fm_iter->getIntensity();
					overall_qualities[idx] = fm_iter->getOverallQuality();
					rt_qualities[idx] = fm_iter->getQuality(Feature::RT);
					mz_qualities[idx] = fm_iter->getQuality(Feature::MZ);
					peak_widths[idx] = fm_iter->getWidth();
				}

				os.precision(writtenDigits<>(Feature::IntensityType() ));
				os << "Intensities:" << endl << some_statistics(intensities) << endl;

				os.precision(writtenDigits<>(Feature::QualityType()));
				os << "Feature FWHM in RT dimension:" << endl << some_statistics(peak_widths) << endl;

				os.precision(writtenDigits<>(Feature::QualityType() ));
				os << "Overall qualities:" << endl << some_statistics(overall_qualities) << endl;

				os.precision(writtenDigits<>(Feature::QualityType()));
				os << "Qualities in retention time dimension:" << endl << some_statistics(rt_qualities) << endl;

				os.precision(writtenDigits<>(Feature::QualityType()));
				os << "Qualities in mass-to-charge dimension:" << endl << some_statistics(mz_qualities) << endl;

			}
			else if (in_type == FileTypes::CONSENSUSXML) //consensus features
			{
				Size size = cons.size();

				vector<double> intensities;
				intensities.reserve(size);
				vector<double> qualities(size);
				qualities.reserve(size);
				vector<double> widths(size);
				widths.reserve(size);

				vector<double> rt_delta_by_elems;
				vector<double> rt_aad_by_elems;
				vector<double> rt_aad_by_cfs;
				rt_aad_by_cfs.reserve(size);

				vector<double> mz_delta_by_elems;
				vector<double> mz_aad_by_elems;
				vector<double> mz_aad_by_cfs;
				mz_aad_by_cfs.reserve(size);

				vector<double> it_delta_by_elems;
				vector<double> it_aad_by_elems;
				vector<double> it_aad_by_cfs;
				it_aad_by_cfs.reserve(size);

				for (ConsensusMap::const_iterator cm_iter = cons.begin();
						 cm_iter != cons.end(); ++cm_iter)
				{
					double rt_aad = 0;
					double mz_aad = 0;
					double it_aad = 0;
					intensities.push_back(cm_iter->getIntensity());
					qualities.push_back(cm_iter->getQuality());
					widths.push_back(cm_iter->getWidth());
					for (ConsensusFeature::HandleSetType::const_iterator hs_iter = cm_iter->begin();
								hs_iter != cm_iter->end(); ++hs_iter)
					{
						double rt_diff = hs_iter->getRT() - cm_iter->getRT();
						rt_delta_by_elems.push_back(rt_diff);
						if (rt_diff < 0)
						{	
							rt_diff = -rt_diff;
						}
						rt_aad_by_elems.push_back(rt_diff);
						rt_aad += rt_diff;
						double mz_diff = hs_iter->getMZ() - cm_iter->getMZ();
						mz_delta_by_elems.push_back(mz_diff);
						if (mz_diff < 0)
						{
							mz_diff = -mz_diff;
						}
						mz_aad_by_elems.push_back(mz_diff);
						mz_aad += mz_diff;
						double it_ratio = hs_iter->getIntensity() / ( cm_iter->getIntensity() ? cm_iter->getIntensity() : 1. );
						it_delta_by_elems.push_back(it_ratio);
						if (it_ratio < 1.)
						{
							it_ratio = 1./it_ratio;
						}
						it_aad_by_elems.push_back(it_ratio);
						it_aad += it_ratio;
					}
					if (!cm_iter->empty())
					{
						rt_aad /= cm_iter->size();
						mz_aad /= cm_iter->size();
						it_aad /= cm_iter->size();
					} // otherwise rt_aad etc. are 0 anyway
					rt_aad_by_cfs.push_back(rt_aad);
					mz_aad_by_cfs.push_back(mz_aad);
					it_aad_by_cfs.push_back(it_aad);
				}

				os.precision(writtenDigits(ConsensusFeature::IntensityType()));
				os << "Intensities of consensus features:" << endl << some_statistics(intensities) << endl;

				os.precision(writtenDigits(ConsensusFeature::QualityType()));
				os << "Qualities of consensus features:" << endl << some_statistics(qualities) << endl;

				os.precision(writtenDigits(ConsensusFeature::CoordinateType()));
				os << "Retention time differences ( element-center, weight 1 per element):" << endl << some_statistics(rt_delta_by_elems) << endl;
				os << "Absolute retention time differences ( |element-center|, weight 1 per element):" << endl << some_statistics(rt_aad_by_elems) << endl;
				os << "Average absolute differences of retention time within consensus features ( |element-center|, weight 1 per consensus features):" << endl << some_statistics(rt_aad_by_cfs) << endl;

				os.precision(writtenDigits(ConsensusFeature::CoordinateType()));
				os << "Mass-to-charge differences ( element-center, weight 1 per element):" << endl << some_statistics(mz_delta_by_elems) << endl;
				os << "Absolute differences of mass-to-charge ( |element-center|, weight 1 per element):" << endl << some_statistics(mz_aad_by_elems) << endl;
				os << "Average absolute differences of mass-to-charge within consensus features ( |element-center|, weight 1 per consensus features):" << endl << some_statistics(mz_aad_by_cfs) << endl;

				os.precision(writtenDigits(ConsensusFeature::IntensityType()));
				os << "Intensity ratios ( element/center, weight 1 per element):" << endl << some_statistics(it_delta_by_elems) << endl;
				os << "Relative intensity error ( max{(element/center),(center/element)}, weight 1 per element):" << endl << some_statistics(it_aad_by_elems) << endl;
				os << "Average relative intensity error within consensus features ( max{(element/center),(center/element)}, weight 1 per consensus features):" << endl << some_statistics(it_aad_by_cfs) << endl;

			}
			else if (in_type == FileTypes::IDXML) //identifications
			{
				//TODO
			}
			else if (in_type == FileTypes::PEPXML)
			{
				// TODO
			}
			else //peaks
			{
				//copy intensities of  MS-level 1 peaks
				exp.updateRanges(1);
				Size size = exp.getSize();
				vector<double> intensities;
				intensities.reserve(size);
				for (MSExperiment<Peak1D>::const_iterator spec = exp.begin(); spec != exp.end(); ++spec)
				{
					if (spec->getMSLevel() != 1)
					{
						continue;
					}
					for (MSExperiment<Peak1D>::SpectrumType::const_iterator it = spec->begin(); it!=spec->end(); ++it)
					{
						intensities.push_back(it->getIntensity());
					}
				}

				sort(intensities.begin(),intensities.end());
				os.precision(writtenDigits(Peak1D::IntensityType()));
				os << "Intensities:" << endl << some_statistics(intensities) << endl;

				//Statistics for meta information
				for (Map<String,int>::ConstIterator it=meta_names.begin();it!=meta_names.end();++it)
				{
					String name = it->first;
					os << "Meta data: " << name << endl;
					vector<Real> m_values;
					DoubleReal sum = 0.0;
					for (MSExperiment<Peak1D>::const_iterator spec = exp.begin(); spec != exp.end(); ++spec)
					{
						for (Size meta = 0; meta < spec->getFloatDataArrays().size(); ++meta)
						{
							if (spec->getFloatDataArrays()[meta].getName()!=name) continue;
							for (Size peak = 0; peak < spec->getFloatDataArrays()[meta].size(); ++peak)
							{
								m_values.push_back(spec->getFloatDataArrays()[meta][peak]);
								sum += spec->getFloatDataArrays()[meta][peak];
							}
						}
						for (Size meta = 0; meta < spec->getIntegerDataArrays().size(); ++meta)
						{
							if (spec->getIntegerDataArrays()[meta].getName()!=name)
							{
								continue;
							}
							for (Size peak = 0; peak < spec->getIntegerDataArrays()[meta].size(); ++peak)
							{
								m_values.push_back(spec->getIntegerDataArrays()[meta][peak]);
								sum += spec->getIntegerDataArrays()[meta][peak];
							}
						}
					}
					sort(m_values.begin(),m_values.end());
					os << "  count: " << m_values.size() << endl
							 << "  min: " << QString::number(m_values.front(),'f',2).toStdString() << endl
							 << "  max: " << QString::number(m_values.back(),'f',2).toStdString() << endl
							 << "  mean: " << QString::number(sum/m_values.size(),'f',2).toStdString() << endl
							 << "  median: " << QString::number(m_values[m_values.size()/2],'f',2).toStdString() << endl
							 << endl;
				}
			}
		}

		os << endl << endl;

		return EXECUTION_OK;
	}

	ExitCodes main_(int, const char**)
	{
		String out = getStringOption_("out");

		//output to command line
		if (out == "")
		{
			return outputTo(cout);
		}
		//output to file
		else
		{
			ofstream os(out.c_str());
			return outputTo(os);
		}
	}
};

int main(int argc, const char** argv)
{
	TOPPFileInfo tool;
	return tool.main(argc, argv);
}

/// @endcond
