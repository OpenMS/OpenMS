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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page OMSSAAdapter OMSSAAdapter
	
	@brief Identifies peptides in MS/MS spectra via OMSSA (Open Mass Spectrometry Search Algorithm).
	
	todo
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPOMSSAAdapter
	: public TOPPBase
{
	public:
		TOPPOMSSAAdapter()
			: TOPPBase("OMSSAAdapter","annotates MS/MS spectra using OMSSA")
		{
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerStringOption_("out", "<file>", "", "output file in analysisXML format.");
      registerStringOption_("in", "<file>", "", "input file in mzData format.");
			
			//Sequence library
			//-d <String> Blast sequence library to search.  Do not include .p* filename suffixes.
			//-pc <Integer> The number of pseudocounts to add to each precursor mass bin.
			registerStringOption_("d", "<file>", "", "Blast sequence library to search.  Do not include .p* filename suffixes", true);
			registerIntOption_("pc", "<Integer>", 1, "The number of pseudocounts to add to each precursor mass bin", false);
			registerStringOption_("omssa_dir", "<Directory>", "", "The directory of the OMSSA installation", true);
			
			//Input format and filename
			//-f <String> single dta file to search
			//-fx <String> multiple xml-encapsulated dta files to search
			//-fb <String> multiple dta files separated by blank lines to search
			//-fm <String> mgf formatted file
			//-fp <String> pkl formatted file
			//-hs <Integer> the minimum number of m/z values a spectrum must have to be searched
			//-fxml <String> omssa xml search request file (contains search parameters and spectra. overrides command line)
			//-pm <String> search parameter input in xml format (contains search parameters but no spectra. overrides command line except for name of file containing spectra)
			// input options are not all necessary as TOPP tools only accept mzData
			registerIntOption_("hs", "<Integer>", 4, "the minimum number of m/z values a spectrum must have to be searched", false);
			registerStringOption_("pm", "<file>", "", "search parameter input in xml format", false);
			
			//Output results
			//-o <String> filename for text asn.1 formatted search results
			//-ob <String> filename for binary asn.1 formatted search results
			//-ox <String> filename for xml formatted search results
			//-oc <String> filename for comma separated value (excel .csv) formatted search results
			// output options of OMSSA are not necessary 
			
			//The following options output the search parameters and search spectra in the output results. This is necessary for viewing result in the OMSSA browser:
			//-w include spectra and search params in search results

			//To turn off informational messages (but not error messages), use:
			//-ni don't print informational messages
 
			//Mass type and tolerance
			//-to <Real> product ion mass tolerance in Da
			//-te <Real> precursor ion mass tolerance in Da
			//-tez <Integer> scaling of precursor mass tolerance with charge (0 = none, 1= linear)
			registerDoubleOption_("to", "<Real>", 0.8, "product ion mass tolerance in Da", false);
			registerDoubleOption_("te", "<Real>", 2.0, "precursor ion mass tolerance in Da", false);
			registerIntOption_("tez", "<Integer>", 1, "scaling of precursor mass tolerance with charge (0 = none, 1= linear)", false);
			
			//A precursor ion is the ion before fragmentation and the product ions are the ions generated after fragmentation. These values are specified in Daltons +/- the measured value, e.g. a value of 2.0 means +/- 2.0 Daltons of the measured value.
			//The tez value allows you to specify how the mass tolerance scales with the charge of the precursor. For example, you may search a precursor assuming that it has a charge state of 2+ and 3+. If you set tez to 1, then the mass tolerance for the +2 charge state will be 2 times the precursor mass tolerance, and for the 3+ charge state it will be 3 times the precursor mass tolerance. If you set tez to 0, the mass tolerance will always be equal to the precursor mass tolerance, irrespective of charge state.

			//-tom <Integer> product ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact.
			//-tem <Integer> precursor ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact.
			registerIntOption_("tom", "<Integer>", 0, "product ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact", false);
			registerIntOption_("tem", "<Integer>", 0, "precursor ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact", false);
			
			//Monoisotopic searching searches spectral peaks that correspond to peptides consisting entirely of carbon-12. Average mass searching searches on the average natural isotopic mass of peptides. Exact mass searches on the most abundant isotopic peak for a given mass range.

			//-tex <Double> threshold in Da above which the mass of a neutron should be added in an exact mass search.
			registerDoubleOption_("tex", "<Real>", 1446.94, "threshold in Da above which the mass of a neutron should be added in an exact mass search", false);

			//Preprocessing
			//Preprocessing is the process of eliminating noise from a spectrum. Normally, you do not need to adjust these options as OMSSA automatically adjusts its preprocessing for best results.

			//-cl <Real> low intensity cutoff as a fraction of max peak
			//-ch <Real> high intensity cutoff as a fraction of max peak
			//-ci <Real> intensity cutoff increment as a fraction of max peak
			//-w1 <Integer> single charge window in Da
			//-w2 <Integer> double charge window in Da
			//-h1 <Integer> number of peaks allowed in single charge window
			//-h2 <Integer> number of peaks allowed in double charge window
			//-cp <Integer> eliminate charge reduced precursors in spectra (0=no, 1=yes). Typically turned on for ETD spectra.
			// TODO	

			//Charge Handling
			//Determination of precursor charge and product ion charges.  Presently, OMSSA estimates which precursors are 1+.  All other precursors are searched with charge from the minimum to maximum precursor charge specified.
			//-zl <Integer> minimum precursor charge to search when not 1+
			//-zh <Integer> maximum precursor charge to search when not 1+
			//-zt <Integer> minimum precursor charge to start considering multiply charged products 
			//-z1 <Double> the fraction of peaks below the precursor used to determine if the spectrum is charge +1 
			//-zc <Integer> should charge +1 be determined algorithmically (1=yes)
			//-zcc <Integer> how should precursor charges be determined? (1=believe the input file,2=use the specified range)
			//-zoh <Integer> set the maximum product charge to search
			registerIntOption_("zl", "<Integer>", 1, "minimum precursor charge to search when not 1+", false);
			registerIntOption_("zh", "<Integer>", 3, "maximum precursor charge to search when not 1+", false);
			registerIntOption_("zt", "<Integer>", 3, "minimum precursor charge to start considering multiply charged products", false);
			registerDoubleOption_("z1", "<Real>", 0.95, "the fraction of peaks below the precursor used to determine if the spectrum is charge +1", false);
			registerIntOption_("zc", "<Integer>", 1, "should charge +1 be determined algorithmically (1=yes)", false);
			registerIntOption_("zcc", "<Integer>", 2, "how should precursor charges be determined? (1=believe the input file,2=use the specified range)", false);
			registerIntOption_("zoh", "<Integer>", 2, "set the maximum product charge to search", false);

			//Enzyme specification
			//Additional enzymes can be added upon request.
			//-v <Integer> number of missed cleavages allowed
			//-e <Integer> id number of enzyme to use (trypsin is the default)
			//-el print a list of enzymes and their corresponding id number
			//-no <Integer> minimum size of peptides for no-enzyme and semi-tryptic searches
			//-nox <Integer> maximum size of peptides for no-enzyme and semi-tryptic searches
			registerIntOption_("v", "<Integer>", 1, "number of missed cleavages allowed", false);
			registerIntOption_("e", "<Integer>", 0, "id number of enzyme to use (trypsin is the default)", false);
			registerIntOption_("no", "<Integer>", 4, "minimum size of peptides for no-enzyme and semi-tryptic searches", false);
			registerIntOption_("nox", "<Integer>", 40, "maximum size of peptides for no-enzyme and semi-tryptic searches", false);
			
			//Ions to search
			//OMSSA searches two ions series, both of which can be specified.  Normally one of the ion series specified is a forward ion series and the other is a reverse ion series.
			//-il print a list of ions and their corresponding id number
			//-i comma delimited list of id numbers of ions to search 
			//-sp <Integer> number of product ions to search 
			//-sb1 <Integer> should first forward (e.g. b1) product ions be searched (1 = no, 0 = yes)
			//-sct <Integer> should c terminus ions (e.g. y1) be searched (1 = no, 0 = yes)
			registerStringOption_("i", "<Num>,<Num>,<Num>", "1,4", "comma delimited list of id numbers of ions to search", false);
			registerIntOption_("sp", "<Integer>", 100, "number of product ions to search", false);
			registerIntOption_("sb1", "<Integer>", 1, "should first forward (e.g. b1) product ions be searched (1 = no, 0 = yes)", false);
			registerIntOption_("sct", "<Integer>", 0, "should c terminus ions (e.g. y1) be searched (1 = no, 0 = yes)", false);
			
			//Taxonomy
			//By default, OMSSA searches without limiting by taxonomy.  By specifying an NCBI taxonomy id, you can limit your search to a particular organism.  The taxonomy id can by found by searching the NCBI tax browser (enter the scientific name of the organism of interest in the search box and then click the correct search result and then the scientific name in the taxonomy browser to get the numeric taxonomy id).
			//-x comma delimited list of NCBI taxonomy ids to search (0 = all.  This is the default)
			registerStringOption_("x", "<Num>,<Num>,<Num>", "0", "comma delimited list of NCBI taxonomy ids to search (0 = all.  This is the default)", false);
			
			//Search heuristic parameters
			//These are options that can speed up the search.  They can result in decreased sensitivity
			//-hm <Integer> the minimum number of m/z matches a sequence library peptide must have for the hit to the peptide to be recorded
			//-ht <Integer> number of m/z values corresponding to the most intense peaks that must include one match to the theoretical peptide
			registerIntOption_("hm", "<Integer>", 2, "the minimum number of m/z matches a sequence library peptide must have for the hit to the peptide to be recorded", false);
			registerIntOption_("ht", "<Integer>", 6, "number of m/z values corresponding to the most intense peaks that must include one match to the theoretical peptide", false);

			//Results
			//-hl <Integer> maximum number of hits retained for one spectrum
			//-he <Double> the maximum e-value allowed in the hit list
			registerIntOption_("hl", "<Integer>", 30, "maximum number of hits retained for one spectrum", false);
			registerDoubleOption_("he", "<Real>", 1, "the maximum e-value allowed in the hit list", false);

			//Post translational modifications
			//To specify modifications, first type in "omssacl -ml" to see a list of modifications available and their corresponding id number.  Then when running the search, specify the id numbers of the modification you wish to apply, e.g. "omssacl -mf 5 -mv 1,8 ...". Multiple PTMs can be specified by placing commas between the numbers without any spaces.  At the present time, the list of allowed post translational modifications will be expanded over time.
			//-mf  comma delimited list of id numbers for fixed modifications
			//-mv  comma delimited list of id numbers for variable modifications
			//-ml  print a list of modifications and their corresponding id number
			registerStringOption_("mf", "<Num>,<Num>,<Num>", "", "comma delimited list of id numbers for fixed modifications", false);
			registerStringOption_("mv", "<Num>,<Num>,<Num>", "", "comma delimited list of id numbers for variable modifications", false);
			
			//To add your own user defined modifications, edit the usermod0-29 entries in the mods.xml file. If it is common modification, please contact NCBI so that it can be added to the standard list.
			//To reduce the combinatorial expansion that results when specifying multiple variable modifications, you can put an upper bound on the number of mass ladders generated per peptide using the -mm option.  The ladders are generated in the order of the least number of modification to the most number of modifications.
			//-mm <Integer> the maximum number of mass ladders to generate per database peptide
			registerIntOption_("mm", "<Integer>", 128, "the maximum number of mass ladders to generate per database peptide", false);

			//There is an upper bound on the number of combinations of variable mods that can be applied to a peptide from the sequence library. The hard upper bound is 1024, which effectively limits the number of variable modification sites per peptide for an exhaustive search to 10. If you set this number too low, you will miss highly modified peptides. If you set it too high, it will make the e-values less significant by searching for too many possible modifications.
			//OMSSA treats cleavage of the initial methionine in each protein record as a variable modification by default. To turn off this behavior use the command line option
			//-mnm n-term methionine should not be cleaved
			registerFlag_("mnm", "n-term methionine should not be cleaved");
			
			//Iterative searching
			//-is <Double> evalue threshold to include a sequence in the iterative search, 0 = all
			//-ir <Double> evalue threshold to replace a hit, 0 = only if better
			//-ii <Double> evalue threshold to iteratively search a spectrum again, 0 = always
			//registerDoubleOption_("is", "<Real>", "0.0", "evalue threshold to include a sequence in the iterative search, 0 = all", false);
			//registerDoubleOption_("ir", "<Real>", "0.0", "evalue threshold to replace a hit, 0 = only if better", false);
			//registerDoubleOption_("ii", "<Real>", "0.0", "evalue threshold to iteratively search a spectrum again, 0 = always", false);
			
			//-foms <String> read in search result in .oms format (binary asn.1). 
			//-fomx <Double> read in search result in .omx format (xml). 
			//Iterative searching is the ability to re-search search results in hopes of increasing the number of spectra identified. To accomplish this, an iterative search may change search parameters, such as using a no-enzyme search, or restrict the sequence search library to sequences already hit.
			
		}

		ExitCodes main_(int , char**)
		{
			// instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
			String ini_location;
			// path to the log file
			String logfile(getStringOption_("log"));
			String omssa_dir(getStringOption_("omssa_dir"));
			// log filestream (as long as the real logfile is not determined yet)
			ofstream log;
			String inputfile_name;
			String outputfile_name;
			String omssa_outfile_name("omssa_tmp_output.xml");
			PeakMap map;
			
			String parameters;
			
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			
			inputfile_name = getStringOption_("in");			
			writeDebug_(String("Input file: ") + inputfile_name, 1);
			if (inputfile_name == "")
			{
				writeLog_("No input file specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
	
			outputfile_name = getStringOption_("out");
			writeDebug_(String("Output file: ") + outputfile_name, 1);
			if (outputfile_name == "")
			{
				writeLog_("No output file specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}				
		
			parameters += " -d ";
			parameters += getStringOption_("d");
			parameters += " -to ";
			parameters += String(getDoubleOption_("to"));
			parameters += " -te ";
			parameters += String(getDoubleOption_("te"));
			parameters += " -zl ";
			parameters += String(getIntOption_("zl"));
			parameters += " -zh ";
			parameters += String(getIntOption_("zh"));
			parameters += " -f ";
			String omssa_tmp_filename("omssa_tmp.dta");
			parameters += omssa_tmp_filename;
			parameters += " -ox ";
			parameters += omssa_outfile_name;
			parameters += " -ni ";
			parameters += " -he ";
			parameters += String(getDoubleOption_("he"));
			if (getStringOption_("mf") != "")
			{
				parameters += " -mf " + getStringOption_("mf");
			}
			if (getStringOption_("mv") != "")
			{
				parameters += " -mv " + getStringOption_("mv");
			}

			
			//-------------------------------------------------------------
			// testing whether input and output files are accessible
			//-------------------------------------------------------------
			
			inputFileReadable_(inputfile_name);
			outputFileWritable_(outputfile_name);
	
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			MzDataFile mzdata_infile;
			ProteinIdentification protein_identification;
			vector<Identification> identifications;
			mzdata_infile.load(inputfile_name, map);
				
			vector<IdentificationData> peptide_ids;
			vector<ProteinIdentification> protein_identifications;
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
	
			Size i(0);
			for (PeakMap::ConstIterator it = map.begin(); it != map.end(); ++it)
			{
				PeakSpectrum spec(*it);
				spec.getPrecursorPeak().setCharge(1);
				DTAFile().store(omssa_tmp_filename, spec);
				
				String call = omssa_dir + "/omssacl " + parameters;

				writeDebug_(String(++i) + "/" + String(map.size()), 2);
				int status = system(call.c_str());
		
				if (status != 0)
				{
					writeLog_("OMSSA problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
					//return EXTERNAL_PROGRAM_ERROR;
				}

				// read OMSSA output
				OMSSAXMLFile omssa_out_file;
				vector<IdentificationData> tmp_peptide_ids;
				ProteinIdentification tmp_protein_id;
				omssa_out_file.load(omssa_outfile_name, tmp_protein_id, tmp_peptide_ids);

				if (tmp_peptide_ids.size() == 1)
				{
					peptide_ids.push_back(tmp_peptide_ids[0]);
					protein_identifications.push_back(tmp_protein_id);
				}
				else
				{
					peptide_ids.push_back(IdentificationData());
					protein_identifications.push_back(tmp_protein_id);
				}
			}
				
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			AnalysisXMLFile().store(outputfile_name, protein_identifications, peptide_ids);
													 		 												 		 
			/// Deletion of temporary files
			
			return EXECUTION_OK;	
		}
};


int main( int argc, char ** argv )
{
	TOPPOMSSAAdapter tool;

	return tool.main(argc,argv);
}

/// @endcond
