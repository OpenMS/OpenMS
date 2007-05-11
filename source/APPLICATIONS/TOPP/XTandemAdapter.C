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
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
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
	@page XTandemAdapter XTandemAdapter
	
	@brief Identifies peptides in MS/MS spectra via XTandem (Open Mass Spectrometry Search Algorithm).
	
	todo
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPXTandemAdapter
	: public TOPPBase
{
	public:
		TOPPXTandemAdapter()
			: TOPPBase("XTandemAdapter","annotates MS/MS spectra using XTandem")
		{
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerStringOption_("out", "<file>", "", "output file in IdXML format.");
			registerStringOption_("in", "<file>", "", "input file in mzData format.");


			
			//<note type="input" label="list path, default parameters">default_input.xml</note>
			// needed???
			registerStringOption_("default_input_file", "<file>", "", "default parameters input file", false);
			
			registerStringOption_("taxonomy_file", "<file>", "", "taxonomy information file");
		  //<note type="input" label="list path, taxonomy information">taxonomy.xml</note>
			
			registerDoubleOption_("fragment_mass_error", "<error>", 0.4, "fragment monoisotopic mass error", false);
  		//<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
			
			// 100??? was ist das 
  		registerDoubleOption_("parent_mass_error_plus", "<error>", 100, "parent monoisotopic mass error plus", false);
			registerDoubleOption_("parent_mass_error_minus", "<error>", 100, "parent monoisotopic mass error plus", false);
			//<note type="input" label="spectrum, parent monoisotopic mass error plus">100</note>
  		//<note type="input" label="spectrum, parent monoisotopic mass error minus">100</note>
  		
			registerStringOption_("parent_mono_mass_error", "<yes|no>", "yes", "parent monoisotopic mass isotope error", false); 
			//<note type="input" label="spectrum, parent monoisotopic mass isotope error">yes</note>
  		//<note type="input" label="spectrum, fragment monoisotopic mass error units">Daltons</note>
			registerStringOption_("fragment_mono_mass_error_units", "<unit>", "Daltons", "fragment monoisotopic mass error units", false);

  		//<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>
			registerStringOption_("parent_mono_mass_error_units", "<ppm>", "ppm", "parent monoisotopic mass error units", false);


			//<note type="input" label="spectrum, fragment mass type">monoisotopic</note>
			registerStringOption_("fragment_mass_type", "<type>", "monoisotopic", "spectrum, fragment mass type", false);

			//<note type="input" label="spectrum, dynamic range">100.0</note>
			registerDoubleOption_("dynamic_range", "<range>", 100.0, "dynamic range", false);

			registerIntOption_("total_peaks", "<number>", 50, "total peaks", false);

			//<note type="input" label="spectrum, total peaks">50</note>
			

 			//<note type="input" label="spectrum, maximum parent charge">4</note>
			registerIntOption_("maximum_parent_charge", "<charge>", 4, "maximum parent charge", false);
  
	
			//<note type="input" label="spectrum, use noise suppression">yes</note>
			registerFlag_("use_noise_supression", "uses noise supression");


		 	// <note type="input" label="spectrum, minimum parent m+h">500.0</note>
			registerDoubleOption_("minimum_parent_mh", "<num>", 500.0, "minimum parent m+h", false);

			// <note type="input" label="spectrum, minimum fragment mz">150.0</note>
			registerDoubleOption_("minimum_fragment_mz", "<num>", 150.0, "minimum fragment mz", false);

			// <note type="input" label="spectrum, minimum peaks">15</note>
			registerIntOption_("minimum_peaks", "<num>", 15, "minimum number of peaks", false);

  		// <note type="input" label="spectrum, threads">1</note>
			registerIntOption_("threads", "<num>", 1, "number of threads", false);

			registerIntOption_("sequence_batch_size", "<num>", 1000, "sequence batch size", false);
			// <note type="input" label="spectrum, sequence batch size">1000</note>

  		// <note type="input" label="residue, modification mass">57.022@C</note>
			registerStringOption_("modification_mass", "<mod>", "57.022@C", "modification mass", false);

  		// <note type="input" label="residue, potential modification mass"></note>
			registerStringOption_("potential_modification_mass", "<mass>", "", "modification mass", false);

  		// <note type="input" label="residue, potential modification motif"></note>
			registerStringOption_("potential_modification_motif", "<motif>", "", "potential modification motif", false);

  		// <note type="input" label="protein, taxon">yeast</note>
			registerStringOption_("taxon", "<taxon>", "yeast", "taxon", false);

  		// <note type="input" label="protein, cleavage site">[RK]|{P}</note>
			registerStringOption_("cleavage_site", "<cleavage site>", "[RK]|{P}", "cleavage site", false);

  		// <note type="input" label="protein, modified residue mass file"></note>
			registerStringOption_("modified_residue_mass_file", "<file>", "", "modified residue mass file", false);

			registerDoubleOption_("cleavage_C_terminal_mass_change", "<mass diff>", 17.002735, "cleavage C-terminal mass change", false);

  		// <note type="input" label="protein, cleavage C-terminal mass change">+17.002735</note>
	
	
  		// <note type="input" label="protein, cleavage N-terminal mass change">+1.007825</note>
			registerDoubleOption_("cleavage_N_terminal_mass_change", "<mass diff>", 1.007825, "cleavage N-terminal mass change", false);

  		// <note type="input" label="protein, N-terminal residue modification mass">0.0</note>
			registerDoubleOption_("N_terminal_residue_modification_mass", "<mass diff>", 0.0, "protein, N-terminal residue modification mass", false);

	  	// <note type="input" label="protein, C-terminal residue modification mass">0.0</note>
			registerDoubleOption_("C_terminal_residue_modification_mass", "<mass diff>", 0.0, "protein, C-terminal residue modification mass", false);

  		// <note type="input" label="protein, homolog management">no</note>
			registerFlag_("homolog_management", "homolog management");

			registerStringOption_("XTandem_path", "<path>", "", "Path to X!Tandem, ending with '/bin'");

		}

		ExitCodes main_(int , char**)
		{
			// instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
			String ini_location;
			// path to the log file
			String logfile(getStringOption_("log"));
			String tandem_path(getStringOption_("XTandem_path"));
			// log filestream (as long as the real logfile is not determined yet)
			ofstream log;
			String inputfile_name;
			String outputfile_name;
			String tandem_outfile_name("tandem_tmp_output.xml");
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
		
			String tandem_tmp_filename("tandem_tmp.dta");
			
			//-------------------------------------------------------------
			// testing whether input and output files are accessible
			//-------------------------------------------------------------
			
			inputFileReadable_(inputfile_name);
			outputFileWritable_(outputfile_name);
	
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			MzDataFile mzdata_infile;
			mzdata_infile.setLogType(log_type_);
			Identification protein_identification;
			vector<PeptideIdentification> peptide_ids;
			mzdata_infile.load(inputfile_name, map);
				
			vector<Identification> protein_identifications;
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
	
			UInt i(0);
			for (PeakMap::ConstIterator it = map.begin(); it != map.end(); ++it)
			{
				PeakSpectrum spec(*it);
				spec.getPrecursorPeak().setCharge(1);
				DTAFile().store(tandem_tmp_filename, spec);
				
				String call = tandem_path + "/tandem.exe " + parameters;

				writeDebug_(String(++i) + "/" + String(map.size()), 2);
				int status = system(call.c_str());
		
				if (status != 0)
				{
					writeLog_("XTandem problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
					//return EXTERNAL_PROGRAM_ERROR;
				}

				// read XTandem output
				XTandemXMLFile tandem_out_file;
				vector<PeptideIdentification> tmp_peptide_ids;
				Identification tmp_protein_id;
				tandem_out_file.load(tandem_outfile_name, tmp_protein_id, tmp_peptide_ids);

				if (tmp_peptide_ids.size() == 1)
				{
					peptide_ids.push_back(tmp_peptide_ids[0]);
					protein_identifications.push_back(tmp_protein_id);
				}
				else
				{
					peptide_ids.push_back(PeptideIdentification());
					protein_identifications.push_back(tmp_protein_id);
				}
			}
				
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			IdXMLFile().store(outputfile_name, protein_identifications, peptide_ids);
													 		 												 		 
			/// Deletion of temporary files
			
			return EXECUTION_OK;	
		}
};


int main( int argc, char ** argv )
{
	TOPPXTandemAdapter tool;

	return tool.main(argc,argv);
}

/// @endcond
