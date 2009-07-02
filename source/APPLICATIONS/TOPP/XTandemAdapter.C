// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/FORMAT/MascotInfile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

#include <QtCore/QFile>

#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_XTandemAdapter XTandemAdapter
	
	@brief Identifies peptides in MS/MS spectra via XTandem.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_XTandemAdapter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPXTandemAdapter
	: public TOPPBase
{
	public:
		TOPPXTandemAdapter()
			: TOPPBase("XTandemAdapter","Annotates MS/MS spectra using XTandem.")
		{
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			addEmptyLine_();
			addText_("This adapter to X!Tandem provides a small interface with only "
							 "a small number of parameters. Other parameters need to be set "
							 "via the default file. This file is read and the parameters, are "
							 "used to generate the input file for X!Tandem itself. The results "
							 "are converted from the X!Tandem format into the idXML format.");
			addEmptyLine_();
			addText_("Common Identification engine options");
			
			registerInputFile_("in", "<file>", "", "input file ");
      setValidFormats_("in",StringList::create("mzData"));
      registerOutputFile_("out", "<file>", "", "output file ");
      setValidFormats_("out",StringList::create("IdXML"));
			registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 1.5, "precursor mass tolerance", false);
			registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.3, "fragment mass error", false);
			
			registerStringOption_("precursor_error_units", "<unit>", "ppm", "parent monoisotopic mass error units", false);
      registerStringOption_("fragment_error_units", "<unit>", "Da", "fragment monoisotopic mass error units", false);
			registerStringOption_("database", "<file>", "", "FASTA file or related which contains the sequences");
      vector<String> valid_strings;
      valid_strings.push_back("ppm");
      valid_strings.push_back("Da");
      setValidStrings_("precursor_error_units", valid_strings);
      setValidStrings_("fragment_error_units", valid_strings);
			registerIntOption_("min_precursor_charge", "<charge>", 1, "minimum precursor charge", false);
			registerIntOption_("max_precursor_charge", "<charge>", 4, "maximum precursor charge", false);
			
			registerStringOption_("fixed_modifications", "<mods>", "", "fixed modifications, specified using PSI-MOD terms, e.g. MOD:01214,MOD:00048", false);
      registerStringOption_("variable_modifications", "<mods>", "", "variable modifications, specified using PSI-MOD terms, e.g. MOD:01214,MOD:00048", false);
			registerIntOption_("missed_cleavages", "<num>", 1, "Number of possible cleavage sites missed by the enzyme", false);
	
			addEmptyLine_();
			addText_("X!Tandem specific options");
			registerStringOption_("XTandem_path", "<path>", "", "Path to X!Tandem, ending with '/bin'");
			registerInputFile_("default_input_file", "<file>", "default_input.xml", "default parameters input file, if not given default parameters are used", false);			
			registerDoubleOption_("minimum_fragment_mz", "<num>", 150.0, "minimum fragment mz", false);
			registerStringOption_("cleavage_site", "<cleavage site>", "[RK]|{P}", "cleavage site", false);
			registerDoubleOption_("max_valid_expect", "<E-Value>", 0.1, "maximal E-Value of a hit to be reported", false);
			registerFlag_("no_refinement", "Disable the refinement, especially useful for matching only peptides without proteins");

			registerStringOption_("temp_directory", "<dir>", "", "Directory were temporary data can be stored. If not set the directory were startet is used.", false, true);
		}

		ExitCodes main_(int , const char**)
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
			PeakMap exp;
		
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
	
      // write input xml file
			String parameters;
			XTandemInfile infile;


			String unique_name = File::getUniqueName(); // body for the tmp files
			String temp_directory(getStringOption_("temp_directory"));
			if (temp_directory != "")
			{
				temp_directory.ensureLastChar('/');
			}

			String input_filename(temp_directory + unique_name + "_tandem_input_file.xml");
			String tandem_input_filename(temp_directory + unique_name + "_tandem_input_file.mzData");
			String tandem_output_filename(temp_directory + unique_name + "_tandem_output_file.xml");
			String tandem_taxonomy_filename(temp_directory + unique_name + "_tandem_taxonomy_file.xml");
	
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			// only load msLevel 2
			MzDataFile mzdata_infile;
			mzdata_infile.getOptions().addMSLevel(2);
			mzdata_infile.setLogType(log_type_);
			mzdata_infile.load(inputfile_name, exp);

			// we need to replace the native id with a simple numbering schema, to be able to
			// map the IDs back to the spectra (RT, and MZ infomration)
			Size native_id(1);
			for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
			{
				it->setNativeID(native_id++);
			}

			// We store the file in mzData file format, because mgf file somehow produce in most 
			// of the cases ids with charge 2+. We do not use the input file of this TOPP-tools
			// because XTandem sometimes stumbles over misleading substrings in the filename,
			// e.g. mzXML ...
			MzDataFile mzdata_outfile;
			mzdata_outfile.store(tandem_input_filename, exp);

			infile.setInputFilename(tandem_input_filename);
			infile.setOutputFilename(tandem_output_filename);

			
			String fasta_file(getStringOption_("database"));
			
			ofstream tax_out(tandem_taxonomy_filename.c_str());
			tax_out << "<?xml version=\"1.0\"?>" << endl;
			tax_out << "\t<bioml label=\"x! taxon-to-file matching list\">" << endl;
  		tax_out << "\t\t<taxon label=\"OpenMS_dummy_taxonomy\">" << endl;
    	tax_out << "\t\t\t<file format=\"peptide\" URL=\"" << fasta_file << "\" />" << endl;
  		tax_out << "\t</taxon>" << endl;
			tax_out << "</bioml>" << endl;
			tax_out.close();

			infile.setTaxonomyFilename(tandem_taxonomy_filename);

			if (getStringOption_("precursor_error_units") == "Da")
			{
				infile.setPrecursorMassErrorUnit(XTandemInfile::DALTONS);
			}
			else
			{
				infile.setPrecursorMassErrorUnit(XTandemInfile::PPM);
			}
			
			if (getStringOption_("fragment_error_units") == "Da")
			{
				infile.setFragmentMassErrorUnit(XTandemInfile::DALTONS);
			}
			else
			{
				infile.setFragmentMassErrorUnit(XTandemInfile::PPM);
			}

			if (setByUser_("default_input_file"))
			{
				infile.load(getStringOption_("default_input_file"));
				infile.setDefaultParametersFilename(getStringOption_("default_input_file"));
			}
			else
			{
				String default_file = File::find("CHEMISTRY/XTandem_default_input.xml");
				infile.load(default_file);
				infile.setDefaultParametersFilename(default_file);
			}

			infile.setPrecursorMassTolerancePlus(getDoubleOption_("precursor_mass_tolerance"));
			infile.setPrecursorMassToleranceMinus(getDoubleOption_("precursor_mass_tolerance"));
			infile.setFragmentMassTolerance(getDoubleOption_("fragment_mass_tolerance"));
			infile.setMaxPrecursorCharge(getIntOption_("max_precursor_charge"));
			infile.setNumberOfThreads(getIntOption_("threads"));
			infile.setModifications(ModificationDefinitionsSet(getStringOption_("fixed_modifications"), getStringOption_("variable_modifications")));
			infile.setTaxon("OpenMS_dummy_taxonomy");
			infile.setMaxValidEValue(getDoubleOption_("max_valid_expect"));
			infile.setNumberOfMissedCleavages(getIntOption_("missed_cleavages"));
			infile.write(input_filename);
			
			vector<ProteinIdentification> protein_identifications;
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------

			// @todo translate call to windows
			String call = tandem_path + "/./tandem.exe " + input_filename;
			int status = system(call.c_str());

			if (status != 0)
			{
				writeLog_("XTandem problem. Aborting! (Details can be seen in the logfile: \"" + logfile + "\")");
				// clean temporary files
				QFile(input_filename.toQString()).remove();
      	QFile(tandem_input_filename.toQString()).remove();
      	QFile(tandem_taxonomy_filename.toQString()).remove();
				return EXTERNAL_PROGRAM_ERROR;
			}

			vector<ProteinIdentification> protein_ids;
			ProteinIdentification protein_id;
			vector<PeptideIdentification> peptide_ids;

			// read the output of X!Tandem and write it to IdXML
			XTandemXMLFile tandem_output;
			tandem_output.setModificationDefinitionsSet(ModificationDefinitionsSet(getStringOption_("fixed_modifications"), getStringOption_("variable_modifications")));
			// find the file, because XTandem extends the filename with a timestamp we do not know (exactly)
			StringList files;
			File::fileList(temp_directory, unique_name + "_tandem_output_file*.xml", files);
			if (files.size() != 1)
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, tandem_output_filename);
			}
			tandem_output.load(temp_directory + files[0], protein_id, peptide_ids);
			
			// now put the RTs into the peptide_ids from the spectrum ids
			for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
			{
				UInt id = (Int)it->getMetaValue("spectrum_id");
				id -= 1; // native IDs were written 1-based
				if (id < exp.size())
				{
					it->setMetaValue("RT", exp[id].getRT());
					DoubleReal pre_mz = 0.0;
					if (!exp[id].getPrecursors().empty()) pre_mz = exp[id].getPrecursors()[0].getMZ();
					it->setMetaValue("MZ", pre_mz);
					it->removeMetaValue("spectrum_id");
				}
				else
				{
					cerr << "XTandemAdapter: Error: id '" << id << "' not found in peak map!" << endl;
				}
			}

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
	
			// handle the search parameters
      ProteinIdentification::SearchParameters search_parameters;
      search_parameters.db = getStringOption_("database");
      search_parameters.charges = "+" + String(getIntOption_("min_precursor_charge")) + "-+" + String(getIntOption_("max_precursor_charge"));

      ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;

      search_parameters.mass_type = mass_type;

      vector<String> fixed_mods, var_mods;
      getStringOption_("fixed_modifications").split(',', fixed_mods);
      if (fixed_mods.size() == 0)
      {
        if (getStringOption_("fixed_modifications") != "")
        {
          fixed_mods.push_back(getStringOption_("fixed_modifications"));
        }
      }
      getStringOption_("variable_modifications").split(',', var_mods);
      if (var_mods.size() == 0)
      {
        if (getStringOption_("variable_modifications") != "")
        {
          var_mods.push_back(getStringOption_("variable_modifications"));
        }
      }

      search_parameters.fixed_modifications = fixed_mods;
      search_parameters.variable_modifications = var_mods;

      search_parameters.missed_cleavages = getIntOption_("missed_cleavages");
      search_parameters.peak_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");
      search_parameters.precursor_tolerance = getDoubleOption_("precursor_mass_tolerance");


      protein_id.setSearchParameters(search_parameters);
      protein_id.setSearchEngineVersion("");
      protein_id.setSearchEngine("XTandem");

			protein_ids.push_back(protein_id);
			
			IdXMLFile id_output;
			id_output.store(outputfile_name, protein_ids, peptide_ids);

			/// Deletion of temporary files	
			QFile(input_filename.toQString()).remove();
			QFile((temp_directory + files[0]).toQString()).remove();
			QFile(tandem_input_filename.toQString()).remove();
			QFile(tandem_taxonomy_filename.toQString()).remove();
			
			return EXECUTION_OK;	
		}
};


int main( int argc, const char** argv )
{
	TOPPXTandemAdapter tool;

	return tool.main(argc,argv);
}

/// @endcond
