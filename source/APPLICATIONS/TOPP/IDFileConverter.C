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
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_IDFileConverter IDFileConverter

 @brief Converts identification engine file formats.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_IDFileConverter.cli

 @todo Write tests (Clemens, Chris, Hendrik)
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDFileConverter : public TOPPBase
{
public:
  TOPPIDFileConverter() :
    TOPPBase("IDFileConverter", "Converts identification engine file formats.", true)
  {
  }

protected:
  void
  registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<path>", "", "Input file/directory containing the output of the search engine.\n"
      "Sequest: Directory containing the .out files\n"
      "PepXML: Single PepXML file.\n"
      "idXML: Single idXML file.\n", true);
    registerOutputFile_("out", "<file>", "", "Output file", true);
    setValidFormats_("out", StringList::create("idXML,PepXML"));
    registerStringOption_("out_type", "<type>", "", "output file type -- default: determined from file extension or content\n", false);
    setValidStrings_("out_type", StringList::create("idXML,pepXML"));

    addEmptyLine_();
    addText_("Sequest options:");
    registerStringOption_("mz_file", "<file>", "", "Retention times will be looked up in this file, if supplied.\n"
      "Note: Sequest .out files do not contain retention times, but only scan numbers.", false);
    // Please contact the maintainers if you know more about Sequest .out files and might help to resolve this issue
    registerFlag_("ignore_proteins_per_peptide", "Workaround to deal with .out files that contain e.g. \"+1\" in references column,\n"
      "but do not list extra references in subsequent lines (try -debug 3 or 4)", true);

    addEmptyLine_();
    addText_("PepXML options:");
    registerStringOption_("mz_file", "<file>", "", "Retention times will be looked up in this file, if supplied.\n"
    																							 "Note: PepXML files do not contain retention times, but only scan numbers.", false);
		registerStringOption_("mz_name", "<file>", "", "Experiment filename/path to match in the PepXML file ('base_name' attribute);\nonly necessary if different from 'mz_file'.", false);
  }

  ExitCodes
  main_(int, const char**)
  {
    //-------------------------------------------------------------
    // general variables and data
    //-------------------------------------------------------------
    FileHandler fh;
    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
		const String in = getStringOption_("in");
		FileTypes::Type in_type = fh.getType(in);

    if (in_type==FileTypes::PEPXML)
  	{
  		String exp_name = getStringOption_("mz_file"),
				orig_name =	getStringOption_("mz_name");

			// no extension present => add one (will be removed by PepXMLFile)
			if (!orig_name.empty() && !orig_name.has('.')) {
				orig_name = orig_name + ".mzXML";
			}

  		protein_identifications.resize(1);
			if (exp_name.empty()) {
				PepXMLFile().load(in, protein_identifications[0],
													peptide_identifications, orig_name);
			}
			else {
				MSExperiment<> exp;
				fh.loadExperiment(exp_name, exp);
				if (!orig_name.empty()) {
					exp_name = orig_name;
				}
				PepXMLFile().load(in, protein_identifications[0],
													peptide_identifications, exp_name, exp);
			}
  	}
    else if ( in_type==FileTypes::IDXML)
  	{
  		IdXMLFile().load(in, protein_identifications, peptide_identifications);
  	}
    else if ( in_type==FileTypes::UNKNOWN && File::isDirectory(in) )
    {

      const String in_directory = File::absolutePath(in).ensureLastChar('/');
      const String mz_file = getStringOption_("mz_file");
      const bool ignore_proteins_per_peptide = getFlag_("ignore_proteins_per_peptide");

      UInt i = 0;
      FileHandler fh;
      FileTypes::Type type;
      MSExperiment<Peak1D> msexperiment;
      // Note: we had issues with leading zeroes, so let us represent scan numbers as Int (next line used to be map<String, Real> num_and_rt;)  However, now String::toInt() might throw.
      map<Int, Real> num_and_rt;
      vector<String> NativeID;

      // The mz-File (if given)
      if ( !mz_file.empty() )
      {
        type = fh.getTypeByFileName(mz_file);
        fh.loadExperiment(mz_file, msexperiment, type);

        for ( MSExperiment<Peak1D>::Iterator spectra_it = msexperiment.begin(); spectra_it != msexperiment.end(); ++spectra_it )
        {
          String(spectra_it->getNativeID()).split('=', NativeID);
          try
          {
            num_and_rt[NativeID[1].toInt()] = spectra_it->getRT();
            // std::cout << "num_and_rt: " << NativeID[1] << " = " << NativeID[1].toInt() << " : " << num_and_rt[NativeID[1].toInt()] << std::endl; // CG debuggging 2009-07-01
          }
          catch ( Exception::ConversionError &e )
          {
            writeLog_(String("Error: Cannot read scan number as integer. '") + e.getMessage());
          }
        }
      }

      // Get list of the actual Sequest .out-Files
      StringList in_files;
      if ( !File::fileList(in_directory, String("*.out"), in_files) )
      {
        writeLog_(String("Error: No .out files found in '") + in_directory + "'. Aborting!");
      }

      // Now get to work ...
      for ( vector<String>::const_iterator in_files_it = in_files.begin(); in_files_it != in_files.end(); ++in_files_it )
      {
        vector<PeptideIdentification> peptide_ids_seq;
        ProteinIdentification protein_id_seq;
        vector<DoubleReal> pvalues_seq;

        vector<String> in_file_vec;

        SequestOutfile sequest_outfile;

        writeDebug_(String("Reading file ") + *in_files_it, 3);

        try
        {
          sequest_outfile.load((String) (in_directory + *in_files_it), peptide_ids_seq, protein_id_seq, 1.0, pvalues_seq, "Sequest",
              ignore_proteins_per_peptide);

          in_files_it->split('.', in_file_vec);

          for ( Size j = 0; j < peptide_ids_seq.size(); ++j )
          {

            // We have to explicitly set the identifiers, because the normal set ones are composed of search engine name and date, which is the same for a bunch of sequest out-files.
            peptide_ids_seq[j].setIdentifier(*in_files_it + "_" + i);

            Int scan_number = 0;
            if ( !mz_file.empty() )
            {
              try
              {
                scan_number = in_file_vec.at(2).toInt();
                peptide_ids_seq[j].setMetaValue("RT", num_and_rt[scan_number]);
              }
              catch ( Exception::ConversionError &e )
              {
                writeLog_(String("Error: Cannot read scan number as integer. '") + e.getMessage());
              }
              catch ( std::exception &e )
              {
                writeLog_(String("Error: Cannot read scan number as integer. '") + e.what());
              }
              //	DoubleReal real_mz = ( (DoubleReal)peptide_ids_seq[j].getMetaValue("MZ") - hydrogen_mass )/ (DoubleReal)peptide_ids_seq[j].getHits()[0].getCharge(); // ???? semantics of mz
              const DoubleReal real_mz = (DoubleReal) peptide_ids_seq[j].getMetaValue("MZ") / (DoubleReal) peptide_ids_seq[j].getHits()[0].getCharge();
              peptide_ids_seq[j].setMetaValue("MZ", real_mz);
            }

            writeDebug_(String("scan: ") + String(scan_number) + String("  RT: ") + String(peptide_ids_seq[j].getMetaValue("RT")) + "  MZ: " + String(
                peptide_ids_seq[j].getMetaValue("MZ")) + "  Ident: " + peptide_ids_seq[j].getIdentifier(), 4);

            peptide_identifications.push_back(peptide_ids_seq[j]);
          }

          protein_id_seq.setIdentifier(*in_files_it + "_" + i);
          protein_identifications.push_back(protein_id_seq);
          ++i;
        }
        catch ( Exception::ParseError & pe )
        {
          writeLog_(pe.getMessage() + String("(file: ") + *in_files_it + ")");
          throw ;
        }
        catch (...)
        {
          writeLog_(String("Error reading file: ") + *in_files_it);
          throw;
        }
      }

      writeDebug_("All files processed.", 3);

    }
    else
    {
      writeLog_("Unknown input file type given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    const String out = getStringOption_("out");
    FileTypes::Type out_type = fh.nameToType(getStringOption_("out_type"));
    if (out_type==FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(out);
    }
    if (out_type==FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    if (out_type==FileTypes::PEPXML)
    {
      PepXMLFile().store(out, protein_identifications, peptide_identifications);
    }
    else if (out_type==FileTypes::IDXML)
    {
      IdXMLFile().store(out, protein_identifications, peptide_identifications);
    }
    else
    {
      writeLog_("Unknown output file type given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    return EXECUTION_OK;
  }
};

int
main(int argc, const char** argv)
{
  TOPPIDFileConverter tool;
  return tool.main(argc, argv);
}
