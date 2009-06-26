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

	@brief Converts Sequest's .out-Files to .IdXML-Files.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDFileConverter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDFileConverter
	: public TOPPBase
{
public:
  TOPPIDFileConverter() :
    TOPPBase("IDFileConverter", "Converts identification engine output to an idXML file.", true)
  {
  }

protected:
  void
  registerOptionsAndFlags_()
  {

    registerStringOption_("in", "<dir>", "", "Input directory containing the .out-Files", true); // TODO description should be generic for all search engines. Is -in always a directory?
    registerStringOption_("out", "<file>", "", "Output file in IdXML-Format", true);

    registerStringOption_("type", "<type>", "", "Search engine", true);
    setValidStrings_("type", StringList::create("sequest")); // TODO implement support for mascot and more ...

    registerStringOption_("mz_file", "<file>", "", "Retention times will be looked up in this file, if supplied.\n[Note: Sequest .out files do not record retention times, but scan numbers.]", false);
    registerDoubleOption_("p_value", "<prob>", 1.0, "Filtering: Annotations with inferior p-value are ignored", false);

    registerFlag_("ignore_proteins_per_peptide", "Workaround to deal with .out files that contain e.g. \"+1\" in references column but do not list extra references in subsequent lines", true);
    // Please contact the maintainers if you know more about Sequest .out files and might help to resolve this issue
  }

  ExitCodes
  main_(int, const char**)
  {
    //-------------------------------------------------------------
    // variables
    //-------------------------------------------------------------

    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;
    const String outputfile_name = getStringOption_("out");

    const String type = getStringOption_("type");

    if ( type == "sequest" )
    {

      const String in_directory = File::absolutePath(getStringOption_("in")).ensureLastChar('/');
      const String mz_file = getStringOption_("mz_file");
      const DoubleReal p_value = getDoubleOption_("p_value");
      const bool ignore_proteins_per_peptide = getFlag_("ignore_proteins_per_peptide");

      UInt i = 0;
      FileHandler fh;
      FileTypes::Type type;
      MSExperiment<Peak1D> msexperiment;
      map<String, Real> num_and_rt;
      vector<String> NativeID;

      //-------------------------------------------------------------
      // reading input
      //-------------------------------------------------------------

      // The mz-File (if given)
      if ( !mz_file.empty() )
      {
        type = fh.getTypeByFileName(mz_file);
        fh.loadExperiment(mz_file, msexperiment, type);

        for ( MSExperiment<Peak1D>::Iterator spectra_it = msexperiment.begin(); spectra_it != msexperiment.end(); ++spectra_it )
        {
          String(spectra_it->getNativeID()).split('=', NativeID);
          num_and_rt[NativeID[1]] = spectra_it->getRT();

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
          sequest_outfile.load((String) (in_directory + *in_files_it), peptide_ids_seq, protein_id_seq, p_value, pvalues_seq, "Sequest",
              ignore_proteins_per_peptide);

          in_files_it->split('.', in_file_vec);

          for ( Size j = 0; j < peptide_ids_seq.size(); ++j )
          {

            // We have to explicitly set the identifiers, because the normal set ones are composed of search engine name and date, which is the same for a bunch of sequest out-files.
            peptide_ids_seq[j].setIdentifier(*in_files_it + "_" + i);

            writeDebug_(String("RT: ") + String(peptide_ids_seq[j].getMetaValue("RT")) + "  MZ: " + String(peptide_ids_seq[j].getMetaValue("MZ")) + "  Ident: "
                + peptide_ids_seq[j].getIdentifier(), 4);

            if ( !mz_file.empty() )
            {
              peptide_ids_seq[j].setMetaValue("RT", num_and_rt[in_file_vec[2]]);
              //	DoubleReal real_mz = ( (DoubleReal)peptide_ids_seq[j].getMetaValue("MZ") - hydrogen_mass )/ (DoubleReal)peptide_ids_seq[j].getHits()[0].getCharge(); // ???? semantics of mz
              const DoubleReal real_mz = (DoubleReal) peptide_ids_seq[j].getMetaValue("MZ") / (DoubleReal) peptide_ids_seq[j].getHits()[0].getCharge();
              peptide_ids_seq[j].setMetaValue("MZ", real_mz);
            }
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
      // Actually, this should be trapped far above when the valid strings are checked, but who knows...
      writeLog_(String("Support for -type ") + type + " is not yet implemented, sorry!");
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    IdXMLFile().store(outputfile_name, protein_identifications, peptide_identifications);

    return EXECUTION_OK;
  }
};

int
main(int argc, const char** argv)
{
  TOPPIDFileConverter tool;
  return tool.main(argc, argv);
}
