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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_ConvertTSVToTraML ConvertTSVToTraML

  @brief Converts special TSV files to TraML files

  The TSV files needs to have the following headers, all fields need to be separated by tabs:

  PrecursorMz (float)
  ProductMz (float)
  Tr_calibrated (float)
  transition_name (free text, needs to be unique for each transition [in this file])
  CE (float)
  LibraryIntensity (float)
  transition_group_id (free text, designates the transition group [e.g. peptide] to which this transition belongs)
  decoy (1==decoy, 0== no decoy; determines whether the transition is a decoy transition or not)
  PeptideSequence  (free text, sequence only (no modifications) )
  ProteinName  (free text)
  Annotation  (free text, e.g. y7)
  FullPeptideName  (free text, should contain modifications*)  
  MissedCleavages
  Replicates
  NrModifications
  Charge (integer)
  Labelgroup (free text, e.g. heavy or light)

* modifications should be supplied inside the sequence using UniMod
  identifiers or freetext identifiers that are understood by OpenMS. Please do
  not use the ambigous bracket notation (e.g. PEPT[+80]IDE or PEPT[181]IDE)
  since this is ambigous and will NOT be interpreted correctly!
  example: PEPT(Phosphorylation)IDE(UniMod:27)A )

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPConvertTSVToTraML : public TOPPBase
{
public:

  TOPPConvertTSVToTraML() :
  TOPPBase("ConvertTSVToTraML", "Converts a csv into a TraML file", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "transition file ('csv')");

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", StringList::create("TraML"));

  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    const char * tr_file = in.c_str();

    TraMLFile traml;
    TargetedExperiment targeted_exp;

    TransitionTSVReader tsv_reader = TransitionTSVReader();
    std::cout << "Reading " << in << std::endl;
    tsv_reader.setLogType(log_type_);
    tsv_reader.convertTSVToTargetedExperiment(tr_file, targeted_exp);
    tsv_reader.validateTargetedExperiment(targeted_exp);

    std::cout << "Writing " << out << std::endl;
    traml.store(out, targeted_exp);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{

  TOPPConvertTSVToTraML tool;
  return tool.main(argc, argv);
}
