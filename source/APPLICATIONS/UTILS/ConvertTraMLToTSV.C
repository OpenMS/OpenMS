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
  @page TOPP_ConvertTraMLToTSV ConvertTraMLToTSV

  @brief Converts a TraML file to a csv file
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPConvertTraMLToTSV : public TOPPBase
{
public:

  TOPPConvertTraMLToTSV() :
  TOPPBase("ConvertTraMLToTSV", "Converts a TraML file into TSV", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "transition file");
    setValidFormats_("in", StringList::create("TraML"));

    registerOutputFile_("out", "<file>", "", "output file csv");
  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    const char * tr_file = out.c_str();

    TraMLFile traml;
    TargetedExperiment targeted_exp;

    TransitionTSVReader tsv_reader = TransitionTSVReader();
    tsv_reader.setLogType(log_type_);
    traml.load(in, targeted_exp);
    tsv_reader.convertTargetedExperimentToTSV(tr_file, targeted_exp);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{

  TOPPConvertTraMLToTSV tool;
  return tool.main(argc, argv);
}
