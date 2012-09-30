// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
//#include "ConvertTSVToTraML.h"
//#include "ChromatogramExtractor.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/TraMLFile.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MRMDecoyGenerator MRMDecoyGenerator

  @brief Generates decoys according to different models for a specific TraML

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPMRMDecoyGenerator
  : public TOPPBase
{
 public:

  TOPPMRMDecoyGenerator()
    : TOPPBase("MRMDecoyGenerator","Generates decoys according to different models for a specific TraML", false)
  {
  }

 protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in","<file>","","input file ('TraML')");

    registerOutputFile_("out","<file>","","output file");
    setValidFormats_("out",StringList::create("TraML"));

    registerStringOption_("method","<type>","shuffle","decoy generation method ('shuffle','reverse','trypticreverse')", false);
    registerDoubleOption_("identity_threshold","<double>",0.7,"identity threshold", false);
    registerDoubleOption_("mz_threshold","<double>",0.8,"MZ threshold in Thomson", false);
    registerStringOption_("decoy_tag","<type>","DECOY_","decoy tag", false);
    registerIntOption_("min_transitions","<int>",2,"minimal number of transitions",false);
    registerIntOption_("max_transitions","<int>",6,"maximal number of transitions",false);
    registerFlag_("theoretical","Set this flag if only annotated transitions should be used and be corrected to the theoretical mz.");
    registerFlag_("append","Set this flag if non-decoy TraML should be appended to the output.");
  }

  ExitCodes main_(int , const char**)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String method = getStringOption_("method");
    DoubleReal identity_threshold = getDoubleOption_("identity_threshold");
    DoubleReal mz_threshold = getDoubleOption_("mz_threshold");
    String decoy_tag = getStringOption_("decoy_tag");
    Int min_transitions = getIntOption_("min_transitions");
    Int max_transitions = getIntOption_("max_transitions");
    bool theoretical = getFlag_("theoretical");
    bool append = getFlag_("append");

    TraMLFile traml;
    TargetedExperiment targeted_exp;
    TargetedExperiment targeted_decoy;
        
    std::cout << "Loading " << in << std::endl;
    traml.load(in, targeted_exp);
    
    MRMDecoy decoys = MRMDecoy();

    std::cout << "Restricting transitions" << std::endl;
    decoys.restrictTransitions(targeted_exp, min_transitions, max_transitions);
    decoys.generateDecoys(targeted_exp, targeted_decoy, method, decoy_tag, identity_threshold, mz_threshold, theoretical);
    
    if (append)
    {
       TargetedExperiment targeted_merged;
       targeted_merged += targeted_exp + targeted_decoy;
       //TransitionTSVReader tsv_reader = TransitionTSVReader();
       //tsv_reader.validateTargetedExperiment(targeted_merged);
       traml.store(out, targeted_merged);
    }
    else
    {
       //TransitionTSVReader tsv_reader = TransitionTSVReader();
       //tsv_reader.validateTargetedExperiment(targeted_decoy);
       traml.store(out, targeted_decoy);
    }

    return EXECUTION_OK;
  }

};

int main( int argc, const char** argv )
{
  TOPPMRMDecoyGenerator gen;
  return gen.main(argc,argv);
}

