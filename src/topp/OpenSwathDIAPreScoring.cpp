// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataFrameWriter.h>

#include <OpenMS/FORMAT/FileHandler.h>

#include <OpenMS/SYSTEM/File.h>

#include <iostream>


/**
@page TOPP_OpenSwathDIAPreScoring OpenSwathDIAPreScoring

@brief ...

SWATH specific parameters only apply if you have full MS2 spectra maps.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathDIAPreScoring.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathDIAPreScoring.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
using namespace OpenMS;

class DIAPreScoring :
  public TOPPBase
{
public:

  DIAPreScoring() :
    TOPPBase("OpenSwathDIAPreScoring", "Scoring spectra using the DIA scores.")
  {
  }

protected:

  typedef PeakMap MapType;
  typedef boost::shared_ptr<PeakMap> MapTypePtr;

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("tr", "<file>", "", "transition file");
    setValidFormats_("tr", ListUtils::create<String>("traML"));
    registerInputFileList_("swath_files", "<files>", StringList(),
                           "Swath files that were used to extract the transitions. If present, SWATH specific scoring will be applied.",
                           true);
    setValidFormats_("swath_files", ListUtils::create<String>("mzML"));
    registerOutputFileList_("output_files", "<files>", StringList(),
                           "Output files. One per Swath input file.",
                           true);
    setValidFormats_("output_files", ListUtils::create<String>("tsv"));

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0,
                          "Minimal distance to the edge to still consider a precursor, in Thomson (only in SWATH)",
                          false);
  }

  Param getSubsectionDefaults_(const String&) const override
  {
    return OpenMS::DiaPrescore().getDefaults();
  }

  ExitCodes main_(int, const char**) override
  {
    OpenMS::StringList file_list = getStringList_("swath_files");
    OpenMS::StringList outfile_list = getStringList_("output_files");
    std::string tr_file = getStringOption_("tr");
    std::cout << tr_file << std::endl;
    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");

    // If we have a transformation file, trafo will transform the RT in the
    // scoring according to the model. If we don't have one, it will apply the
    // null transformation.
    Param feature_finder_param = getParam_().copy("algorithm:", true);

    // Create the output map, load the input TraML file and the chromatograms
    MapType exp;
    OpenSwath::LightTargetedExperiment transition_exp;

    std::cout << "Loading TraML file" << std::endl;
    {
      OpenMS::TargetedExperiment transition_exp_;
      FileHandler().loadTransitions(tr_file, transition_exp_, {FileTypes::TRAML});
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transition_exp);
      int ltrans = transition_exp.transitions.size();
      std::cout << ltrans << std::endl;
    }
    // Here we deal with SWATH files (can be multiple files)

    for (Size i = 0; i < file_list.size(); ++i)
    {
      MapTypePtr swath_map (new MapType);
      FeatureMap featureFile;
      std::cout << "Loading file " << file_list[i] << std::endl;
      String fname = outfile_list[i];
      FileHandler().loadExperiment(file_list[i], *swath_map, {FileTypes::MZML}, log_type_);
      if (swath_map->empty() || (*swath_map)[0].getPrecursors().empty())
      {
        std::cerr << "WARNING: File " << swath_map->getLoadedFilePath()
                  << " does not have any experiments or any precursors. Is it a SWATH map?"
                  << std::endl;
        continue;
      }
      // Find the transitions to extract and extract them
      OpenSwath::LightTargetedExperiment transition_exp_used;
      double upper, lower;
      const std::vector<Precursor> prec = (*swath_map)[0].getPrecursors();
      lower = prec[0].getMZ() - prec[0].getIsolationWindowLowerOffset();
      upper = prec[0].getMZ() + prec[0].getIsolationWindowUpperOffset();
      OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used,
                                              min_upper_edge_dist, lower, upper);
      if (transition_exp_used.getTransitions().empty())
      {
        std::cerr << "WARNING: For file " << swath_map->getLoadedFilePath()
                  << " there are no transitions to extract." << std::endl;
        continue;
      }
      std::cout << "Using Spectrum Interface!" << std::endl;
      OpenSwath::SpectrumAccessPtr  spectrumAccess = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(
        swath_map);
      OpenSwath::IDataFrameWriter* dfw = new OpenSwath::CSVWriter(fname);
      OpenMS::DiaPrescore dp;
      OpenMS::RangeMobility im_range; // create empty IM range object
      dp.operator()(spectrumAccess, transition_exp_used, im_range, dfw); //note IM not supported here yet
      delete dfw;
    }         //end of for loop
    return EXECUTION_OK;
  }       //end of _main

};

int main(int argc, const char** argv)
{
  DIAPreScoring tool;
  int code = tool.main(argc, argv);
  return code;

}

/// @endcond
