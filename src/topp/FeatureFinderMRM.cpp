// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors:  Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMRM.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FeatureFinderMRM FeatureFinderMRM

 @brief The feature detection application for quantitation.

<CENTER>
 <table>
  <tr>
   <th ALIGN = "center" ROWSPAN=1> pot. predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=3> &rarr; FeatureFinderMRM &rarr;</td>
   <th ALIGN = "center"> pot. successor tools </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=2>  </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MapAlignerPoseClustering @n (or another alignment tool) </td>
  </tr>
 </table>
</CENTER>

 This module identifies "features" in a LC/MS map. By feature, we understand a peptide in a MS sample that
 reveals a characteristic isotope distribution. The algorithm
 computes positions in rt and m/z dimension and a charge estimate
 of each peptide.

 How to find suitable parameters and details of the different algorithms implemented are described
 in the @ref TOPP_example_featuredetection "TOPP tutorial".

 Specialized tools are available for some experimental techniques: @ref TOPP_IsobaricAnalyzer.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_FeatureFinderMRM.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FeatureFinderMRM.html

 For the parameters of the algorithm section see the algorithms documentation: @n

 @ref OpenMS::FeatureFinderAlgorithmMRM FeatureFinderAlgorithmMRM

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMRM :
  public TOPPBase
{
public:
  TOPPFeatureFinderMRM() :
    TOPPBase("FeatureFinderMRM", "Detects two-dimensional features in LC-MS data.")
  {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return FeatureFinder().getParameters(FeatureFinderAlgorithmMRM::getProductName());
  }

  ExitCodes main_(int, const char**) override
  {
    //input file names and types
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    Param feafi_param = getParam_().copy("algorithm:", true);

    writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3);

    //setup of FeatureFinder
    FeatureFinder ff;
    ff.setLogType(log_type_);

    //reading input data
    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    //no seeds supported
    FeatureMap seeds;

    //prevent loading of everything except MRM MS/MS spectra
    //exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), HasScanMode<PeakMap::SpectrumType>(InstrumentSettings::SRM, true)), exp.end());
    // erase the spectra, we just need the chromatograms for the feature finder
    exp.getSpectra().erase(exp.begin(), exp.end());

    // A map for the resulting features
    FeatureMap features;

    if (getFlag_("test"))
    {
      // if test mode set, add file without path so we can compare it
      features.setPrimaryMSRunPath({"file://" + File::basename(in)});
    }
    else
    {
      features.setPrimaryMSRunPath({in}, exp);
    }

    // Apply the feature finder
    ff.run(FeatureFinderAlgorithmMRM::getProductName(), exp, features, feafi_param, seeds);
    features.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    // DEBUG
    if (debug_level_ > 10)
    {
      for (const Feature& ft : features)
      {
        if (!ft.isMetaEmpty())
        {
          vector<String> keys;
          ft.getKeys(keys);
          OPENMS_LOG_INFO << "Feature " << ft.getUniqueId() << endl;
          for (Size i = 0; i < keys.size(); i++)
          {
            OPENMS_LOG_INFO << "  " << keys[i] << " = " << ft.getMetaValue(keys[i]) << endl;
          }
        }
      }
    }

    //-------------------------------------------------------------
    // writing files
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(features, getProcessingInfo_(DataProcessing::QUANTITATION));

    FileHandler().storeFeatures(out, features, {FileTypes::FEATUREXML});

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderMRM tool;
  return tool.main(argc, argv);
}

/// @endcond
