// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors:  Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
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
   <td ALIGN = "center" BGCOLOR="#EBEBEB" ROWSPAN=1> pot. predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderMRM \f$ \longrightarrow \f$</td>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
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

 Specialized tools are available for some experimental techniques: @ref TOPP_SILACAnalyzer, @ref TOPP_ITRAQAnalyzer.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_FeatureFinderMRM.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_FeatureFinderMRM.html

 For the parameters of the algorithm section see the algorithms documentation: @n

 @ref OpenMS::FeatureFinderAlgorithmMRM FeatureFinderAlgorithmMRM

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMRM
  : public TOPPBase
{
public:
  TOPPFeatureFinderMRM()
    : TOPPBase("FeatureFinderMRM", "Detects two-dimensional features in LC-MS data.")
  {}

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", StringList::create("mzML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", StringList::create("featureXML"));
    addEmptyLine_();
    addText_("All other options of the FeatureFinder are set in the 'algorithm' section of the INI file.\n");

    registerSubsection_("algorithm", "Algorithm section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    return FeatureFinder().getParameters(FeatureFinderAlgorithmMRM<Peak1D, Feature>::getProductName());
  }

  ExitCodes main_(int, const char **)
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
    MzMLFile f;
    f.setLogType(log_type_);

    f.load(in, exp);

    //no seeds supported
    FeatureMap<> seeds;

    //prevent loading of everything except MRM MS/MS spectra
    //exp.erase(remove_if(exp.begin(), exp.end(), HasScanMode<PeakMap::SpectrumType>(InstrumentSettings::SRM, true)), exp.end());
    // erase the spectra, we just need the chromatograms for the feature finder
    exp.erase(exp.begin(), exp.end());

    // A map for the resulting features
    FeatureMap<> features;

    // Apply the feature finder
    ff.run(FeatureFinderAlgorithmMRM<Peak1D, Feature>::getProductName(), exp, features, feafi_param, seeds);
    features.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    // DEBUG
    if (debug_level_ > 10)
    {
      FeatureMap<>::Iterator it;
      for (it = features.begin(); it != features.end(); ++it)
      {
        if (!it->isMetaEmpty())
        {
          vector<String> keys;
          it->getKeys(keys);
          LOG_INFO << "Feature " << it->getUniqueId() << endl;
          for (Size i = 0; i < keys.size(); i++)
          {
            LOG_INFO << "  " << keys[i] << " = " << it->getMetaValue(keys[i]) << endl;
          }
        }
      }
    }

    //-------------------------------------------------------------
    // writing files
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(features, getProcessingInfo_(DataProcessing::QUANTITATION));

    FeatureXMLFile map_file;
    map_file.store(out, features);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPFeatureFinderMRM tool;
  return tool.main(argc, argv);
}

/// @endcond
