// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors:  Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>

#include <limits>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FeatureFinderCentroided FeatureFinderCentroided

@brief The feature detection application for quantitation (centroided).

<CENTER>
 <table>
  <tr>
   <th ALIGN = "center"> pot. predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=3> &rarr; FeatureFinderCentroided &rarr;</td>
   <th ALIGN = "center"> pot. successor tools </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_SeedListGenerator </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MapAlignerPoseClustering @n (or another alignment tool) </td>
  </tr>
 </table>
</CENTER>

Reference:\n
Weisser <em>et al.</em>: <a href="https://doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

This module identifies "features" in a LC/MS map. By feature, we understand a peptide in a MS sample that
reveals a characteristic isotope distribution. The algorithm
computes positions in rt and m/z dimension and a charge estimate
of each peptide.

The algorithm identifies pronounced regions of the data around so-called <tt>seeds</tt>.
In the next step, we iteratively fit a model of the isotope profile and the retention time to
these data points. Data points with a low probability under this model are removed from the
feature region. The intensity of the feature is then given by the sum of the data points included
in its regions.

How to find suitable parameters and details of the different algorithms implemented are described
in the "TOPP tutorial" (on https://openms.readthedocs.io/).

Specialized tools are available for some experimental techniques: @ref TOPP_IsobaricAnalyzer.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FeatureFinderCentroided.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FeatureFinderCentroided.html

For the parameters of the algorithm section see the algorithms documentation: @n
@ref OpenMS::FeatureFinderAlgorithmPicked "centroided" @n

In the following table you can find example values of the most important parameters for
different instrument types. @n These parameters are not valid for all instruments of that type,
but can be used as a starting point for finding suitable parameters.

<b>'centroided' algorithm</b>:
<table>
<tr>
 <td>&nbsp;</td>
 <td><b>Q-TOF</b></td>
 <td><b>LTQ Orbitrap</b></td>
</tr>
<tr>
 <td><b>intensity:bins</b></td>
 <td>10</td>
 <td>10</td>
</tr>
<tr>
 <td><b>mass_trace:mz_tolerance</b></td>
 <td>0.02</td>
 <td>0.004</td>
</tr>
<tr>
 <td><b>isotopic_pattern:mz_tolerance</b></td>
 <td>0.04</td>
 <td>0.005</td>
</tr>
</table>

For the @em centroided algorithm centroided data is needed. In order to create centroided data from profile data use the @ref TOPP_PeakPickerHiRes.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderCentroided :
  public TOPPBase
{
public:
  TOPPFeatureFinderCentroided() :
    TOPPBase("FeatureFinderCentroided", 
             "Detects two-dimensional features in LC-MS data.",
             true,
             {
               Citation{ "Sturm M",
                         "A novel feature detection algorithm for centroided data",
                         "Dissertation, 2010-09-15, p.37 ff",
                         "https://publikationen.uni-tuebingen.de/xmlui/bitstream/handle/10900/49453/pdf/Dissertation_Marc_Sturm.pdf"
                       }
             })
  {}

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerInputFile_("seeds", "<file>", "", "User specified seed list", false);
    setValidFormats_("seeds", ListUtils::create<String>("featureXML"));

    addEmptyLine_();

    registerSubsection_("algorithm", "Algorithm section");
  }


  Param getSubsectionDefaults_(const String& ) const override
  {
    return FeatureFinderAlgorithmPicked().getDefaultParameters();
  }


  ExitCodes main_(int, const char**) override
  {
    //input file names
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    // prevent loading of fragment spectra
    PeakFileOptions options;
    options.setMSLevels(vector<Int>(1, 1));

    // filter out zero (and negative) intensities
    using RP_TYPE = DRange<1>::PositionType;
    options.setIntensityRange({std::numeric_limits<RP_TYPE>::min(), RP_TYPE::maxPositive()});

    // reading input data
    FileHandler f;
    f.getOptions() = options;

    PeakMap exp;
    f.loadExperiment(in, exp, {FileTypes::MZML}, log_type_);
    exp.updateRanges();

    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS1 spectra in input file.");
    }

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType  spectrum_type = exp[0].getType();

    if (spectrum_type == SpectrumSettings::PROFILE)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    //load seeds
    FeatureMap seeds;
    if (!getStringOption_("seeds").empty())
    {
      FileHandler().loadFeatures(getStringOption_("seeds"), seeds, {FileTypes::FEATUREXML});
    }

    //setup of FeatureFinder
    FeatureFinderAlgorithmPicked ff;
    //ff.setLogType(log_type_); TODO

    // A map for the resulting features
    FeatureMap features;

    if (getFlag_("test"))
    {
      // if test mode set, add file without path so we can compare it
      features.setPrimaryMSRunPath({"file://" + File::basename(in)}, exp);
    }
    else
    {
      features.setPrimaryMSRunPath({in}, exp);
    }    
    
    // get parameters specific for the feature finder
    Param feafi_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3);

    // Apply the feature finder
    ff.run(exp, features, feafi_param, seeds);
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

    // write features to user specified output file
    FileHandler map_file;

    // Remove detailed convex hull information and subordinate features
    // (unless requested otherwise) to reduce file size of feature files
    // unless debugging is turned on.
    if (debug_level_ < 5)
    {
      for (Feature& ft : features)
      {
        ft.getConvexHull().expandToBoundingBox();
        for (Size i = 0; i < ft.getConvexHulls().size(); ++i)
        {
          ft.getConvexHulls()[i].expandToBoundingBox();
        }
        ft.getSubordinates().clear();
      }
    }

    map_file.storeFeatures(out, features, {FileTypes::FEATUREXML});

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderCentroided tool;
  return tool.main(argc, argv);
}

/// @endcond
