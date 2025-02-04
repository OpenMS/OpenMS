// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_HighResPrecursorMassCorrector HighResPrecursorMassCorrector

@brief Corrects the precursor mz of high resolution data.

<CENTER>
<table>
 <tr>
   <th ALIGN = "center"> pot. predecessor tools </td>
   <td VALIGN= "middle" ROWSPAN=2> &rarr; HighResPrecursorMassCorrector &rarr;</td>
   <th ALIGN = "center"> pot. successor tools </td>
 </tr>
 <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
 </tr>
</table>
</CENTER>

This tool performs precursor m/z correction on picked (=centroided) high resolution data.

Three methods are available: 'nearest_peak', 'highest_intensity_peak' and 'feature'.
  - nearest_peak: Use nearest centroided MS1 peak for precursor mass correction.
  - highest_intensity_peak: Use highest intensity centroided MS1 peak in a given mass range for precursor mass correction.
  - feature: Use features for precursor mass correction, which allows for charge correction.

The method hightest_intensity_peak searches in a specific m/z-window of the precursor information for the peak with the highest intensity.
Suggestioned value 1/maximal expected charge. E.g maximal expected charge 5, m/z-window = +/- 0.2 Da

See the corresponding parameter subsection for details.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_HighResPrecursorMassCorrector.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_HighResPrecursorMassCorrector.html
*/

/// @cond TOPPCLASSES

#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>

class TOPPHiResPrecursorMassCorrector :
    public TOPPBase
{
  public:
    TOPPHiResPrecursorMassCorrector() :
      TOPPBase("HighResPrecursorMassCorrector", 
      "Corrects the precursor mass and charge determined by the instrument software.")
    {
    }

  protected:
    void registerOptionsAndFlags_() override
    {
      // input files
      registerInputFile_("in", "<file>", "", "Input file (centroided data)");
      setValidFormats_("in", ListUtils::create<String>("mzML"));

      registerOutputFile_("out", "<file>", "", "Output file");
      setValidFormats_("out", ListUtils::create<String>("mzML"));
      
      registerTOPPSubsection_("feature", "Use features for precursor mass correction.");
      registerInputFile_("feature:in", "<file>", "", "Features used to correct precursor masses.", false);
      setValidFormats_("feature:in", ListUtils::create<String>("featureXML"));
      registerDoubleOption_("feature:mz_tolerance", "<num>", 5.0, "The precursor mass tolerance. Used to determine matching to feature mass traces.", false);
      registerStringOption_("feature:mz_tolerance_unit", "<choice>", "ppm", "Unit of precursor mass tolerance", false);
      setValidStrings_("feature:mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));
      registerDoubleOption_("feature:rt_tolerance", "<num>", 0.0, "Additional retention time tolerance added to feature boundaries.", false);
      registerIntOption_("feature:max_trace", "<num>", 2, "Maximum isotopic trace considered in matching a precursor to a feature.", false, true);
      registerFlag_("feature:believe_charge", "Assume precursor charge to be correct.");
      registerFlag_("feature:keep_original", "Make a copy of the precursor and MS2 (true) or discard the original (false).");
      registerFlag_("feature:assign_all_matching", "Correct a precursor using all matching features (true) or only the nearest (false). Only evaluated if copies are created (feature:keep_original).");

      registerTOPPSubsection_("nearest_peak", "Use nearest centroided MS1 peak for precursor mass correction.");
      registerDoubleOption_("nearest_peak:mz_tolerance", "<num>", 0.0, "The precursor mass tolerance to find the closest MS1 peak. (Disable method by setting value to 0.0)", false);
      registerStringOption_("nearest_peak:mz_tolerance_unit", "<choice>", "ppm", "Unit of precursor mass tolerance", false);
      setValidStrings_("nearest_peak:mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));

      registerTOPPSubsection_("highest_intensity_peak", "Use centroided MS1 peak with the highest intensity in a certrain mass range - for precursor mass correction");
      registerDoubleOption_("highest_intensity_peak:mz_tolerance", "<num>", 0.0, "The precursor mass tolerance to find the highest intensity MS1 peak. Suggested value 1/max. expected charge. (Disable method by setting value to 0.0)", false);
      registerStringOption_("highest_intensity_peak:mz_tolerance_unit", "<choice>", "ppm", "Unit of precursor mass tolerance", false);
      setValidStrings_("highest_intensity_peak:mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));

      registerOutputFile_("out_csv", "<file>", "", "Optional CSV output file for results on 'nearest_peak' or 'highest_intensity_peak' algorithm (see corresponding subsection) containing columns: " + ListUtils::concatenate(ListUtils::create<String>(PrecursorCorrection::csv_header), ", ") + ".", false);
      setValidFormats_("out_csv", ListUtils::create<String>("csv"));
    }

    ExitCodes main_(int, const char **) override
    {
      const string in_mzml(getStringOption_("in"));
      const string in_feature(getStringOption_("feature:in"));
      const string out_mzml(getStringOption_("out"));
      const string out_csv = getStringOption_("out_csv");

      const double mz_tolerance = getDoubleOption_("feature:mz_tolerance");
      const bool mz_unit_ppm = getStringOption_("feature:mz_tolerance_unit") == "ppm" ? true : false;
      const double rt_tolerance = getDoubleOption_("feature:rt_tolerance");
      const int max_trace = getIntOption_("feature:max_trace");
      bool keep_original = getFlag_("feature:keep_original");
      bool assign_all_matching = getFlag_("feature:assign_all_matching");
      bool believe_charge = getFlag_("feature:believe_charge");

      const double nearest_peak_mz_tolerance = getDoubleOption_("nearest_peak:mz_tolerance");
      const bool nearest_peak_ppm = getStringOption_("nearest_peak:mz_tolerance_unit") == "ppm" ? true : false;

      const double highest_intensity_peak_mz_tolerance = getDoubleOption_("highest_intensity_peak:mz_tolerance");
      const bool highest_intensity_peak_ppm = getStringOption_("highest_intensity_peak:mz_tolerance_unit") == "ppm" ? true : false;

      PeakMap exp;
      FileHandler().loadExperiment(in_mzml, exp, {FileTypes::MZML});

      cout << setprecision(12);

      // determine accuracy
      vector<double> deltaMZs;
      vector<double> mzs;
      vector<double> rts;
      set<Size> corrected_precursors; // spectrum index of corrected precursors

      if ((nearest_peak_mz_tolerance <= 0.0) && (highest_intensity_peak_mz_tolerance <= 0.0) && in_feature.empty())
      {
        OPENMS_LOG_ERROR << "No method for PC correction requested. Either provide featureXML input files or set 'nearest_peak:mz_tolerance' > 0 or specify a 'highest_intensity_peak:mz_tolerance' > 0" << std::endl;
        return MISSING_PARAMETERS;
      }

      // perform correction to closest MS1 peak
      set<Size> corrected_to_nearest_peak;
      if (nearest_peak_mz_tolerance > 0.0 && highest_intensity_peak_mz_tolerance <= 0.0)
      {
        corrected_to_nearest_peak = PrecursorCorrection::correctToNearestMS1Peak(exp, nearest_peak_mz_tolerance, nearest_peak_ppm, deltaMZs, mzs, rts);
      }

      //perform correction to highest intensity MS1 peak
      set<Size> corrected_to_highest_intensity_peak;
      if (highest_intensity_peak_mz_tolerance > 0.0)
      {
        corrected_to_highest_intensity_peak = PrecursorCorrection::correctToHighestIntensityMS1Peak(exp, highest_intensity_peak_mz_tolerance, highest_intensity_peak_ppm, deltaMZs, mzs, rts);
      }
 
      // perform correction to closest feature (also corrects charge if not disabled)
      set<Size> corrected_to_nearest_feature;      
      if (!in_feature.empty())
      {
        FeatureMap features;
        FileHandler().loadFeatures(in_feature, features);
        corrected_to_nearest_feature = PrecursorCorrection::correctToNearestFeature(features, exp, rt_tolerance, mz_tolerance, mz_unit_ppm, believe_charge, keep_original, assign_all_matching, max_trace, debug_level_);
        corrected_precursors.insert(corrected_to_nearest_feature.begin(), corrected_to_nearest_feature.end());
      }

      FileHandler().storeExperiment(out_mzml, exp, {FileTypes::MZML},log_type_);

      if (!out_csv.empty())
      {
        if (nearest_peak_mz_tolerance > 0.0 && highest_intensity_peak_mz_tolerance <= 0.0)
        {
          OPENMS_LOG_INFO << "Corrected " << corrected_to_nearest_peak.size() << " precursor to a MS1 peak." << endl;
        }
        else if (highest_intensity_peak_mz_tolerance > 0.0)
        {
          OPENMS_LOG_INFO << "Corrected " << corrected_to_highest_intensity_peak.size() << " precursor to a MS1 peak." << endl;
        }
        else
        {
          OPENMS_LOG_WARN << "Output file 'out_csv': No data collected since 'nearest_peak:mz_tolerance' was not enabled. CSV will be empty." << endl;
        }
        PrecursorCorrection::writeHist(out_csv, deltaMZs, mzs, rts);
      }

      if (!in_feature.empty())
      {
        OPENMS_LOG_INFO << "Corrected " << corrected_to_nearest_feature.size() << " precursors to a feature." << endl;
      }

      return EXECUTION_OK;
    }

};

int main(int argc, const char ** argv)
{
  TOPPHiResPrecursorMassCorrector tool;
  return tool.main(argc, argv);
}

/// @endcond

