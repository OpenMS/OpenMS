// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_PrecursorMassCorrector PrecursorMassCorrector

  @brief Corrects the precursor entries of MS/MS spectra, by using MS1 information.

  <center>
  <table>
  <tr>
  <th ALIGN = "center"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> &rarr; PrecursorMassCorrector &rarr;</td>
  <th ALIGN = "center"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
  </tr>
  </table>
  </center>

  @experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

  This tool corrects the m/z entries of MS/MS spectra by using MS1 information. Therefore,
  MS1 spectra must be supplied as profile mode spectra. The isotope distribution of the
  peptide in the MS1 level information are then used to determine the exact position
  of the monoisotopic peak. If no isotope distribution can be found the original
  entry is kept. As a side effect of determining the exact position of the monoisotopic
  peak is that the charge state is also annotated.

  This implementation uses the isotopewavelet featurefinder and sets the monoisotopic
  peak (and the charge) to the nearest feature.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_PrecursorMassCorrector.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_PrecursorMassCorrector.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPrecursorMassCorrector :
  public TOPPBase
{
public:
  TOPPPrecursorMassCorrector() :
    TOPPBase("PrecursorMassCorrector", "Corrects the precursor entries of MS/MS spectra, by using MS1 information.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file containing the spectra.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Output mzML file.");
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerInputFile_("feature_in", "<file>", "", "Input featureXML file, containing features; if set, the MS/MS spectra precursor entries \n"
                                                   "will be matched to the feature m/z values if possible.", false);
    setValidFormats_("feature_in", ListUtils::create<String>("featureXML"));

    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 1.5, "Maximal deviation in Th which is acceptable to be corrected;\n"
                                                                          "this value should be set to the instruments selection window.", false);
    setMinFloat_("precursor_mass_tolerance", 0);

    registerIntOption_("max_charge", "<charge>", 3, "Maximal charge that should be assumed for precursor peaks", false, true);
    registerDoubleOption_("intensity_threshold", "<threshold>", -1.0, "Intensity threshold value for isotope wavelet feature finder, please look at the documentation of the class for details.", false, true);
  }

  ExitCodes main_(int, const char**) override
  {
    // parsing parameters
    String in(getStringOption_("in"));
    String feature_in(getStringOption_("feature_in"));
    String out(getStringOption_("out"));
    double precursor_mass_tolerance(getDoubleOption_("precursor_mass_tolerance"));

    // reading input
    FileHandler fh;
    FileTypes::Type in_type = fh.getType(in);

    PeakMap exp;
    fh.loadExperiment(in, exp, in_type, log_type_, false, false);
    exp.sortSpectra();

    FeatureMap feature_map;
    if (!feature_in.empty())
    {
      FeatureXMLFile().load(feature_in, feature_map);
    }

    // calculations
    FeatureFinderAlgorithmIsotopeWavelet iso_ff;
    Param ff_param(iso_ff.getParameters());
    ff_param.setValue("max_charge", abs(getIntOption_("max_charge")));
    ff_param.setValue("intensity_threshold", getDoubleOption_("intensity_threshold"));
    iso_ff.setParameters(ff_param);

    FeatureFinder ff;
    ff.setLogType(ProgressLogger::NONE);

    PeakMap exp2 = exp;
    exp2.clear(false);
    for (const MSSpectrum& ms : exp)
    {
      if (!ms.empty())
      {
        exp2.addSpectrum(ms);
      }
    }

    exp = exp2;
    exp.updateRanges();

    // TODO check MS2 and MS1 counts
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, exp.size(), "Correcting precursor masses");
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      progresslogger.setProgress(exp.end() - it);
      if (it->getMSLevel() != 2)
      {
        continue;
      }
      // find first MS1 scan of the MS/MS scan
      PeakMap::Iterator ms1_it = it;
      while (ms1_it != exp.begin() && ms1_it->getMSLevel() != 1)
      {
        --ms1_it;
      }
      if (ms1_it == exp.begin() && ms1_it->getMSLevel() != 1)
      {
        writeLogWarn_("Did not find a MS1 scan to the MS/MS scan at RT=" + String(it->getRT()));
        continue;
      }
      if (ms1_it->empty())
      {
        writeDebug_("No peaks in scan at RT=" + String(ms1_it->getRT()) + String(", skipping"), 1);
        continue;
      }

      PeakMap::Iterator ms2_it = ms1_it;
      ++ms2_it;

      while (ms2_it != exp.end() && ms2_it->getMSLevel() == 2)
      {
        // first: error checks
        if (ms2_it->getPrecursors().empty())
        {
          writeDebug_("Warning: found no precursors of spectrum RT=" + String(ms2_it->getRT()) + ", skipping it.", 1);
          ++ms2_it;
          continue;
        }
        else if (ms2_it->getPrecursors().size() > 1)
        {
          writeLogWarn_("Warning: found more than one precursor of spectrum RT=" + String(ms2_it->getRT()) + ", using first one.");
        }

        Precursor prec = *ms2_it->getPrecursors().begin();
        double prec_pos = prec.getMZ();

        PeakMap new_exp;
        // now excise small region from the MS1 spec for the feature finder (isotope pattern must be covered...)
        PeakSpectrum zoom_spec;
        for (PeakSpectrum::ConstIterator pit = ms1_it->begin(); pit != ms1_it->end(); ++pit)
        {
          if (pit->getMZ() > prec_pos - 3 && pit->getMZ() < prec_pos + 3)
          {
            zoom_spec.push_back(*pit);
          }
        }
        new_exp.addSpectrum(zoom_spec);
        new_exp.updateRanges();
        FeatureMap features, seeds;
        ff.run("isotope_wavelet", new_exp, features, ff_param, seeds);
        if (features.empty())
        {
          writeDebug_("No features found for scan RT=" + String(ms1_it->getRT()), 1);
          ++ms2_it;
          continue;
        }

        double max_int(numeric_limits<double>::min());
        double min_dist(numeric_limits<double>::max());
        Size max_int_feat_idx(0);

        for (Size i = 0; i != features.size(); ++i)
        {
          if (fabs(features[i].getMZ() - prec_pos) < precursor_mass_tolerance &&
              features[i].getIntensity() > max_int)
          {
            max_int_feat_idx = i;
            max_int = features[i].getIntensity();
            min_dist = fabs(features[i].getMZ() - prec_pos);
          }
        }


        writeDebug_(" max_int=" + String(max_int) + " mz=" + String(features[max_int_feat_idx].getMZ()) + " charge=" + String(features[max_int_feat_idx].getCharge()), 5);
        if (min_dist < precursor_mass_tolerance)
        {
          prec.setMZ(features[max_int_feat_idx].getMZ());
          prec.setCharge(features[max_int_feat_idx].getCharge());
          vector<Precursor> precs;
          precs.push_back(prec);
          ms2_it->setPrecursors(precs);
          writeDebug_("Correcting precursor mass of spectrum RT=" + String(ms2_it->getRT()) + " from " + String(prec_pos) + " to " + String(prec.getMZ()) + " (z=" + String(prec.getCharge()) + ")", 1);
        }

        ++ms2_it;
      }
      it = --ms2_it;
    }
    progresslogger.endProgress();

    // writing output
    fh.storeExperiment(out, exp, log_type_);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPPrecursorMassCorrector tool;
  return tool.main(argc, argv);
}

/// @endcond
