// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <iostream>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_OpenSwathAssayGenerator OpenSwathAssayGenerator

  @brief Generates filtered and optimized assays using TraML files.

  <CENTER>
      <table>
          <tr>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ OpenSwathAssayGenerator \f$ \longrightarrow \f$</td>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathDecoyGenerator </td>
          </tr>
      </table>
  </CENTER>

  This module generates assays for targeted proteomics using a set of rules
  that was found to improve the sensitivity and selectivity for detection
  of typical peptides (Schubert et al., 2015). The tool operates on TraML
  files, which can come from TargetedFileConverter or any other tool. If the
  TraML is annotated with the CV terms for fragment ion annotation, it can
  directly filter the transitions according to the set rules. If this is not
  the case (e.g. if an older version of TargetedFileConverter was used), the
  option -enable_reannotation can do the reannotation.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_OpenSwathAssayGenerator.cli

  <B>The algorithm parameters for the Analyzer filter are:</B>
  @htmlinclude TOPP_OpenSwathAssayGenerator.html


*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathAssayGenerator :
  public TOPPBase
{
public:

  TOPPOpenSwathAssayGenerator() :
    TOPPBase("OpenSwathAssayGenerator", "Generates assays according to different models for a specific TraML", true)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ('traML')");
    setValidFormats_("in", ListUtils::create<String>("traML"));

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("traML"));

    registerIntOption_("min_transitions", "<int>", 6, "minimal number of transitions", false);
    registerIntOption_("max_transitions", "<int>", 6, "maximal number of transitions", false);
    registerStringOption_("allowed_fragment_types", "<type>", "b,y", "allowed fragment types", false);
    registerStringOption_("allowed_fragment_charges", "<type>", "1,2,3,4", "allowed fragment charge states", false);
    registerFlag_("enable_detection_specific_losses", "set this flag if specific neutral losses for detection fragment ions should be allowed");
    registerFlag_("enable_detection_unspecific_losses", "set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for detection fragment ions should be allowed");
    registerFlag_("enable_identification_specific_losses", "set this flag if specific neutral losses for identification fragment ions should be allowed");
    registerFlag_("enable_identification_unspecific_losses", "set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for identification fragment ions should be allowed");
    registerFlag_("enable_identification_ms2_precursors", "set this flag if MS2-level precursor ions for identification should be allowed to enable extraction of the precursor signal from the fragment ion data (MS2-level). This may help in identification if the MS1 signal is weak.");
    registerFlag_("enable_ms1_uis_scoring", "set this flag if MS1-UIS assays for UIS scoring should be generated");
    registerFlag_("enable_ms2_uis_scoring", "set this flag if MS2-UIS assays for UIS scoring should be generated");
    registerIntOption_("max_num_alternative_localizations", "<int>", 20, "maximum number of site-localization permutations", false);
    registerFlag_("enable_reannotation", "set this flag if reannotation of fragment ions should be allowed.");
    registerDoubleOption_("precursor_mz_threshold", "<double>", 0.025, "MZ threshold in Thomson for precursor ion selection", false);
    registerDoubleOption_("precursor_lower_mz_limit", "<double>", 400, "lower MZ limit for precursor ions", false);
    registerDoubleOption_("precursor_upper_mz_limit", "<double>", 1200, "upper MZ limit for precursor ions", false);
    registerDoubleOption_("product_mz_threshold", "<double>", 0.025, "MZ threshold in Thomson for fragment ion annotation", false);
    registerDoubleOption_("product_lower_mz_limit", "<double>", 350, "lower MZ limit for fragment ions", false);
    registerDoubleOption_("product_upper_mz_limit", "<double>", 2000, "upper MZ limit for fragment ions", false);

    registerInputFile_("swath_windows_file", "<file>", "", "Tab separated file containing the SWATH windows for exclusion of fragment ions falling into the precursor isolation window: lower_offset upper_offset \\newline 400 425 \\newline ... Note that the first line is a header and will be skipped.", false, true);
    setValidFormats_("swath_windows_file", ListUtils::create<String>("txt"));
  }

  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    Int min_transitions = getIntOption_("min_transitions");
    Int max_transitions = getIntOption_("max_transitions");
    String allowed_fragment_types_string = getStringOption_("allowed_fragment_types");
    String allowed_fragment_charges_string = getStringOption_("allowed_fragment_charges");
    bool enable_detection_specific_losses = getFlag_("enable_detection_specific_losses");
    bool enable_detection_unspecific_losses = getFlag_("enable_detection_unspecific_losses");
    bool enable_identification_specific_losses = getFlag_("enable_identification_specific_losses");
    bool enable_identification_unspecific_losses = getFlag_("enable_identification_unspecific_losses");
    bool enable_identification_ms2_precursors = getFlag_("enable_identification_ms2_precursors");
    bool enable_ms1_uis_scoring = getFlag_("enable_ms1_uis_scoring");
    bool enable_ms2_uis_scoring = getFlag_("enable_ms2_uis_scoring");
    size_t max_num_alternative_localizations = getIntOption_("max_num_alternative_localizations");
    double precursor_mz_threshold = getDoubleOption_("precursor_mz_threshold");
    double precursor_lower_mz_limit = getDoubleOption_("precursor_lower_mz_limit");
    double precursor_upper_mz_limit = getDoubleOption_("precursor_upper_mz_limit");
    double product_mz_threshold = getDoubleOption_("product_mz_threshold");
    double product_lower_mz_limit = getDoubleOption_("product_lower_mz_limit");
    double product_upper_mz_limit = getDoubleOption_("product_upper_mz_limit");
    String swath_windows_file = getStringOption_("swath_windows_file");
    bool enable_reannotation = getFlag_("enable_reannotation");

    std::vector<String> allowed_fragment_types;
    allowed_fragment_types_string.split(",", allowed_fragment_types);

    std::vector<String> allowed_fragment_charges_string_vector;
    std::vector<size_t> allowed_fragment_charges;
    allowed_fragment_charges_string.split(",", allowed_fragment_charges_string_vector);
    for (size_t i = 0; i < allowed_fragment_charges_string_vector.size(); i++)
    {
      size_t charge = std::atoi(allowed_fragment_charges_string_vector.at(i).c_str());
      allowed_fragment_charges.push_back(charge);
    }

    std::vector<std::pair<double, double> > swathes;
    // Check swath window input
    if (!swath_windows_file.empty())
    {
      LOG_INFO << "Validate provided Swath windows file:" << std::endl;
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);

      LOG_INFO << "Read Swath maps file with " << swath_prec_lower.size() << " windows." << std::endl;
      for (Size i = 0; i < swath_prec_lower.size(); i++)
      {
        swathes.push_back(std::make_pair(swath_prec_lower[i], swath_prec_upper[i]));
        LOG_DEBUG << "Read lower swath window " << swath_prec_lower[i] << " and upper window " << swath_prec_upper[i] << std::endl;
      }
    }

    TraMLFile traml;
    TargetedExperiment targeted_exp;

    LOG_INFO << "Loading " << in << std::endl;
    traml.load(in, targeted_exp);

    MRMAssay assays = MRMAssay();
    assays.setLogType(ProgressLogger::CMD);

    LOG_INFO << "Annotating transitions" << std::endl;
    assays.reannotateTransitions(targeted_exp, precursor_mz_threshold, product_mz_threshold, allowed_fragment_types, allowed_fragment_charges, enable_reannotation, enable_detection_specific_losses, enable_detection_unspecific_losses);

    LOG_INFO << "Annotating detecting transitions" << std::endl;
    assays.restrictTransitions(targeted_exp, product_lower_mz_limit, product_upper_mz_limit, swathes);
    assays.detectingTransitions(targeted_exp, min_transitions, max_transitions);

    if (enable_ms1_uis_scoring || enable_ms2_uis_scoring)
    {
      std::vector<std::pair<double, double> > uis_swathes;

      if (enable_ms1_uis_scoring)
      {
        int num_precursor_windows = static_cast<int>(Math::round((precursor_upper_mz_limit - precursor_lower_mz_limit) / precursor_mz_threshold));
        for (int i = 0; i < num_precursor_windows; i++)
        {
          uis_swathes.push_back(std::make_pair((precursor_lower_mz_limit+(i*precursor_mz_threshold)),(precursor_lower_mz_limit+((i+1)*precursor_mz_threshold))));
        }
      }
      else {uis_swathes = swathes;}
      
      std::cout << "Generating identifying (UIS) transitions" << std::endl;
      assays.uisTransitions(targeted_exp, allowed_fragment_types, allowed_fragment_charges, enable_identification_specific_losses, enable_identification_unspecific_losses, enable_identification_ms2_precursors, product_mz_threshold, uis_swathes, -4, max_num_alternative_localizations, -1);
      std::vector<std::pair<double, double> > empty_swathes;
      assays.restrictTransitions(targeted_exp, product_lower_mz_limit, product_upper_mz_limit, empty_swathes);
    }

    std::cout << "Writing assays " << out << std::endl;
    traml.store(out, targeted_exp);
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPOpenSwathAssayGenerator gen;
  return gen.main(argc, argv);
}

/// @endcond
