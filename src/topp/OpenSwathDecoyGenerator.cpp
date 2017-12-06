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
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/TraMLFile.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_OpenSwathDecoyGenerator OpenSwathDecoyGenerator

  @brief Generates decoys according to different models for a specific TraML

  <CENTER>
      <table>
          <tr>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ OpenSwathDecoyGenerator \f$ \longrightarrow \f$</td>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
          </tr>
      </table>
  </CENTER>

  This module generates "decoy" transitions from a set of real or "target"
  transitions. The idea is to use the decoy transitions in a statistical scoring
  process to estimate the false hits in an SRM / SWATH experiment.

  There are multiple methods to create the decoy transitions, the simplest ones
  are reverse and pseudo-reverse which reverse the sequence either completely or
  leaving the last (tryptic) AA untouched respectively.

  Another decoy generation method is "shuffle" which uses an algorithm similar
  to the one described in Lam, Henry, et al. (2010). "Artificial decoy spectral
  libraries for false discovery rate estimation in spectral library searching in
  proteomics".  Journal of Proteome Research 9, 605-610. It shuffles the amino
  acid sequence and shuffles the fragment ion intensities accordingly, however
  for this to work the fragment ions need to be matched to annotated before.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_OpenSwathDecoyGenerator.cli

  <B>The algorithm parameters for the Analyzer filter are:</B>
  @htmlinclude TOPP_OpenSwathDecoyGenerator.html


*/
// TODO: could theoretical also produce an annotation in the TraML of what it thinks the ion is?

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathDecoyGenerator
: public TOPPBase
{
public:

  TOPPOpenSwathDecoyGenerator() :
    TOPPBase("OpenSwathDecoyGenerator", "Generates decoys according to different models for a specific TraML", true)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ('traML')");
    setValidFormats_("in", ListUtils::create<String>("traML"));
    
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("traML"));

    registerStringOption_("method", "<type>", "shuffle", "decoy generation method ('shuffle','pseudo-reverse','reverse','shift')", false);
    registerStringOption_("decoy_tag", "<type>", "DECOY_", "decoy tag", false);
    registerDoubleOption_("mz_threshold", "<double>", 0.05, "MZ threshold in Thomson for fragment ion annotation", false);
    registerFlag_("exclude_similar", "set this flag if decoy assays with similarity of the peptide sequence to the target assays higher than the identity_threshold should be excluded. If similarity_threshold is over 0, decoy assays with an absolute difference of the decoy and target product mz smaller than similarity_threshold are further excluded.");
    registerDoubleOption_("similarity_threshold", "<double>", -1, "Similarity threshold for absolute difference of the product mz of target and decoy assays for exclusion in Dalton. Suggested value: 0.05", false);
    registerFlag_("append", "set this flag if non-decoy TraML should be appended to the output.");
    registerFlag_("remove_CNterm_mods", "set this flag to remove decoy peptides with C/N terminal modifications (may be necessary depending on the decoy generation method).");
    registerFlag_("remove_unannotated", "set this flag if target assays with unannotated ions should be ignored from decoy generation.");
    registerDoubleOption_("identity_threshold", "<double>", 0.7, "shuffle: identity threshold for the shuffle algorithm", false);
    registerIntOption_("max_attempts", "<int>", 10, "shuffle: maximum attempts to lower the sequence identity between target and decoy for the shuffle algorithm", false);
    registerDoubleOption_("mz_shift", "<double>", 20, "shift: MZ shift in Thomson for shift decoy method", false);
    registerDoubleOption_("precursor_mass_shift", "<double>", 0.0, "Mass shift to apply to the precursor ion", false);
    registerStringOption_("allowed_fragment_types", "<type>", "b,y", "allowed fragment types", false);
    registerStringOption_("allowed_fragment_charges", "<type>", "1,2,3,4", "allowed fragment charge states", false);
    registerFlag_("enable_detection_specific_losses", "set this flag if specific neutral losses for detection fragment ions should be allowed");
    registerFlag_("enable_detection_unspecific_losses", "set this flag if unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) for detection fragment ions should be allowed");

  }

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String method = getStringOption_("method");
    String decoy_tag = getStringOption_("decoy_tag");
    double mz_threshold = getDoubleOption_("mz_threshold");
    bool exclude_similar = getFlag_("exclude_similar");
    double similarity_threshold = getDoubleOption_("similarity_threshold");
    bool append = getFlag_("append");
    bool remove_CNterm_mods = getFlag_("remove_CNterm_mods");
    bool remove_unannotated = getFlag_("remove_unannotated");
    double identity_threshold = getDoubleOption_("identity_threshold");
    Int max_attempts = getIntOption_("max_attempts");
    double mz_shift = getDoubleOption_("mz_shift");
    double precursor_mass_shift = getDoubleOption_("precursor_mass_shift");
    String allowed_fragment_types_string = getStringOption_("allowed_fragment_types");
    String allowed_fragment_charges_string = getStringOption_("allowed_fragment_charges");
    bool enable_detection_specific_losses = getFlag_("enable_detection_specific_losses");
    bool enable_detection_unspecific_losses = getFlag_("enable_detection_unspecific_losses");

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

    if (method != "shuffle" && method != "pseudo-reverse" && method != "reverse" && method != "shift")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No valid decoy generation method selected!");
    }

    TraMLFile traml;
    TargetedExperiment targeted_exp;
    TargetedExperiment targeted_decoy;

    std::cout << "Loading " << in << std::endl;
    traml.load(in, targeted_exp);

    MRMDecoy decoys = MRMDecoy();

    std::cout << "Generate decoys" << std::endl;
    decoys.generateDecoys(targeted_exp, targeted_decoy, method, decoy_tag, identity_threshold, max_attempts, mz_threshold, mz_shift, exclude_similar, similarity_threshold, remove_CNterm_mods, precursor_mass_shift, allowed_fragment_types, allowed_fragment_charges, enable_detection_specific_losses, enable_detection_unspecific_losses, remove_unannotated);

    if (append)
    {
      TargetedExperiment targeted_merged;
      targeted_merged += targeted_exp + targeted_decoy;
      traml.store(out, targeted_merged);
    }
    else
    {
      traml.store(out, targeted_decoy);
    }
    return EXECUTION_OK;
  }

};

int main(int argc, const char **argv)
{
  TOPPOpenSwathDecoyGenerator gen;
  return gen.main(argc, argv);
}

/// @endcond
