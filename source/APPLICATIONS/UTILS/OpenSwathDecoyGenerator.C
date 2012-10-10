// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
#include <OpenMS/FORMAT/TraMLFile.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_OpenSwathDecoyGenerator OpenSwathDecoyGenerator

  @brief Generates decoys according to different models for a specific TraML

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


TODO: could theoretical also produce an annotation in the TraML of what it thinks the ion is?



*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathDecoyGenerator
: public TOPPBase
{
public:

  TOPPOpenSwathDecoyGenerator() :
    TOPPBase("OpenSwathDecoyGenerator", "Generates decoys according to different models for a specific TraML", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file ('TraML')");

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", StringList::create("TraML"));

    registerStringOption_("method", "<type>", "shuffle", "decoy generation method ('shuffle','pseudo-reverse','reverse','shift')", false);
    registerDoubleOption_("identity_threshold", "<double>", 0.7, "identity threshold for the shuffle algorithm", false);
    registerIntOption_("max_attempts", "<int>", 10, "maximum attempts to lower the sequence identity between target and decoy for the shuffle algorithm", false);
    registerDoubleOption_("mz_threshold", "<double>", 0.8, "MZ threshold in Thomson", false);
    registerDoubleOption_("mz_shift", "<double>", 20, "MZ shift in Thomson for shift decoy method", false);
    registerStringOption_("decoy_tag", "<type>", "DECOY_", "decoy tag", false);
    registerIntOption_("min_transitions", "<int>", 2, "minimal number of transitions", false);
    registerIntOption_("max_transitions", "<int>", 6, "maximal number of transitions", false);
    registerFlag_("theoretical", "Set this flag if only annotated transitions should be used and be corrected to the theoretical mz.");
    registerFlag_("append", "Set this flag if non-decoy TraML should be appended to the output.");
  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String method = getStringOption_("method");
    DoubleReal identity_threshold = getDoubleOption_("identity_threshold");
    Int max_attempts = getIntOption_("max_attempts");
    DoubleReal mz_threshold = getDoubleOption_("mz_threshold");
    DoubleReal mz_shift = getDoubleOption_("mz_shift");
    String decoy_tag = getStringOption_("decoy_tag");
    Int min_transitions = getIntOption_("min_transitions");
    Int max_transitions = getIntOption_("max_transitions");
    bool theoretical = getFlag_("theoretical");
    bool append = getFlag_("append");

    if (method != "shuffle" && method != "pseudo-reverse" && method != "reverse" && method != "shift")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No valid decoy generation method selected!");
    }

    TraMLFile traml;
    TargetedExperiment targeted_exp;
    TargetedExperiment targeted_decoy;

    std::cout << "Loading " << in << std::endl;
    traml.load(in, targeted_exp);

    MRMDecoy decoys = MRMDecoy();

    std::cout << "Restricting transitions" << std::endl;
    decoys.restrictTransitions(targeted_exp, min_transitions, max_transitions);
    decoys.generateDecoys(targeted_exp, targeted_decoy, method, decoy_tag, identity_threshold, max_attempts, mz_threshold, theoretical, mz_shift);

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
