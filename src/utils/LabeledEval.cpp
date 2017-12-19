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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_LabeledEval LabeledEval

    @brief Evaluation tool for isotope-labeled quantitation experiments.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_LabeledEval.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_LabeledEval.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPLabeledEval :
  public TOPPBase
{
public:
  TOPPLabeledEval() :
    TOPPBase("LabeledEval", " Evaluation tool for isotope-labeled quantitation experiments.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Feature result file");
    setValidFormats_("in", ListUtils::create<String>("featureXML"));
    registerInputFile_("truth", "<file>", "", "Expected result file.");
    setValidFormats_("truth", ListUtils::create<String>("consensusXML"));
    registerDoubleOption_("rt_tol", "<tol>", 20.0, "Maximum allowed retention time deviation", false);
    registerDoubleOption_("mz_tol", "<tol>", 0.25, "Maximum allowed m/z deviation (divided by charge)", false);
  }

  String fiveNumbers(vector<double> a, UInt decimal_places)
  {
    sort(a.begin(), a.end());
    return String::number(a[0], decimal_places) + " " + String::number(a[a.size() / 4], decimal_places) + " " + String::number(a[a.size() / 2], decimal_places) + " " + String::number(a[(3 * a.size()) / 4], decimal_places) + " " + String::number(a.back(), decimal_places);
  }

  String fiveNumberQuotients(vector<double> a, vector<double> b, UInt decimal_places)
  {
    vector<double> errors;
    for (Size i = 0; i < a.size(); ++i) errors.push_back(a[i] / b[i]);
    return fiveNumbers(errors, decimal_places);
  }

  ExitCodes main_(int, const char **) override
  {
    //load input features
    FeatureMap input;
    FeatureXMLFile().load(getStringOption_("in"), input);

    //load truth consensusXML
    ConsensusMap truth;
    ConsensusXMLFile().load(getStringOption_("truth"), truth);

    //parameters
    double mz_tol = getDoubleOption_("mz_tol");
    double rt_tol = getDoubleOption_("rt_tol");

    //seek manual feature in automatic feature map
    UInt matched_pairs = 0;
    UInt half_matched_pairs = 0;
    vector<double> t_ratio, i_ratio, rt_diffs, mz_diffs;
    for (Size t = 0; t < truth.size(); ++t)
    {
      if (truth[t].size() != 2)
      {
        cerr << "Error: consensus feature must contain exactly two elements!" << endl;
        continue;
      }
      vector<Feature> best_matches(2);
      vector<UInt> match_counts(2, 0);
      vector<Peak2D> elements(2);
      elements[0] = *(truth[t].getFeatures().begin());
      elements[1] = *(++(truth[t].getFeatures().begin()));
      double mz_tol_charged = mz_tol / truth[t].getCharge();
      for (Size e = 0; e < 2; ++e)
      {
        double best_score = 0.0;
        for (Size i = 0; i < input.size(); ++i)
        {
          const Feature & f_i = input[i];
          if (fabs(f_i.getRT() - elements[e].getRT()) < rt_tol
             && fabs(f_i.getMZ() - elements[e].getMZ()) < mz_tol_charged)
          {
            ++match_counts[e];
            double score = (1.0 - fabs(f_i.getMZ() - elements[e].getMZ()) / mz_tol_charged) * (1.0 - fabs(f_i.getRT() - elements[e].getRT()) / rt_tol);
            if (score > best_score)
            {
              best_score = score;
              best_matches[e] = f_i;
            }
          }
        }
      }

      //not matched
      if (match_counts[0] == 0 && match_counts[1] == 0)
      {
      }
      //half matched
      else if ((match_counts[0] > 0 && match_counts[1] == 0) || (match_counts[0] == 0 && match_counts[1] > 0))
      {
        ++half_matched_pairs;
      }
      //matched
      else
      {
        ++matched_pairs;
        double a_r = best_matches[0].getIntensity() / best_matches[1].getIntensity();
        t_ratio.push_back(a_r);
        double m_r = elements[0].getIntensity() / elements[1].getIntensity();
        i_ratio.push_back(m_r);
        rt_diffs.push_back(best_matches[1].getRT() - best_matches[0].getRT());
        mz_diffs.push_back((best_matches[1].getMZ() - best_matches[0].getMZ()) * truth[t].getCharge());
      }
    }

    cout << endl;
    cout << "pair detection statistics:" << endl;
    cout << "==========================" << endl;
    cout << "truth pairs: " << truth.size() << endl;
    cout << "input features: " << input.size() << endl;
    cout << endl;
    cout << "found: " << matched_pairs << " (" << String::number(100.0 * matched_pairs / truth.size(), 2) << "%)" << endl;
    cout << "half found : " << half_matched_pairs << " (" << String::number(100.0 * half_matched_pairs / truth.size(), 2) << "%)" << endl;
    cout << "not found : " << truth.size() - (matched_pairs + half_matched_pairs) << " (" << String::number(100.0 - 100.0 * (matched_pairs + half_matched_pairs) / truth.size(), 2) << "%)" << endl;
    cout << endl;
    cout << "relative pair ratios: " << fiveNumberQuotients(i_ratio, t_ratio, 3) << endl;
    cout << "pair distance RT : " << fiveNumbers(rt_diffs, 2) << endl;
    cout << "pair distance m/z: " << fiveNumbers(mz_diffs, 2) << endl;

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPLabeledEval tool;
  return tool.main(argc, argv);
}

/// @endcond
