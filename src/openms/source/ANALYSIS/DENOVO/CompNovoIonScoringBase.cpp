// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

#include <numeric>

//#define ION_SCORING_DEBUG

using namespace std;

namespace OpenMS
{
  CompNovoIonScoringBase::IonScore::IonScore() :
    score(0),
    s_bion(0),
    s_yion(0),
    s_witness(0),
    position(0),
    s_isotope_pattern_1(0),
    is_isotope_1_mono(0),
    s_isotope_pattern_2(0)
  {
  }

  CompNovoIonScoringBase::IonScore::IonScore(const IonScore & rhs) :
    score(rhs.score),
    s_bion(rhs.s_bion),
    s_yion(rhs.s_yion),
    s_witness(rhs.s_witness),
    position(rhs.position),
    s_isotope_pattern_1(rhs.s_isotope_pattern_1),
    is_isotope_1_mono(rhs.is_isotope_1_mono),
    s_isotope_pattern_2(rhs.s_isotope_pattern_2)
  {
  }

  CompNovoIonScoringBase::IonScore::~IonScore()
  {
  }

  CompNovoIonScoringBase::IonScore & CompNovoIonScoringBase::IonScore::operator=(const IonScore & rhs)
  {
    if (this != &rhs)
    {
      score = rhs.score;
      s_bion = rhs.s_bion;
      s_yion = rhs.s_yion;
      s_witness = rhs.s_witness;
      position = rhs.position;
      s_isotope_pattern_1 = rhs.s_isotope_pattern_1;
      is_isotope_1_mono = rhs.is_isotope_1_mono;
      s_isotope_pattern_2 = rhs.s_isotope_pattern_2;
    }
    return *this;
  }

  CompNovoIonScoringBase::CompNovoIonScoringBase() :
    DefaultParamHandler("CompNovoIonScoringBase"),
    fragment_mass_tolerance_(0)
  {
    defaults_.setValue("fragment_mass_tolerance", 0.4, "fragment mass tolerance");
    defaults_.setValue("decomp_weights_precision", 0.01, "precision used to calculate the decompositions, this only affects cache usage!", ListUtils::create<String>("advanced"));
    defaults_.setValue("double_charged_iso_threshold", 0.9, "minimal isotope intensity correlation of doubly charged ions to be used to score the single scored ions", ListUtils::create<String>("advanced"));
    defaults_.setValue("double_charged_iso_threshold_single", 0.99, "Isotope scoring threshold used for doubly charged ions to infer singly charged variants", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_isotope_to_score", 3, "max isotope peak to be considered in the scoring", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_decomp_weight", 600, "maximal m/z difference used to calculate the decompositions", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_isotope", 3, "max isotope used in the theoretical spectra to score", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_mz", 2000.0, "maximal m/z value used to calculate isotope distributions", ListUtils::create<String>("advanced"));

    defaultsToParam_();
  }

  CompNovoIonScoringBase::CompNovoIonScoringBase(const CompNovoIonScoringBase & rhs) :
    DefaultParamHandler(rhs)
  {
    updateMembers_();
  }

  CompNovoIonScoringBase & CompNovoIonScoringBase::operator=(const CompNovoIonScoringBase & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
      updateMembers_();
      // TODO
    }
    return *this;
  }

  CompNovoIonScoringBase::~CompNovoIonScoringBase()
  {
  }

  void CompNovoIonScoringBase::addSingleChargedIons_(Map<double, IonScore> & ion_scores, PeakSpectrum & CID_spec)
  {
    double double_charged_iso_threshold_single((double)param_.getValue("double_charged_iso_threshold_single"));
    PeakSpectrum CID_spec_new = CID_spec;
    for (PeakSpectrum::ConstIterator it = CID_spec.begin(); it != CID_spec.end(); ++it)
    {
      if (it->getPosition()[0] < CID_spec.getPrecursors().begin()->getMZ() / 2.0)
      {
        double score = scoreIsotopes_(CID_spec, it, ion_scores, 2);
        if (score > double_charged_iso_threshold_single)
        {
          // infer this peak as single charged variant
          double mz_comp = it->getPosition()[0] * 2.0 - Constants::PROTON_MASS_U;
          bool found(false);
          for (PeakSpectrum::ConstIterator it1 = CID_spec.begin(); it1 != CID_spec.end(); ++it1)
          {
            if (fabs(mz_comp - it1->getPosition()[0]) < fragment_mass_tolerance_)
            {
              found = true;
              break;
            }
          }

          if (!found)
          {
            Peak1D p;
            p.setIntensity(it->getIntensity());
            p.setPosition(mz_comp);
            CID_spec_new.push_back(p);
          }
        }
      }
      else
      {
        break;
      }
    }

    CID_spec = CID_spec_new;
  }

  CompNovoIonScoringBase::IsotopeType CompNovoIonScoringBase::classifyIsotopes_(const PeakSpectrum & spec, PeakSpectrum::ConstIterator it)
  {
    double it_pos(it->getPosition()[0]);

    // is there a peak left of it with diff 1Da?
    for (PeakSpectrum::ConstIterator it1 = it; it1 != spec.end(); --it1)
    {
      double it1_pos(it1->getPosition()[0]);

      if (it1 == spec.begin() || fabs(it_pos - it1_pos) > 1.5)
      {
        break;
      }

      if (fabs(fabs(it_pos - it1_pos) - 1.0) < fragment_mass_tolerance_)
      {
        return CHILD;
      }
    }

    // is there a peak right of it with diff 1Da?
    for (PeakSpectrum::ConstIterator it1 = it; it1 != spec.end(); ++it1)
    {
      double it1_pos(it1->getPosition()[0]);
      if (fabs(fabs(it_pos - it1_pos) - 1.0) < fragment_mass_tolerance_)
      {
        return PARENT;
      }

      if (fabs(it_pos - it1_pos) > 1.5)
      {
        break;
      }
    }

    return LONE;
  }

  double CompNovoIonScoringBase::scoreIsotopes_(const PeakSpectrum & CID_spec, PeakSpectrum::ConstIterator it, Map<double, IonScore> & ion_scores, Size charge)
  {
    double it_pos(it->getMZ());  // ~ weight of the fragment
    UInt max_isotope_to_score(param_.getValue("max_isotope_to_score"));
    double double_charged_iso_threshold(param_.getValue("double_charged_iso_threshold"));
    double actual_pos = it_pos;

    vector<double> iso_pattern;
    vector<PeakSpectrum::ConstIterator> iso_pattern_its;
    iso_pattern.push_back(it->getIntensity());
    iso_pattern_its.push_back(it);
    // get all peaks that have the right distance right of the given peak
    for (PeakSpectrum::ConstIterator it1 = it; it1 != CID_spec.end(); ++it1)
    {
      double it1_pos(it1->getPosition()[0]);
      if (fabs(fabs(actual_pos - it1_pos) - Constants::C13C12_MASSDIFF_U / (double)charge) < fragment_mass_tolerance_)
      {
        iso_pattern.push_back(it1->getIntensity());
        actual_pos = it1_pos;
        iso_pattern_its.push_back(it1);
      }

      if (iso_pattern.size() == max_isotope_to_score)
      {
        break;
      }
    }

    if (iso_pattern.size() == 1)
    {
      return -1;
    }

    // normalize the intensity to a sum of one
    double sum(0);
    for (vector<double>::const_iterator it1 = iso_pattern.begin(); it1 != iso_pattern.end(); ++it1)
    {
      sum += *it1;
    }

    for (vector<double>::iterator it1 = iso_pattern.begin(); it1 != iso_pattern.end(); ++it1)
    {
      *it1 /= sum;
    }

    // get the theoretical isotope distribution
    CoarseIsotopePatternGenerator solver(iso_pattern.size());
    auto iso_dist = solver.estimateFromPeptideWeight((it_pos - charge * Constants::PROTON_MASS_U) * charge + Constants::PROTON_MASS_U);

    // compare the distribution sizes
    if (iso_dist.size() != iso_pattern.size())
    {
      cerr << "scoreIsotopes: error istope distributions have differing sizes" << endl;
      return -1;
    }

    // calculate simple correlation score
    double score(0.0);

    double numerator(0), auto1(0), auto2(0);
    for (Size i = 0; i != iso_dist.size(); ++i)
    {
      numerator += iso_dist.getContainer()[i].getIntensity() * iso_pattern[i];
      auto1 += iso_dist.getContainer()[i].getIntensity() * iso_dist.getContainer()[i].getIntensity();
      auto2 += iso_pattern[i] * iso_pattern[i];
    }

    score = numerator * numerator / auto1 / auto2;

    // if score is great enough, we accept it
    if (score > double_charged_iso_threshold)
    {
      if (ion_scores[it_pos].is_isotope_1_mono == 0)
      {
        ion_scores[it_pos].is_isotope_1_mono = 1;
      }

      for (Size i = 1; i < iso_pattern_its.size(); ++i)
      {
        ion_scores[iso_pattern_its[i]->getPosition()[0]].is_isotope_1_mono = -1;
#ifdef ION_SCORING_DEBUG
        cerr << "scoreIsotopes: disabling " << iso_pattern_its[i]->getPosition()[0] << endl;
#endif
      }
    }
#ifdef ION_SCORING_DEBUG
    cerr << "IsotopeScore: " << it_pos << " " << score << " " << iso_dist.size() << " " << ion_scores[it->getPosition()[0]].is_isotope_1_mono << " z=" << charge << endl;
#endif
    return score;
  }

  double CompNovoIonScoringBase::scoreIsotopes(const PeakSpectrum & spec, PeakSpectrum::ConstIterator it, Size charge)
  {
#ifdef ION_SCORING_DEBUG
    cerr << "scoreIsotopes: " << spec.size() << " " << it->getPosition()[0] << " " << it->getIntensity() << " " << charge << endl;
#endif
    double it_pos(it->getMZ()); // ~ weight of the fragment
    double actual_pos = it_pos;
    UInt max_isotope_to_score = (UInt)param_.getValue("max_isotope_to_score");

    vector<double> iso_pattern;
    iso_pattern.push_back(it->getIntensity());

    // get all peaks that have the right distance right of the given peak
    //cerr << "Scoring peaks...";
    for (PeakSpectrum::ConstIterator it1 = it; it1 != spec.end(); ++it1)
    {
      double it1_pos(it1->getMZ());

      //cerr << "PRE: " << actual_pos << " " << it1_pos << " " << Constants::NEUTRON_MASS_U << " " << charge << " " << fragment_mass_tolerance_ << endl;

      if (fabs(fabs(actual_pos - it1_pos) - Constants::NEUTRON_MASS_U / (double)charge) < fragment_mass_tolerance_ / (double)charge)
      {
#ifdef ION_SCORING_DEBUG
        cerr << actual_pos << " " << it1_pos << " " << charge << " " << fragment_mass_tolerance_ << endl;
#endif
        iso_pattern.push_back(it1->getIntensity());
        actual_pos = it1_pos;
      }

      if (iso_pattern.size() == max_isotope_to_score)
      {
        break;
      }
    }
    //cerr << "ended" << endl;

    if (iso_pattern.size() == 1)
    {
      return 0;
    }

    // normalize the intensity to a sum of one
    /*
  double sum(0);
  for (vector<double>::const_iterator it1 = iso_pattern.begin(); it1 != iso_pattern.end(); ++it1)
  {
    sum += *it1;
  }

  for (vector<double>::iterator it1 = iso_pattern.begin(); it1 != iso_pattern.end(); ++it1)
  {
    *it1 /= sum;
  }*/


    // get the theoretical isotope distribution
    CoarseIsotopePatternGenerator solver(UInt(iso_pattern.size()));
    auto iso_dist = solver.estimateFromPeptideWeight(it_pos * (double)charge - (double)(charge - 1) * Constants::PROTON_MASS_U);

    // compare the distribution sizes
    if (iso_dist.size() != iso_pattern.size())
    {
      cerr << "scoreIsotopes: error istope distributions have differing sizes" << endl;
      return -1;
    }

    // calculate simple correlation score
    double score(0.0);

    double numerator(0), auto1(0), auto2(0);
    for (Size i = 0; i != iso_dist.size(); ++i)
    {
      numerator += iso_dist.getContainer()[i].getIntensity() * iso_pattern[i];
      auto1 += iso_dist.getContainer()[i].getIntensity() * iso_dist.getContainer()[i].getIntensity();
      auto2 += iso_pattern[i] * iso_pattern[i];
    }

    //score *= accumulate(iso_pattern.begin(), iso_pattern.end(), 0.0);

    score = numerator * numerator / auto1 / auto2;
    score *= accumulate(iso_pattern.begin(), iso_pattern.end(), 0.0);
#ifdef ION_SCORING_DEBUG
    cerr << "IsotopeScore: " << it_pos << " " << score << " " << iso_dist.size() << " z=" << charge << endl;
#endif
    return score;
  }

  void CompNovoIonScoringBase::initIsotopeDistributions_()
  {
    double max_mz(param_.getValue("max_mz"));
    Size max_isotope(param_.getValue("max_isotope"));
    CoarseIsotopePatternGenerator solver(max_isotope);
    for (Size i = 1; i <= max_mz; ++i)
    {
      auto iso_dist = solver.estimateFromPeptideWeight((double)i);
      iso_dist.renormalize();
      vector<double> iso(max_isotope, 0.0);

      for (Size j = 0; j != iso_dist.size(); ++j)
      {
        iso[j] = iso_dist.getContainer()[j].getIntensity();
      }
      isotope_distributions_[i] = iso;
    }
  }

  void CompNovoIonScoringBase::updateMembers_()
  {
    fragment_mass_tolerance_ = (double)param_.getValue("fragment_mass_tolerance");

    initIsotopeDistributions_();

    return;
  }

}
