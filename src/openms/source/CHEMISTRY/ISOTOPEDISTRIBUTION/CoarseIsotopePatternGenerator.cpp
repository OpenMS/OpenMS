// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <include/OpenMS/CONCEPT/Constants.h>

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <functional>
#include <numeric>

using namespace std;

namespace OpenMS
{
  CoarseIsotopePatternGenerator::CoarseIsotopePatternGenerator(const Size max_isotope, const bool round_masses) :
    IsotopePatternGenerator(),
    max_isotope_(max_isotope),
    round_masses_(round_masses)
  {
  }

  CoarseIsotopePatternGenerator::~CoarseIsotopePatternGenerator() = default;


  void CoarseIsotopePatternGenerator::setMaxIsotope(const Size& max_isotope)
  {
    max_isotope_ = max_isotope;
  }

  Size CoarseIsotopePatternGenerator::getMaxIsotope() const
  {
    return max_isotope_;
  }

  void CoarseIsotopePatternGenerator::setRoundMasses(const bool round_masses)
  {
    round_masses_ = round_masses;
  }

  bool CoarseIsotopePatternGenerator::getRoundMasses() const
  {
    return round_masses_;
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::run(const EmpiricalFormula& formula) const
  {
    IsotopeDistribution result;

    auto it = formula.begin();
    for (; it != formula.end(); ++it)
    {
      IsotopeDistribution tmp = it->first->getIsotopeDistribution();
      result.set(convolve(result.getContainer(),
                          convolvePow_(tmp.getContainer(), it->second)));
    }

    // replace atomic numbers with masses.
    result.set(correctMass_(result.getContainer(), formula.getLightestIsotopeWeight()));

    result.renormalize();

    return result;
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromPeptideWeight(double average_weight)
  {
    // Element counts are from Senko's Averagine model
    return estimateFromWeightAndComp(average_weight, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromPeptideMonoWeight(double mono_weight)
  {
    // Element counts are from Senko's Averagine model
    return estimateFromMonoWeightAndComp(mono_weight, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
  }



  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromPeptideWeightAndS(double average_weight, UInt S)
  {
    // Element counts are from Senko's Averagine model, excluding sulfur.
    return estimateFromWeightAndCompAndS(average_weight, S, 4.9384, 7.7583, 1.3577, 1.4773, 0);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromRNAWeight(double average_weight)
  {
    return estimateFromWeightAndComp(average_weight, 9.75, 12.25, 3.75, 7, 0, 1);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromRNAMonoWeight(double mono_weight)
  {
    return estimateFromMonoWeightAndComp(mono_weight, 9.75, 12.25, 3.75, 7, 0, 1);
  }


  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromDNAWeight(double average_weight)
  {
    return estimateFromWeightAndComp(average_weight, 9.75, 12.25, 3.75, 6, 0, 1);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P)
  {
    EmpiricalFormula ef;
    ef.estimateFromWeightAndComp(average_weight, C, H, N, O, S, P);
    return ef.getIsotopeDistribution(*this);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromMonoWeightAndComp(double mono_weight, double C, double H, double N, double O, double S, double P)
  {
    EmpiricalFormula ef;
    ef.estimateFromMonoWeightAndComp(mono_weight, C, H, N, O, S, P);
    return ef.getIsotopeDistribution(*this);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P)
  {
    EmpiricalFormula ef;
    ef.estimateFromWeightAndCompAndS(average_weight, S, C, H, N, O, P);
    return ef.getIsotopeDistribution(*this);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes)
  {
    // Element counts are from Senko's Averagine model
    return estimateForFragmentFromWeightAndComp(average_weight_precursor, average_weight_fragment, precursor_isotopes, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor, UInt S_precursor, double average_weight_fragment, UInt S_fragment, const std::set<UInt>& precursor_isotopes) const
  {
    UInt max_depth = *std::max_element(precursor_isotopes.begin(), precursor_isotopes.end())+1;

    double average_weight_comp_fragment = average_weight_precursor - average_weight_fragment;
    double S_comp_fragment = S_precursor - S_fragment;

    // We need the solver to return atomic numbers to be compatible with calcFragmentIsotopeDist
    CoarseIsotopePatternGenerator solver(max_depth, false);

    EmpiricalFormula ef_fragment;
    ef_fragment.estimateFromWeightAndCompAndS(average_weight_fragment, S_fragment, 4.9384, 7.7583, 1.3577, 1.4773, 0);

    IsotopeDistribution id_fragment(ef_fragment.getIsotopeDistribution(solver));
    IsotopeDistribution id_comp_fragment(solver.estimateFromPeptideWeightAndS(average_weight_comp_fragment, S_comp_fragment));

    IsotopeDistribution result = calcFragmentIsotopeDist(id_fragment, id_comp_fragment, precursor_isotopes, ef_fragment.getLightestIsotopeWeight());

    return result;
  }

  // static
  IsotopeDistribution CoarseIsotopePatternGenerator::approximateFromPeptideWeight(double mass, UInt num_peaks, UInt charge)
  {
    IsotopeDistribution result;
    result.resize(num_peaks);

    // lambda * mass. Lambda is the parameter for Poisson distribution. Value (1/1800) taken from Bellew et al
    double factor = mass / 1800.0;

    // for k=0, non-normalized value is always 1
    result[0] = Peak1D(mass, 1.0f);

    float curr_intensity = 1.0f;
    for (UInt k = 1; k < num_peaks; ++k) // result[0] is always 1 anyway
    {
      curr_intensity *= factor / k; // represents (m * lambda)^k / k!

      // at some point, curr_intensity will become too small for float (which is the intensity type of Peak1D)
      result[k] = Peak1D(mass + (k * OpenMS::Constants::NEUTRON_MASS_U / charge),
                         curr_intensity != curr_intensity ? 0.0f : curr_intensity);
    }

    result.renormalize();

    return result;
  }

  std::vector<double> CoarseIsotopePatternGenerator::approximateIntensities(double mass, UInt num_peaks)
  {
    std::vector<double> result(num_peaks, 1.0);

    // lambda * mass. Lambda is the parameter of Poisson distribution. Value (1/1800) taken from Bellew et al
    double factor = mass / 1800.0;
    double curr_intensity = 1.0;
    double sum = 1.0; // result[0] is always factor^0/1 = 1, which is the reason why we start the loop at 1

    for (UInt k = 1; k < num_peaks; ++k)
    {
      curr_intensity *= factor / k; // represents (m * lambda)^k / k!

      // at some point, curr_intensity will become too small for float (which is the intensity type of Peak1D)
      result[k] = curr_intensity != curr_intensity ? 0.0: curr_intensity;
      sum += result[k];
    }

    // normalize
    for (UInt k = 0; k != result.size(); ++k)
    {
      result[k] /= sum;
    }

    return result;
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateForFragmentFromRNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes)
  {
    return estimateForFragmentFromWeightAndComp(average_weight_precursor, average_weight_fragment, precursor_isotopes, 9.75, 12.25, 3.75, 7, 0, 1);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateForFragmentFromDNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes)
  {
    return estimateForFragmentFromWeightAndComp(average_weight_precursor, average_weight_fragment, precursor_isotopes, 9.75, 12.25, 3.75, 6, 0, 1);
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::estimateForFragmentFromWeightAndComp(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes, double C, double H, double N, double O, double S, double P) const
  {
    UInt max_depth = *std::max_element(precursor_isotopes.begin(), precursor_isotopes.end()) + 1;

    // We need the solver to return atomic numbers to be compatible with calcFragmentIsotopeDist
    CoarseIsotopePatternGenerator solver(max_depth, false);

    EmpiricalFormula ef_fragment;
    ef_fragment.estimateFromWeightAndComp(average_weight_fragment, C, H, N, O, S, P);
    IsotopeDistribution id_fragment = ef_fragment.getIsotopeDistribution(solver);

    EmpiricalFormula ef_comp_frag;
    ef_comp_frag.estimateFromWeightAndComp(average_weight_precursor-average_weight_fragment, C, H, N, O, S, P);
    IsotopeDistribution id_comp_fragment = ef_comp_frag.getIsotopeDistribution(solver);

    IsotopeDistribution result = calcFragmentIsotopeDist(id_fragment, id_comp_fragment, precursor_isotopes, ef_fragment.getLightestIsotopeWeight());

    return result;
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::calcFragmentIsotopeDist(const IsotopeDistribution& fragment_isotope_dist, const IsotopeDistribution& comp_fragment_isotope_dist, const std::set<UInt>& precursor_isotopes, const double fragment_lightest_iso_mass) const
  {
    IsotopeDistribution result = calcFragmentIsotopeDist_(fragment_isotope_dist.getContainer(), comp_fragment_isotope_dist.getContainer(), precursor_isotopes);

    // replace atomic numbers with masses.
    result.set(correctMass_(result.getContainer(), fragment_lightest_iso_mass));

    return result;
  }

  IsotopeDistribution::ContainerType CoarseIsotopePatternGenerator::convolve(const IsotopeDistribution::ContainerType& left, const IsotopeDistribution::ContainerType& right) const
  {
    IsotopeDistribution::ContainerType result;

    if (left.empty() || right.empty())
    {
      result.clear();
      return result;
    }

    // ensure the isotope cluster has no gaps
    // (e.g. from Bromine there is only Bromine-79 & Bromine-81, so we need to insert Bromine-80 with zero probability)
    IsotopeDistribution::ContainerType left_l = fillGaps_(left);
    IsotopeDistribution::ContainerType right_l = fillGaps_(right);

    IsotopeDistribution::ContainerType::size_type r_max = left_l.size() + right_l.size() - 1;

    if ((IsotopeDistribution::ContainerType::size_type)max_isotope_ != 0 && r_max > (IsotopeDistribution::ContainerType::size_type)max_isotope_)
    {
      r_max = (IsotopeDistribution::ContainerType::size_type)max_isotope_;
    }

    // pre-fill result with masses
    result.resize(r_max);
    for (IsotopeDistribution::ContainerType::size_type i = 0; i != r_max; ++i)
    {
      result[i] = Peak1D(left_l[0].getMZ() + right_l[0].getMZ() + i, 0);
    }

    // fill result with probabilities
    // (we loop backwards because then the small products tend to come first, for better numerics)
    for (SignedSize i = left_l.size() - 1; i >= 0; --i)
    {
      for (SignedSize j = min<SignedSize>(r_max - i, right_l.size()) - 1; j >= 0; --j)
      {
        Peak1D& peak = result[i + j];
        Peak1D::IntensityType p_intensity = peak.getIntensity();
        peak.setIntensity( p_intensity + left_l[i].getIntensity() * right_l[j].getIntensity());
      }
    }
    return result;
  }

  IsotopeDistribution::ContainerType CoarseIsotopePatternGenerator::convolvePow_(const IsotopeDistribution::ContainerType& input, Size n) const
  {
    IsotopeDistribution::ContainerType result;
    if (input.empty())
    {
      return result;
    }
    // TODO: use FFT convolve?
    if (n == 1)
    {
      result = input; // Not needed copy
      return result;
    }

    Size log2n = 0;
    // modification by Chris to prevent infinite loop when n > 2^63
    if (n > (Size(1) << (std::numeric_limits<Size>::digits - 1)))
    {
      log2n = std::numeric_limits<Size>::digits;
    }
    else
    {
      // find binary logarithm of n
      for (; (Size(1) << log2n) < n; ++log2n)
      {
      }
    }

    IsotopeDistribution::ContainerType input_l = fillGaps_(input);

    // get started
    if (n & 1)
    {
      result = input_l;
    }
    else
    {
      result.clear();
      result.push_back(IsotopeDistribution::MassAbundance(0, 1.0));
    }


    // to avoid taking unnecessary squares, we check the loop condition
    // somewhere in the middle
    IsotopeDistribution::ContainerType convolution_power = convolveSquare_(input_l);
    for (Size i = 1;; ++i)
    {
      if (n & (Size(1) << i))
      {
        result = convolve(result, convolution_power);
      }
      // check the loop condition
      if (i >= log2n)
      {
        break;
      }
      // prepare next round
      convolution_power = convolveSquare_(convolution_power);
    }
    return result;
  }

  IsotopeDistribution::ContainerType CoarseIsotopePatternGenerator::convolveSquare_(const IsotopeDistribution::ContainerType& input) const
  {
    IsotopeDistribution::ContainerType result;
    IsotopeDistribution::ContainerType::size_type r_max = 2 * input.size() - 1;

    if ((IsotopeDistribution::ContainerType::size_type)max_isotope_ != 0 && (IsotopeDistribution::ContainerType::size_type)(max_isotope_ + 1) < r_max)
    {
      r_max = (IsotopeDistribution::ContainerType::size_type)(max_isotope_ + 1);
    }

    result.resize(r_max);
    for (IsotopeDistribution::ContainerType::size_type i = 0; i != r_max; ++i)
    {
      result[i] = Peak1D(2 * input[0].getMZ() + i, 0);
    }

    // we loop backwards because then the small products tend to come first
    // (for better numerics)
    for (SignedSize i = input.size() - 1; i >= 0; --i)
    {
      for (SignedSize j = min<SignedSize>(r_max - i, input.size()) - 1; j >= 0; --j)
      {
        result[i + j].setIntensity( result[i + j].getIntensity() + input[i].getIntensity() * input[j].getIntensity());
      }
    }

    return result;
  }

  IsotopeDistribution CoarseIsotopePatternGenerator::calcFragmentIsotopeDist_(const IsotopeDistribution::ContainerType& fragment_isotope_dist, const IsotopeDistribution::ContainerType& comp_fragment_isotope_dist, const std::set<UInt>& precursor_isotopes) const
  {

    IsotopeDistribution result;

    if (fragment_isotope_dist.empty() || comp_fragment_isotope_dist.empty())
    {
      result.clear();
      return result;
    }

    IsotopeDistribution::ContainerType::size_type r_max = fragment_isotope_dist.size();

    if ((IsotopeDistribution::ContainerType::size_type)max_isotope_ != 0 && r_max > (IsotopeDistribution::ContainerType::size_type)max_isotope_)
    {
      r_max = (IsotopeDistribution::ContainerType::size_type)max_isotope_;
    }

    // pre-fill result with masses
    result.resize(r_max);
    for (IsotopeDistribution::ContainerType::size_type i = 0; i != r_max; ++i)
    {
      result[i] = Peak1D(fragment_isotope_dist[0].getMZ() + i, 0);
    }

    // Example: Let the Precursor formula be C2, and assume precursors M0, M+1, and M+2 were isolated.
    // Let the fragment formula be C1, and therefore the complementary fragment formula is also C1
    //
    // let fi = fragment formula's isotope, pi = precursor formula's isotope, ci = complementary fragment formula's isotope
    // let P(fi=x) be the probability of the formula existing as isotope x in precursor form (i.e. random sample from the universe)
    //
    // We want to calculate the probability the fragment will be isotope x given that we isolated precursors M0,M+1,M+2
    //
    // P(fi=0|pi=0 or pi=1 or pi=2) = P(fi=0) * P(pi=0 or pi=1 or pi=2|fi=0) / P(pi=0 or pi=1 or pi=2)  // Bayes' theorem
    //        = P(fi=0) * (P(pi=0|fi=0) + P(pi=1|fi=0) + P(pi=2|fi=0)) / (P(pi=0) + P(pi=1) + P(pi=2))  // mutually exclusive events
    //        = P(fi=0) * (P(ci=0) + P(ci=1) + P(ci=2)) / (P(pi=0) + P(pi=1) + P(pi=2))                 // The only way pi=x|fi=y, is if ci=x-y
    //        = P(fi=0) * (P(ci=0) + P(ci=1) + P(ci=2))                                                 // ignore normalization for now
    //          ^this is the form we're calculating in the code, which is technically P(fi=0 and (pi=0 or pi=1 or pi=2)) because we didn't normalize
    //        = 0.9893 * (0.9893 + 0.0107 + 0)
    //          Note: In this example, P(ci=2)=0 because the complementary fragment is just C and cannot exist with 2 extra neutrons
    //
    // P(fi=1|pi=0 or pi=1 or pi=2) = P(fi=1) * P(pi=0 or pi=1 or pi=2|fi=1) / P(pi=0 or pi=1 or pi=2)
    //        = P(fi=1) * (P(pi=0|fi=1) + P(pi=1|fi=1) + P(pi=2|fi=1)) / (P(pi=0) + P(pi=1) + P(pi=2))
    //        = P(fi=1) * (P(ci=-1) + P(ci=0) + P(ci=1)) / (P(pi=0) + P(pi=1) + P(pi=2))
    //          Note: P(ci<0)=0
    //        = P(fi=1) * (P(ci=0) + P(ci=1))
    //          ^this is the form we're calculating in the code
    //        = 0.0107 * (0.9893 + 0.0107)
    //
    // P(fi=2|pi=0 or pi=1 or pi=2) = P(fi=2) * P(pi=0 or pi=1 or pi=2|fi=2) / P(pi=0 or pi=1 or pi=2)
    //        = P(fi=2) * (P(pi=0|fi=2) + P(pi=1|fi=2) + P(pi=2|fi=2)) / (P(pi=0) + P(pi=1) + P(pi=2))
    //        = P(fi=2) * (P(ci=-2) + P(ci=-1) + P(ci=0)) / (P(pi=0) + P(pi=1) + P(pi=2))
    //        = P(fi=2) * P(ci=0)
    //          ^this is the form we're calculating in the code
    //        = 0 * (0.9893)
    //          Note: In this example, P(fi=2)=0 because the fragment is just C and cannot exist with 2 extra neutrons.
    //
    // normalization is needed to get true conditional probabilities if desired.
    //
    for (Size i = 0; i < fragment_isotope_dist.size(); ++i)
    {
      for (std::set<UInt>::const_iterator precursor_itr = precursor_isotopes.begin(); precursor_itr != precursor_isotopes.end(); ++precursor_itr)
      {
        if (*precursor_itr >= i &&
            (*precursor_itr-i) < comp_fragment_isotope_dist.size())
        {
          result[i].setIntensity( result[i].getIntensity() + comp_fragment_isotope_dist[*precursor_itr-i].getIntensity());
        }
      }
      result[i].setIntensity(result[i].getIntensity() * fragment_isotope_dist[i].getIntensity());
    }
    return result;

  }

  IsotopeDistribution::ContainerType CoarseIsotopePatternGenerator::fillGaps_(const IsotopeDistribution::ContainerType& id) const
  {
    if (id.empty()) return id;
    IsotopeDistribution::ContainerType id_gapless;
    Size mass = round(id.begin()->getMZ());
    Size massend = round((--id.end())->getMZ());
    id_gapless.resize(massend - mass + 1);

    Size mass_cnt = mass;
    for (auto& peak : id_gapless)
    {
      peak.setMZ(mass_cnt);
      mass_cnt++;
    }

    for (const auto& peak : id)
    {
      Size to_set = round(peak.getMZ());
      id_gapless[to_set-mass].setIntensity(peak.getIntensity());
    }

    return id_gapless;
  }

  IsotopeDistribution::ContainerType CoarseIsotopePatternGenerator::correctMass_(const IsotopeDistribution::ContainerType& input, const double lighest_iso_weight) const
  {
    IsotopeDistribution::ContainerType result(input.size());

    for (Size i = 0; i < input.size(); ++i)
    {
      // We assume that a coarse isotopic peak is mostly composed of carbon-13's
      // and therefore use the mass difference between carbon-13 and carbon-12
      // to determine the expected mass of a coarse isotopic peak.
      double mass = lighest_iso_weight + (i * Constants::C13C12_MASSDIFF_U);
      if (getRoundMasses())
      {
        mass = round(mass);
      }

      result[i] = Peak1D(mass, input[i].getIntensity() );
    }

    return result;
  }
}
