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
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <limits>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace std;

namespace OpenMS
{
  IsotopeDistribution::IsotopeDistribution() :
    max_isotope_(0)
  {
    distribution_.push_back(make_pair<Size, double>(0, 1));
  }

  IsotopeDistribution::IsotopeDistribution(Size max_isotope) :
    max_isotope_(max_isotope)
  {
    distribution_.push_back(make_pair<Size, double>(0, 1));
  }

  IsotopeDistribution::IsotopeDistribution(const IsotopeDistribution & isotope_distribution) :
    max_isotope_(isotope_distribution.max_isotope_),
    distribution_(isotope_distribution.distribution_)
  {
  }

  IsotopeDistribution::~IsotopeDistribution()
  {
  }

  void IsotopeDistribution::setMaxIsotope(Size max_isotope)
  {
    max_isotope_ = max_isotope;
  }

  Size IsotopeDistribution::getMaxIsotope() const
  {
    return max_isotope_;
  }

  IsotopeDistribution & IsotopeDistribution::operator=(const IsotopeDistribution & iso)
  {
    if (this != &iso)
    {
      distribution_ = iso.distribution_;
      max_isotope_ = iso.max_isotope_;
    }
    return *this;
  }

  IsotopeDistribution IsotopeDistribution::operator+(const IsotopeDistribution & iso) const
  {
    ContainerType result;
    convolve_(result, distribution_, iso.distribution_);
    IsotopeDistribution result_iso;
    result_iso.setMaxIsotope(max_isotope_);
    result_iso.set(result);
    return result_iso;
  }

  IsotopeDistribution & IsotopeDistribution::operator+=(const IsotopeDistribution & iso)
  {
    ContainerType result;
    convolve_(result, distribution_, iso.distribution_);
    distribution_ = result;
    return *this;
  }

  IsotopeDistribution & IsotopeDistribution::operator*=(Size factor)
  {
    ContainerType result;
    convolvePow_(result, distribution_, factor);
    distribution_ = result;
    return *this;
  }

  IsotopeDistribution IsotopeDistribution::operator*(Size factor) const
  {
    ContainerType result;
    convolvePow_(result, distribution_, factor);
    IsotopeDistribution result_iso;
    result_iso.setMaxIsotope(max_isotope_);
    result_iso.set(result);
    return result_iso;
  }

  void IsotopeDistribution::set(const ContainerType & distribution)
  {
    distribution_ = distribution;
  }

  const IsotopeDistribution::ContainerType & IsotopeDistribution::getContainer() const
  {
    return distribution_;
  }

  Size IsotopeDistribution::getMax() const
  {
    if (distribution_.empty())
    {
      return 0;
    }
    return distribution_[distribution_.size() - 1].first;
  }

  Size IsotopeDistribution::getMin() const
  {
    if (distribution_.empty())
    {
      return 0;
    }
    return distribution_[0].first;
  }

  Size IsotopeDistribution::size() const
  {
    return distribution_.size();
  }

  void IsotopeDistribution::clear()
  {
    distribution_.clear();
    max_isotope_ = 0;
  }

  void IsotopeDistribution::estimateFromPeptideWeight(double average_weight)
  {
    // Element counts are from Senko's Averagine model
    estimateFromWeightAndComp(average_weight, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
  }

  void IsotopeDistribution::estimateFromPeptideWeightAndS(double average_weight, UInt S)
  {
    // Element counts are from Senko's Averagine model, excluding sulfur.
    estimateFromWeightAndCompAndS(average_weight, S, 4.9384, 7.7583, 1.3577, 1.4773, 0);
  }

  void IsotopeDistribution::estimateFromRNAWeight(double average_weight)
  {
    estimateFromWeightAndComp(average_weight, 9.75, 12.25, 3.75, 7, 0, 1);
  }

  void IsotopeDistribution::estimateFromDNAWeight(double average_weight)
  {
    estimateFromWeightAndComp(average_weight, 9.75, 12.25, 3.75, 6, 0, 1);
  }

  void IsotopeDistribution::estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P)
  {
      EmpiricalFormula ef;
      ef.estimateFromWeightAndComp(average_weight, C, H, N, O, S, P);
      distribution_ = ef.getIsotopeDistribution(max_isotope_).getContainer();
  }

  void IsotopeDistribution::estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P)
  {
    EmpiricalFormula ef;
    ef.estimateFromWeightAndCompAndS(average_weight, S, C, H, N, O, P);
    distribution_ = ef.getIsotopeDistribution(max_isotope_).getContainer();
  }

  void IsotopeDistribution::estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes)
  {
    // Element counts are from Senko's Averagine model
    estimateForFragmentFromWeightAndComp(average_weight_precursor, average_weight_fragment, precursor_isotopes, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
  }

  void IsotopeDistribution::estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor, UInt S_precursor, double average_weight_fragment, UInt S_fragment, const std::set<UInt>& precursor_isotopes)
  {
    UInt max_depth = *std::max_element(precursor_isotopes.begin(), precursor_isotopes.end())+1;

    double average_weight_comp_fragment = average_weight_precursor - average_weight_fragment;
    double S_comp_fragment = S_precursor - S_fragment;

    IsotopeDistribution id_comp_fragment(max_depth), id_fragment(max_depth);

    id_fragment.estimateFromPeptideWeightAndS(average_weight_fragment, S_fragment);
    id_comp_fragment.estimateFromPeptideWeightAndS(average_weight_comp_fragment, S_comp_fragment);

    calcFragmentIsotopeDist(id_fragment, id_comp_fragment, precursor_isotopes);
  }

  void IsotopeDistribution::estimateForFragmentFromRNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes)
  {
    estimateForFragmentFromWeightAndComp(average_weight_precursor, average_weight_fragment, precursor_isotopes, 9.75, 12.25, 3.75, 7, 0, 1);
  }

  void IsotopeDistribution::estimateForFragmentFromDNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes)
  {
    estimateForFragmentFromWeightAndComp(average_weight_precursor, average_weight_fragment, precursor_isotopes, 9.75, 12.25, 3.75, 6, 0, 1);
  }

  void IsotopeDistribution::estimateForFragmentFromWeightAndComp(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes, double C, double H, double N, double O, double S, double P)
  {
    UInt max_depth = *std::max_element(precursor_isotopes.begin(), precursor_isotopes.end()) + 1;

    EmpiricalFormula ef_fragment;
    ef_fragment.estimateFromWeightAndComp(average_weight_fragment, C, H, N, O, S, P);
    IsotopeDistribution id_fragment = ef_fragment.getIsotopeDistribution(max_depth);

    EmpiricalFormula ef_comp_frag;
    ef_comp_frag.estimateFromWeightAndComp(average_weight_precursor-average_weight_fragment, C, H, N, O, S, P);
    IsotopeDistribution id_comp_fragment = ef_comp_frag.getIsotopeDistribution(max_depth);

    calcFragmentIsotopeDist(id_fragment, id_comp_fragment, precursor_isotopes);
  }

  void IsotopeDistribution::calcFragmentIsotopeDist(const IsotopeDistribution& fragment_isotope_dist, const IsotopeDistribution& comp_fragment_isotope_dist, const std::set<UInt>& precursor_isotopes)
  {
    ContainerType result;
    calcFragmentIsotopeDist_(result, fragment_isotope_dist.distribution_, comp_fragment_isotope_dist.distribution_, precursor_isotopes);
    distribution_ = result;
  }

  bool IsotopeDistribution::operator==(const IsotopeDistribution & isotope_distribution) const
  {
    return max_isotope_ == isotope_distribution.max_isotope_ &&
           distribution_ == isotope_distribution.distribution_;
  }

  bool IsotopeDistribution::operator!=(const IsotopeDistribution & isotope_distribution) const
  {
    return !(isotope_distribution == *this);
  }

  IsotopeDistribution::ContainerType IsotopeDistribution::fillGaps_(const IsotopeDistribution::ContainerType& id) const
  {
    ContainerType id_gapless;
    Size mass = id.begin()->first;
    for (ContainerType::const_iterator it = id.begin(); it < id.end(); ++mass) // go through all masses
    {
      if (it->first != mass) 
      { // missing an entry
        id_gapless.push_back(make_pair(mass, 0.0));
      }
      else 
      { // mass is registered already
        id_gapless.push_back(*it); // copy
        ++it;  // ... and advance
      }
    }
    return id_gapless;
  }

  void IsotopeDistribution::convolve_(ContainerType & result, const ContainerType & left, const ContainerType & right) const
  {
    if (left.empty() || right.empty())
    {
      result.clear();
      return;
    }

    
    // ensure the isotope cluster has no gaps 
    // (e.g. from Bromine there is only Bromine-79 & Bromine-81, so we need to insert Bromine-80 with zero probability)
    ContainerType left_l = fillGaps_(left);
    ContainerType right_l = fillGaps_(right);

    ContainerType::size_type r_max = left_l.size() + right_l.size() - 1;

    if ((ContainerType::size_type)max_isotope_ != 0 && r_max > (ContainerType::size_type)max_isotope_)
    {
      r_max = (ContainerType::size_type)max_isotope_;
    }

    // pre-fill result with masses
    result.resize(r_max);
    for (ContainerType::size_type i = 0; i != r_max; ++i)
    {
      result[i] = make_pair(left_l[0].first + right_l[0].first + i, 0);
    }

    // fill result with probabilities
    // (we loop backwards because then the small products tend to come first, for better numerics)
    for (SignedSize i = left_l.size() - 1; i >= 0; --i)
    {
      for (SignedSize j = min<SignedSize>(r_max - i, right_l.size()) - 1; j >= 0; --j)
      {
        result[i + j].second += left_l[i].second * right_l[j].second;
      }
    }
  }

  void IsotopeDistribution::convolvePow_(ContainerType & result, const ContainerType & input, Size n) const
  {
    // TODO: Maybe use FFT convolve
    if (n == 1)
    {
      result = input;
      return;
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

    ContainerType input_l = fillGaps_(input);

    // get started
    if (n & 1)
    {
      result = input_l;
    }
    else
    {
      result.clear();
      result.push_back(make_pair<Size, double>(0, 1.0));
    }

    ContainerType intermediate;

    // to avoid taking unnecessary squares, we check the loop condition
    // somewhere in the middle
    ContainerType convolution_power;
    convolveSquare_(convolution_power, input_l);
    for (Size i = 1;; ++i)
    {
      if (n & (Size(1) << i))
      {
        convolve_(intermediate, result, convolution_power);
        swap(intermediate, result);
      }
      // check the loop condition
      if (i >= log2n)
        break;

      // prepare next round
      convolveSquare_(intermediate, convolution_power);
      swap(intermediate, convolution_power);
    }
  }

  void IsotopeDistribution::convolveSquare_(ContainerType & result, const ContainerType & input) const
  {
    result.clear();
    ContainerType::size_type r_max = 2 * input.size() - 1;

    if ((ContainerType::size_type)max_isotope_ != 0 && (ContainerType::size_type)(max_isotope_ + 1) < r_max)
    {
      r_max = (ContainerType::size_type)(max_isotope_ + 1);
    }

    result.resize(r_max);
    for (ContainerType::size_type i = 0; i != r_max; ++i)
    {
      result[i] = make_pair(2 * input[0].first + i, 0);
    }

    // we loop backwards because then the small products tend to come first
    // (for better numerics)
    for (SignedSize i = input.size() - 1; i >= 0; --i)
    {
      for (SignedSize j = min<SignedSize>(r_max - i, input.size()) - 1; j >= 0; --j)
      {
        result[i + j].second += input[i].second * input[j].second;
      }
    }

    return;
  }

  void IsotopeDistribution::calcFragmentIsotopeDist_(ContainerType& result, const ContainerType& fragment_isotope_dist, const ContainerType& comp_fragment_isotope_dist, const std::set<UInt>& precursor_isotopes)
  {
    if (fragment_isotope_dist.empty() || comp_fragment_isotope_dist.empty())
    {
      result.clear();
      return;
    }

    // ensure the isotope cluster has no gaps
    // (e.g. from Bromine there is only Bromine-79 & Bromine-81, so we need to insert Bromine-80 with zero probability)
    ContainerType fragment_isotope_dist_l = fillGaps_(fragment_isotope_dist);
    ContainerType comp_fragment_isotope_dist_l = fillGaps_(comp_fragment_isotope_dist);

    ContainerType::size_type r_max = fragment_isotope_dist_l.size();

    if ((ContainerType::size_type)max_isotope_ != 0 && r_max > (ContainerType::size_type)max_isotope_)
    {
      r_max = (ContainerType::size_type)max_isotope_;
    }

    // pre-fill result with masses
    result.resize(r_max);
    for (ContainerType::size_type i = 0; i != r_max; ++i)
    {
      result[i] = make_pair(fragment_isotope_dist_l[0].first + i, 0);
    }

    // Example: Let the Precursor formula be C2, and assume precursors M0, M+1, and M+2 were isolated.
    // Let the fragment formula be C1, and therefore the complementary fragment formula is also C1
    //
    // let fi = fragment formula's isotope, pi = precursor formula's isotope, ci = complementary fragment formula's isotope
    // let P(fi=x) be the probability of the formula existing as isotope x in precursor form (i.e. random sample from the universe)
    //
    // We want to calculate the probability the fragment will be isotope x given that we isolated precursors M0,M+1,M+2
    //
    // P(fi=0|pi=0 or pi=1 or pi=2) = P(fi=0) * P(pi=0 or pi=1 or pi=2|fi=0) / P(pi=0 or pi=1 or pi=2)  // bayes' theorem
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
    for (Size i = 0; i < fragment_isotope_dist_l.size(); ++i)
    {
      for (std::set<UInt>::const_iterator precursor_itr = precursor_isotopes.begin(); precursor_itr != precursor_isotopes.end(); ++precursor_itr)
      {
        if (*precursor_itr >= i &&
                (*precursor_itr-i) < comp_fragment_isotope_dist_l.size())
        {
          result[i].second += comp_fragment_isotope_dist_l[*precursor_itr-i].second;
        }
      }
      result[i].second *= fragment_isotope_dist_l[i].second;
    }
  }

  void IsotopeDistribution::renormalize()
  {
    if (distribution_.size() != 0)
    {
      double sum(0);
      // loop backwards as most distributions contains a lot of small values at the end
      for (ContainerType::const_reverse_iterator it = distribution_.rbegin(); it != distribution_.rend(); ++it)
      {
        sum += it->second;
      }

      for (Iterator it = distribution_.begin(); it != distribution_.end(); ++it)
      {
        it->second /= sum;
      }
    }
    return;
  }

  void IsotopeDistribution::trimRight(double cutoff)
  {
    ContainerType::reverse_iterator riter = distribution_.rbegin();

    // loop from right to left until an entry is larger than the cutoff
    for (; riter != distribution_.rend(); ++riter)
    {
      if (riter->second >= cutoff)
        break;
    }
    // trim the container
    distribution_.resize(riter.base() - distribution_.begin());
  }

  void IsotopeDistribution::trimLeft(double cutoff)
  {
    for (ContainerType::iterator iter = distribution_.begin(); iter != distribution_.end(); ++iter)
    {
      if (iter->second >= cutoff)
      {
        distribution_.erase(distribution_.begin(), iter);
        break;
      }
    }
  }

}
