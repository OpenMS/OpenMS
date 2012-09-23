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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  namespace ims
  {
    Weights & Weights::operator=(const Weights & other)
    {
      if (this != &other)
      {
        alphabet_masses_ = other.alphabet_masses_;
        precision_ = other.precision_;
        weights_ = other.weights_;
      }
      return *this;
    }

    void Weights::setPrecision(Weights::alphabet_mass_type precision)
    {
      this->precision_ = precision;
      weights_.clear();
      // convert alphabet masses (double) to integer masses (weights) with the given precision
      for (alphabet_masses_type::size_type i = 0; i < alphabet_masses_.size(); ++i)
      {
        weights_.push_back(static_cast<weight_type>(floor((alphabet_masses_[i] / precision) + 0.5)));
      }
    }

    void Weights::swap(size_type index1, size_type index2)
    {
      weight_type weight = weights_[index1];
      weights_[index1] = weights_[index2];
      weights_[index2] = weight;

      alphabet_mass_type mass = alphabet_masses_[index1];
      alphabet_masses_[index1] = alphabet_masses_[index2];
      alphabet_masses_[index2] = mass;
    }

    Weights::alphabet_mass_type Weights::getParentMass(const std::vector<unsigned int> & decomposition) const
    {
      // checker whether the passed decomposition is applicable
      if (alphabet_masses_.size() != decomposition.size())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("The passed decomposition has the wrong size. Expected ") + String(alphabet_masses_.size()) + String(" but got ") + String(decomposition.size()) + String("."));
      }

      alphabet_mass_type parent_mass = 0;

      for (std::vector<unsigned int>::size_type i = 0; i < decomposition.size(); ++i)
      {
        parent_mass += alphabet_masses_[i] * decomposition[i];
      }
      return parent_mass;
    }

    bool Weights::divideByGCD()
    {
      if (weights_.size() < 2)
      {
        return false;
      }
      weight_type d = Math::gcd(weights_[0], weights_[1]);
      for (weights_type::size_type i = 2; i < weights_.size(); ++i)
      {
        d = Math::gcd(d, weights_[i]);
        if (d == 1)
        {
          return false;
        }
      }
      // if we're here: d != 1
      precision_ *= d;

      // rescales the integer weights. Don't use setPrecision() here since
      // the result could be different due to rounding errors.
      for (weights_type::size_type i = 0; i < weights_.size(); ++i)
      {
        weights_[i] /= d;
      }
      return true;
    }

    Weights::alphabet_mass_type Weights::getMinRoundingError() const
    {
      alphabet_mass_type min_error = 0;
      for (size_type i = 0; i < weights_.size(); ++i)
      {
        alphabet_mass_type error = (precision_ * static_cast<alphabet_mass_type>(weights_[i]) - alphabet_masses_[i]) / alphabet_masses_[i];
        if (error < 0 && error < min_error)
        {
          min_error = error;
        }
      }
      return min_error;
    }

    Weights::alphabet_mass_type Weights::getMaxRoundingError() const
    {
      alphabet_mass_type max_error = 0;
      for (size_type i = 0; i < weights_.size(); ++i)
      {
        alphabet_mass_type error = (precision_ * static_cast<alphabet_mass_type>(weights_[i]) - alphabet_masses_[i]) / alphabet_masses_[i];
        if (error > 0 && error > max_error)
        {
          max_error = error;
        }
      }
      return max_error;
    }

    std::ostream & operator<<(std::ostream & os, const Weights & weights)
    {
      for (Weights::size_type i = 0; i < weights.size(); ++i)
      {
        os << weights.getWeight(i) << std::endl;
      }
      return os;
    }

  } // namespace ims
} // namespace OpenMS
