// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/Isotope.h>

#include <OpenMS/CHEMISTRY/ElementDB.h>

#include <ostream>

using namespace std;

namespace OpenMS
{
  Isotope::Isotope(const std::string & name,
            const std::string & symbol,
            unsigned int atomic_number,
            unsigned int neutrons,
            double mono_weight,
            double abundance,
            double half_life,
            Isotope::DecayMode dm) :
    Element(name, symbol, atomic_number, mono_weight, mono_weight),
    neutrons_(neutrons),
    abundance_(abundance),
    half_life_(half_life),
    decay_mode_(dm)
  {
    IsotopeDistribution iso_isotopes;
    IsotopeDistribution::ContainerType iso_container;
    iso_container.push_back(Peak1D(mono_weight, 1.0));
    iso_isotopes.set(iso_container);
    setIsotopeDistribution(iso_isotopes);
  }

  Isotope::~Isotope()
  {
  }

  /// Get corresponding element
  const Element* Isotope::getElement() const
  {
    return ElementDB::getInstance()->getElement( atomic_number_ );
  }

  void Isotope::setHalfLife(double hl)
  {
    half_life_ = hl;
  }

  double Isotope::getHalfLife() const
  {
    return half_life_;
  }
  void Isotope::setAbundance(double ab)
  {
    abundance_ = ab;
  }
  double Isotope::getAbundance() const
  {
    return abundance_;
  }
  void Isotope::setNeutrons(int n)
  {
    neutrons_ = n;
  }
  int Isotope::getNeutrons() const
  {
    return neutrons_;
  }
  Isotope::DecayMode Isotope::getDecayMode() const
  {
    return decay_mode_;
  }

  void Isotope::setDecayMode(DecayMode dm)
  {
    decay_mode_ = dm;
  }

  std::ostream & operator<<(std::ostream & os, const Isotope & isotope)
  {
    os  << isotope.name_ << " "
    << isotope.symbol_ << " Z="
    << isotope.atomic_number_ << " N="
    << isotope.neutrons_ << " : "
    << isotope.mono_weight_ << " "
    << isotope.abundance_ * 100 << "% : ";
    if (isotope.isStable())
    {
      os << "stable";
    }
    else
    {
      os << "half life: " << isotope.half_life_ << " s "
      << isotope.decay_mode_;
    }

    return os;
  }
} // namespace OpenMS

