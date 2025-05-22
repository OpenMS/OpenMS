// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
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
    std::cerr << "DEBUG: Isotope::getElement() called for " << name_ << " with atomic number " << atomic_number_ << std::endl;
    const Element* e = ElementDB::getInstance()->getElement( atomic_number_ );
    std::cerr << "DEBUG: Isotope::getElement() returning pointer " << e << std::endl;
    return e;
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

