// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/Isotope.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <ostream>
#include <stdexcept>

using namespace std;

namespace OpenMS
{
  Isotope::Isotope() :
    Element(),
    mass_number_(0),
    abundance_(0.0),
    half_life_(-1.0),
    decay_mode_(DecayMode::STABLE)
  {
  }

  Isotope::Isotope(const Isotope& isotope) = default;

  Isotope::Isotope(const string& name,
                   const string& symbol,
                   unsigned int atomic_number,
                   unsigned int mass_number,
                   double atomic_mass,
                   double abundance,
                   double half_life,
                   DecayMode decay_mode) :
    Element(),
    mass_number_(mass_number),
    abundance_(abundance),
    half_life_(half_life),
    decay_mode_(decay_mode)
  {
    // Set Element properties
    setName(name);
    setSymbol(symbol);
    setAtomicNumber(atomic_number);
    setMonoWeight(atomic_mass);
    setAverageWeight(atomic_mass); // For a single isotope, average = mono

    // Create isotope distribution with single peak
    IsotopeDistribution iso_dist;
    IsotopeDistribution::ContainerType container;
    container.push_back(Peak1D(atomic_mass, 1.0));
    iso_dist.set(container);
    setIsotopeDistribution(iso_dist);
  }

  Isotope::~Isotope() = default;

  void Isotope::setMassNumber(unsigned int mass_number)
  {
    mass_number_ = mass_number;
  }

  unsigned int Isotope::getMassNumber() const
  {
    return mass_number_;
  }

  unsigned int Isotope::getNeutrons() const
  {
    return mass_number_ - getAtomicNumber();
  }

  void Isotope::setAbundance(double abundance)
  {
    if (abundance < 0.0 || abundance > 1.0)
    {
      throw std::invalid_argument("Abundance must be between 0.0 and 1.0");
    }
    abundance_ = abundance;
  }

  double Isotope::getAbundance() const
  {
    return abundance_;
  }

  void Isotope::setHalfLife(double half_life)
  {
    half_life_ = half_life;
  }

  double Isotope::getHalfLife() const
  {
    return half_life_;
  }

  void Isotope::setDecayMode(DecayMode decay_mode)
  {
    decay_mode_ = decay_mode;
  }

  Isotope::DecayMode Isotope::getDecayMode() const
  {
    return decay_mode_;
  }

  bool Isotope::isStable() const
  {
    return half_life_ < 0.0;
  }

  bool Isotope::isIsotope() const
  {
    return true;
  }

  Isotope& Isotope::operator=(const Isotope& isotope) = default;

  bool Isotope::operator==(const Isotope& isotope) const
  {
    return Element::operator==(isotope) &&
           mass_number_ == isotope.mass_number_ &&
           abundance_ == isotope.abundance_ &&
           half_life_ == isotope.half_life_ &&
           decay_mode_ == isotope.decay_mode_;
  }

  bool Isotope::operator!=(const Isotope& isotope) const
  {
    return !(*this == isotope);
  }

  bool Isotope::operator<(const Isotope& isotope) const
  {
    if (getAtomicNumber() != isotope.getAtomicNumber())
    {
      return getAtomicNumber() < isotope.getAtomicNumber();
    }
    return mass_number_ < isotope.mass_number_;
  }

  string Isotope::decayModeToString(DecayMode mode)
  {
    switch (mode)
    {
      case DecayMode::STABLE:
        return "STABLE";
      case DecayMode::UNKNOWN:
        return "UNKNOWN";
      case DecayMode::ALPHA:
        return "ALPHA";
      case DecayMode::BETA_PLUS:
        return "BETA_PLUS";
      case DecayMode::BETA_MINUS:
        return "BETA_MINUS";
      case DecayMode::PROTON:
        return "PROTON";
      default:
        return "UNKNOWN";
    }
  }

  Isotope::DecayMode Isotope::stringToDecayMode(const string& mode_str)
  {
    if (mode_str == "STABLE")
      return DecayMode::STABLE;
    if (mode_str == "ALPHA")
      return DecayMode::ALPHA;
    if (mode_str == "BETA_PLUS")
      return DecayMode::BETA_PLUS;
    if (mode_str == "BETA_MINUS")
      return DecayMode::BETA_MINUS;
    if (mode_str == "PROTON")
      return DecayMode::PROTON;
    return DecayMode::UNKNOWN;
  }

  ostream& operator<<(ostream& os, const Isotope& isotope)
  {
    os << isotope.getName() << " "
       << isotope.getSymbol() << " "
       << "Z=" << isotope.getAtomicNumber() << " "
       << "A=" << isotope.getMassNumber() << " "
       << "mass=" << isotope.getMonoWeight() << " "
       << "abundance=" << isotope.getAbundance() * 100 << "% ";

    if (isotope.isStable())
    {
      os << "(stable)";
    }
    else
    {
      os << "half-life=" << isotope.getHalfLife() << "s "
         << "decay=" << Isotope::decayModeToString(isotope.getDecayMode());
    }
    return os;
  }

} // namespace OpenMS
