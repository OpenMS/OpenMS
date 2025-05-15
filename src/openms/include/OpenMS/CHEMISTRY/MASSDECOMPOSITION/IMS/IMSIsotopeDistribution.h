// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <algorithm> // std::min
#include <vector>
#include <iosfwd>

#include <OpenMS/config.h>

namespace OpenMS
{

  namespace ims
  {

    /**
      @brief Represents a distribution of isotopes restricted to the first K elements.

      Represents a distribution of isotopes of chemical elements as a list
      of peaks each as a pair of mass and abundance. @c IsotopeDistribution
      unlike @c IsotopeSpecies has one abundance per a nominal mass.
      Here is an example in the format (mass; abundance %)
      for molecule H2O (values are taken randomly):

      - IsotopeDistribution
        (18.00221; 99.03 %)
        (19.00334; 0.8 %)
        (20.00476; 0.17 %)

      - IsotopeSpecies
        (18.00197; 98.012 %)
        (18.00989; 1.018 %)
        (19.00312; 0.683 %)
        (19.00531; 0.117 %)
        (20.00413; 0.134 %)
        (20.00831; 0.036 %)

      To the sake of faster computations distribution is restricted
      to the first K elements, where K can be set by adjusting size
      @c SIZE of distribution. @note For the elements most abundant in
      living beings (CHNOPS) this restriction is negligible, since abundances
      decrease dramatically in isotopes order and are usually of no interest
      starting from +10 isotope.

      @c IsotopeDistribution implements folding with other distribution using an
      algorithm described in details in paper:
      Boecker et al. "Decomposing metabolic isotope patterns" WABI 2006. doi: 10.1007/11851561_2

      Folding with itself is done using Russian Multiplication Scheme.

      @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>
    */
    class OPENMS_DLLAPI IMSIsotopeDistribution
    {

public:
      /// Type of isotope mass.
      typedef double mass_type;

      /// Type of isotope abundance.
      typedef double abundance_type;

      /// Type of isotope nominal mass.
      typedef unsigned int nominal_mass_type;

      /// @brief Structure that represents an isotope peak - pair of mass and abundance.
      struct Peak
      {
        Peak(mass_type local_mass = 0.0, abundance_type local_abundance = 0.0) :
          mass(local_mass), abundance(local_abundance)
        {}

        bool operator==(const Peak & peak) const
        {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
          return peak.mass == mass && peak.abundance == abundance;
#pragma clang diagnostic pop
        }

        mass_type mass;
        abundance_type abundance;
      };

      /// Type of isotope peak.
      typedef Peak peak_type;

      /// Type of container to store peaks.
      typedef std::vector<peak_type> peaks_container;

      /// Type of iterator over container with peaks.
      typedef peaks_container::iterator peaks_iterator;

      /// Type of const iterator over container with peaks.
      typedef peaks_container::const_iterator const_peaks_iterator;

      /// Type of peaks container's size.
      typedef peaks_container::size_type size_type;

      /// Type of container with isotope masses.
      typedef std::vector<mass_type> masses_container;

      /// Type of iterator over container with isotope masses.
      typedef masses_container::iterator masses_iterator;

      /// Type of const iterator over container with isotope masses.
      typedef masses_container::const_iterator const_masses_iterator;

      /// Type of container with isotope abundances.
      typedef std::vector<abundance_type> abundances_container;

      /// Type of iterator over container with isotope abundances.
      typedef abundances_container::iterator abundances_iterator;

      /// Type of const iterator over container with isotope abundances.
      typedef abundances_container::const_iterator const_abundances_iterator;

      /// Error to be allowed for isotope distribution.
      static abundance_type ABUNDANCES_SUM_ERROR;

      /// Length of isotope distribution.
      static size_type SIZE;

      /// Constructor with nominal mass.
      explicit IMSIsotopeDistribution(nominal_mass_type nominalMass = 0) :
        nominal_mass_(nominalMass)
      {}

      /// Constructor with single isotope.
      explicit IMSIsotopeDistribution(mass_type mass) :
        nominal_mass_(0)
      {
        peaks_.push_back(peaks_container::value_type(mass, 1.0));
      }

      /// Constructor with isotopes and nominal mass.
      IMSIsotopeDistribution(const peaks_container & peaks,
                             nominal_mass_type nominalMass = 0) :
        peaks_(peaks),
        nominal_mass_(nominalMass)
      {}

      /// Copy constructor.
      IMSIsotopeDistribution(const IMSIsotopeDistribution & distribution) :
        peaks_(distribution.peaks_),
        nominal_mass_(distribution.nominal_mass_)
      {}

      /// Destructor.
      ~IMSIsotopeDistribution()
      {}

      /**
        Gets size of isotope distribution. @note Size is not smaller than
        predefined @c SIZE.

        @return Size of isotope distribution.
      */
      size_type size() const { return std::min(peaks_.size(), SIZE); }

      /**
        Assignment operator.

        @param distribution Isotope distribution to be assigned to this one.
        @return Reference to this object.
      */
      IMSIsotopeDistribution & operator=(
        const IMSIsotopeDistribution & distribution);

      /**
        Equality operator. Returns true, if a given @c distribution is equal
        to this one, false - otherwise.

        @return true, if a given distribution is equal to this distribution,
             false - otherwise
      */
      bool operator==(const IMSIsotopeDistribution & distribution) const;

      /**
        Inequality operator. Returns true, if a given @c distribution is
        unequal to this one, false - otherwise.

        @return true, if a given distribution is unequal to this
             distribution, false - otherwise
      */
      bool operator!=(const IMSIsotopeDistribution & distribution) const;

      /**
        Operator for folding this distribution with a given @c distribution.
        @note Operator is unary, so result is stored in this
        object itself.

        @param distribution Distribution to be folded with this one.
        @return Reference to this object.

        @see IsotopeDistribution& operator *=(unsigned int)
      */
      IMSIsotopeDistribution & operator*=(
        const IMSIsotopeDistribution & distribution);

      /**
        Operator for folding this distribution with itself @c pow times.
        @note Operator is unary, so result is stored in this object itself.

        @param pow Number of times this distribution is to be folded with
                                itself.
        @return Reference to this object.

        @see IsotopeDistribution& operator *=(const IsotopeDistribution&)
      */
      IMSIsotopeDistribution & operator*=(unsigned int pow);

      /**
        Gets a mass of isotope @c i.

        @param i An index of isotope.
        @return Mass of isotope @c i.
      */
      mass_type getMass(size_type i) const
      {
        return peaks_[i].mass + nominal_mass_ + i;
      }

      /**
        Gets an abundance of isotope @c i.

        @param i An index of isotope.
        @return An abundance of isotope @c i.
      */
      abundance_type getAbundance(size_type i) const
      {
        return peaks_[i].abundance;
      }

      /**
        Gets an average mass of all isotopes.

        @return An average mass of all isotopes.
      */
      mass_type getAverageMass() const;

      /**
        Gets a nominal mass of distribution.

        @return The nominal mass of the distribution.
      */
      nominal_mass_type getNominalMass() const { return nominal_mass_; }

      /**
        Sets a nominal mass for distribution.

        @param nominalMass The new nominal mass for the distribution.
      */
      void setNominalMass(nominal_mass_type nominalMass)
      {
        this->nominal_mass_ = nominalMass;
      }

      /**
        Gets masses of isotopes.

        @return Masses of isotopes.
      */
      masses_container getMasses() const;

      /**
        Gets abundances of isotopes.

        @return Abundances of isotopes.
      */
      abundances_container getAbundances() const;

      /**
        Normalizes distribution,
        i.e. scaling abundances to be summed up to 1 with an error
        @c ABUNDANCES_SUM_ERROR allowed.
      */
      void normalize();

      /**
        Returns true if the distribution has no peaks, false - otherwise.

        @return True if the distribution has no peaks, false - otherwise.
      */
      bool empty() const { return peaks_.empty(); }

private:
      /// Container for isotopes.
      peaks_container peaks_;

      /// Nominal mass of distribution.
      nominal_mass_type nominal_mass_;

      /// Sets peaks/isotopes container minimum size.
      void setMinimumSize_();
    };

    /**
      Prints isotope distribution to the stream @c os.

      @param os Output stream to which distribution is printed out.
      @param distribution Distribution to be printed out.
    */
    OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os,
                                            const IMSIsotopeDistribution & distribution);

  } // namespace ims
} // namespace OpenMS

