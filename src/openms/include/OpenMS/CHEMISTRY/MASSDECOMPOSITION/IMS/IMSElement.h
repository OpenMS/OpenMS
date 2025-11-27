// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <string>
#include <iosfwd>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>

namespace OpenMS
{
  namespace ims
  {

    /**
      @brief Represents a chemical atom with name and isotope distribution.

      Simulates a chemical atom with name and isotope distribution and can be
      used as a base class for more complex structures that simulate non-trivial
      bio-chemical molecules. @c Element 's name represents the atom's symbol
      in a periodical table. Sequence is by default equal to name and
      introduced for more complex molecules.

      @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>
     */
    class OPENMS_DLLAPI IMSElement
    {
public:
      /// Type of element's name.
      typedef std::string name_type;

      /// Type of element's isotope distribution.
      typedef IMSIsotopeDistribution isotopes_type;

      /// Type of isotope mass.
      typedef isotopes_type::mass_type mass_type;

      /// Type of distribution nominal mass.
      typedef isotopes_type::nominal_mass_type nominal_mass_type;

      /// Type of isotopes size.
      typedef isotopes_type::size_type size_type;

      /// Mass of electron.
      static const mass_type ELECTRON_MASS_IN_U;

      /// Empty constructor.
      IMSElement()
      {}

      /// Copy constructor.
      IMSElement(const IMSElement & element) :
        name_(element.name_),
        sequence_(element.sequence_),
        isotopes_(element.isotopes_)
      {}

      /// Constructor with name and isotope distribution.
      IMSElement(const name_type & name,
                 const isotopes_type & isotopes) :
        name_(name),
        sequence_(name),
        isotopes_(isotopes)
      {}

      /// Constructor with name and mass of single isotope.
      IMSElement(const name_type & name,
                 mass_type mass) :
        name_(name),
        sequence_(name),
        isotopes_(mass)
      {}

      /// Constructor with name and nominal mass.
      IMSElement(const name_type & name,
                 nominal_mass_type nominal_mass = 0) :
        name_(name),
        sequence_(name),
        isotopes_(nominal_mass)
      {}

      /**
        Gets element's name. @note Name represents
        a symbol of element/atom in a periodical table.

        @return Name of element.
      */
      const name_type & getName() const
      {
        return name_;
      }

      /**
        Sets element's name. @note Name represents
        a symbol of element/atom in a periodical table.

        @param name A new name to be set for element.
      */
      void setName(const name_type & name)
      {
        this->name_ = name;
      }

      /**
        Gets element's sequence.

        @return Sequence of element.
      */
      const name_type & getSequence() const
      {
        return sequence_;
      }

      /**
        Sets element's sequence.

        @param sequence A new sequence to be set for element.
      */
      void setSequence(const name_type & sequence)
      {
        this->sequence_ = sequence;
      }

      /**
        Gets element's nominal mass.

        @return A nominal mass of element.
      */
      nominal_mass_type getNominalMass() const
      {
        return isotopes_.getNominalMass();
      }

      /**
        Gets mass of element's isotope @c index.

        @param index Index of element's isotope.
        @return mass of element's isotope with a given index.
      */
      mass_type getMass(size_type index = 0) const
      {
        return isotopes_.getMass(index);
      }

      /**
        Gets element's average mass.

        @return An average mass of element.
      */
      mass_type getAverageMass() const
      {
        return isotopes_.getAverageMass();
      }

      /**
        Gets ion mass of element. By default ion lacks 1 electron,
        but this can be changed by setting other @c electrons_number.

        @param electrons_number Number of electrons lacking in ion.
      */
      mass_type getIonMass(int electrons_number = 1) const
      {
        return this->getMass() - electrons_number * ELECTRON_MASS_IN_U;
      }

      /**
        Gets element's isotope distribution.

        @return Element's isotope distribution.
      */
      const IMSIsotopeDistribution & getIsotopeDistribution() const
      {
        return isotopes_;
      }

      /**
        Sets element's isotope distribution.

        @param isotopes A new isotope distribution to be set for element.
      */
      void setIsotopeDistribution(const IMSIsotopeDistribution & isotopes)
      {
        this->isotopes_ = isotopes;
      }

      /**
        Assignment operator.

        @param element Element to be assigned to this one.
        @return Reference to this object.
      */
      IMSElement & operator=(const IMSElement & element);

      /**
        Equality operator. Returns true, if a given @c element is equal
        to this one, false - otherwise.

        @return true, if a given element is equal to this one,
             false - otherwise
      */
      bool operator==(const IMSElement & element) const;

      /**
        Inequality operator. Returns true, if a given @c element is
        unequal to this one, false - otherwise.

        @return true, if a given element is unequal to this one,
           false - otherwise.
      */
      bool operator!=(const IMSElement & element) const;

      /// Default destructor.
      virtual ~IMSElement() {}

private:
      /// Element's name.
      name_type name_;

      /// Element's sequence.
      name_type sequence_;

      /// Element's isotope distribution.
      isotopes_type isotopes_;
    };

    /**
      Prints element to the stream @c os.

      @param os Output stream to which element is printed out.
      @param element Element to be printed out.
    */
    OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const IMSElement & element);

  } // namespace ims
} // namespace OpenMS

