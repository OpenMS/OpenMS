// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ELEMENT_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ELEMENT_H

#include <string>
#include <ostream>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>
#include <iostream>

namespace OpenMS {
namespace ims {

/**
 * @brief Represents a chemical atom with name and isotope distribution.
 * 
 * Simulates a chemical atom with name and isotope distribution and can be
 * used as a base class for more complex structures that simulate non-trivial
 * bio-chemical molecules. @c Element 's name represents the atom's symbol 
 * in a periodical table. Sequence is by default equal to name and 
 * introduced for more complex molecules.
 * 
 * @see ComposedElement
 * 
 * @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>  
 * 
 * @ingroup alphabet 
 */
class Element {
public:
  /**
   * Type of element's name.
   */
  typedef std::string name_type;

  /**
   * Type of element's isotope distribution.
   */
  typedef IMSIsotopeDistribution isotopes_type;

  /**
   * Type of isotope mass.
   */
  typedef isotopes_type::mass_type mass_type;

  /**
   * Type of distribution nominal mass.
   */
  typedef isotopes_type::nominal_mass_type nominal_mass_type;

  /**
   * Type of isotopes size.
   */
  typedef isotopes_type::size_type size_type;

  /**
   * Mass of electron.
   */
  static const mass_type ELECTRON_MASS_IN_U;

  /**
   * Empty constructor.
   */
  Element()
  {
  }

  /**
   * Copy constructor.
   */
  Element(const Element& element) :
    name(element.name),
    sequence(element.sequence),
    isotopes(element.isotopes)
  {
  }

  /**
   * Constructor with name and isotope distribution.
   */
  Element(const name_type& name,
          const isotopes_type& isotopes) :
    name(name),
    sequence(name),
    isotopes(isotopes)
  {
  }

  /**
   * Constructor with name and mass of single isotope.
   */
  Element(	const name_type& name,
          mass_type mass):
    name(name),
    sequence(name),
    isotopes(mass)
  {
  }

  /**
   * Constructor with name and nominal mass.
   */
  Element(const name_type& name,
          nominal_mass_type nominal_mass = 0)
    : name(name),
    sequence(name),
    isotopes(nominal_mass)
  {
  }

  /**
   * Gets element's name. @note Name represents
   * a symbol of element/atom in a periodical table.
   *
   * @return Name of element.
   */
  const name_type& getName() const
  {
    return name;
  }

  /**
   * Sets element's name. @note Name represents
   * a symbol of element/atom in a periodical table.
   *
   * @param name A new name to be set for element.
   */
  void setName(const name_type& name)
  {
    this->name = name;
  }

  /**
   * Gets element's sequence.
   *
   * @return Sequence of element.
   */
  const name_type& getSequence() const
  {
    return sequence;
  }

  /**
   * Sets element's sequence.
   *
   * @param sequence A new sequence to be set for element.
   */
  void setSequence(const name_type& sequence)
  {
    this->sequence = sequence;
  }

  /**
   * Gets element's nominal mass.
   *
   * @return A nominal mass of element.
   */
  nominal_mass_type getNominalMass() const
  {
    return isotopes.getNominalMass();
  }

  /**
   * Gets mass of element's isotope @c index.
   *
   * @param index Index of element's isotope.
   * @return mass of element's isotope with a given index.
   */
  mass_type getMass(size_type index = 0) const
  {
    return isotopes.getMass(index);
  }

  /**
   * Gets element's average mass.
   *
   * @return An average mass of element.
   */
  mass_type getAverageMass() const
  {
    return isotopes.getAverageMass();
  }

  /**
   * Gets ion mass of element. By default ion lacks 1 electron,
   * but this can be changed by setting other @c electrons_number.
   *
   * @param electrons_number Number of electrons lacking in ion.
   */
  mass_type getIonMass(int electrons_number = 1) const
  {
    return this->getMass() - electrons_number * ELECTRON_MASS_IN_U;
  }

  /**
   * Gets element's isotope distribution.
   *
   * @return Element's isotope distribution.
   */
  const IMSIsotopeDistribution& getIsotopeDistribution() const {
    return isotopes;
  }

  /**
   * Sets element's isotope distribution.
   *
   * @param isotopes A new isotope distribution to be set for element.
   */
  void setIsotopeDistribution(const IMSIsotopeDistribution& isotopes) {
    this->isotopes = isotopes;
  }

  /**
   * Assignment operator.
   *
   * @param element Element to be assigned to this one.
   * @return Reference to this object.
   */
  Element& operator =(const Element& element);

  /**
   * Equality operator. Returns true, if a given @c element is equal
   * to this one, false - otherwise.
   *
   * @return true, if a given element is equal to this one,
   * 		   false - otherwise
   */
  bool operator ==(const Element& element) const;

  /**
   * Inequality operator. Returns true, if a given @c element is
   * unequal to this one, false - otherwise.
   *
   * @return true, if a given element is unequal to this one,
   * 		   false - otherwise.
   */
  bool operator !=(const Element& element) const;

  /**
   * Default destructor.
   */
  virtual ~Element() {}

private:
  /**
   * Element's name.
   */
  name_type name;

  /**
   * Element's sequence.
   */
  name_type sequence;

  /**
   * Element's isotope distribution.
   */
  isotopes_type isotopes;
};

/**
 * Prints element to the stream @c os.
 *
 * @param os Output stream to which element is printed out.
 * @param element Element to be printed out.
 */
std::ostream& operator<<(std::ostream& os, const Element& element);

} // namespace ims
} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ELEMENT_H
