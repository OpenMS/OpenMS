// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>
#include <ostream>

namespace OpenMS::ims
{

  /**
    @note Value for electron mass is taken from
    @link www.mcelwee.net/html/table_of_physical_constants.html
  */
  const IMSElement::mass_type IMSElement::ELECTRON_MASS_IN_U = 0.00054858;

  IMSElement & IMSElement::operator=(const IMSElement & element)
  {
    // if one doesn't assign object to itself,
    // assign all object elements to the elements of the given object
    if (this != &element)
    {
      name_ = element.name_;
      sequence_ = element.sequence_;
      isotopes_ = element.isotopes_;
    }
    return *this;
  }

  bool IMSElement::operator==(const IMSElement & element) const
  {
    return this == &element ||
           (name_ == element.name_ &&
            sequence_ == element.sequence_ &&
            isotopes_ == element.isotopes_);
  }

  bool IMSElement::operator!=(const IMSElement & element) const
  {
    return !this->operator==(element);
  }

  std::ostream & operator<<(std::ostream & os, const IMSElement & element)
  {
    os << "name:\t" << element.getName() << "\nsequence:\t" << element.getSequence()
    << "\nisotope distribution:\n" << element.getIsotopeDistribution() << '\n';
    return os;
  }
} // namespace OpenMS  // namespace ims
