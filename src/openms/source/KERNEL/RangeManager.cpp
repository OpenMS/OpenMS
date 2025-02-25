// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/RangeManager.h>

#include <ostream>

namespace OpenMS
{
  std::ostream& operator<<(std::ostream& out, const RangeBase& b)
  {
    out << "[" << b.getMin() << ", " << b.getMax() << "]";
    return out;
  }
  
  std::ostream& operator<<(std::ostream& out, const RangeRT& range)
  {
    out << "rt: " << (OpenMS::RangeBase) range << "\n";
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const RangeMZ& range)
  {
    out << "mz: " << (RangeBase) range << "\n";
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const RangeIntensity& range)
  {
    out << "intensity: " << (RangeBase) range << "\n";
    return out;
  }

  std::ostream& operator<<(std::ostream& out, const RangeMobility& range)
  {
    out << "mobility: " << (RangeBase) range << "\n";
    return out;
  }

  RangeBase::operator RangeRT() const
  {
    RangeRT out;
    out.RangeBase::operator=(*this);
    return out;
  }
  RangeBase::operator RangeMZ() const
  {
    RangeMZ out;
    out.RangeBase::operator=(*this);
    return out;
  }
  RangeBase::operator RangeIntensity() const
  {
    RangeIntensity out;
    out.RangeBase::operator=(*this);
    return out;
  }
  RangeBase::operator RangeMobility() const
  {
    RangeMobility out;
    out.RangeBase::operator=(*this);
    return out;
  }
} // namespace OpenMS
