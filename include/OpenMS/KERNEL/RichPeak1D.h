// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_RICHPEAK1D_H
#define OPENMS_KERNEL_RICHPEAK1D_H

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{

  /**
    @brief A 1-dimensional raw data point or peak mith meta information.

    This datastructure is intended for continuous data or peak data.
    If wou do not need to annotated single peaks with meta data, use Peak1D instead.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI RichPeak1D :
    public Peak1D,
    public MetaInfoInterface
  {
public:

    /// Default constructor
    inline RichPeak1D() :
      Peak1D(),
      MetaInfoInterface()
    {}

    /// Copy constructor
    inline RichPeak1D(const RichPeak1D & p) :
      Peak1D(p),
      MetaInfoInterface(p)
    {}

    /// Destructor
    ~RichPeak1D()
    {}

    /// Assignment operator
    inline RichPeak1D & operator=(const RichPeak1D & rhs)
    {
      if (this == &rhs) return *this;

      Peak1D::operator=(rhs);
      MetaInfoInterface::operator=(rhs);

      return *this;
    }

    /// Equality operator
    inline bool operator==(const RichPeak1D & rhs) const
    {
      return Peak1D::operator==(rhs) &&
             MetaInfoInterface::operator==(rhs);
    }

    /// Equality operator
    inline bool operator!=(const RichPeak1D & rhs) const
    {
      return !(operator==(rhs));
    }

  };

} // namespace OpenMS

#endif // OPENMS_KERNEL_RICHPEAK1D_H
