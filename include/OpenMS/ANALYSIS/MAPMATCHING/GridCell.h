// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_GRIDCELL_H
#define OPENMS_ANALYSIS_MAPMATCHING_GRIDCELL_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseMapping.h>

#include <vector>

namespace OpenMS
{

  /** @brief 2-dimensional grid cell over a map.
   */
  class GridCell : public DRange<2>
  {
  public:
    enum { DIMENSION = 2 };
    typedef DoubleReal  CoordinateType;
    typedef DPosition<2> PositionType;
    typedef BaseMapping MappingType;
    typedef std::vector<MappingType*> MappingVector;

    /** @name Constructors and Destructor
     */
    //@{
    /// Default constructor
    GridCell() :
        DRange<2>(),
        mappings_()
    { }

    /// "Convenience contructor"
    GridCell(CoordinateType i, CoordinateType j,
             CoordinateType k, CoordinateType l) :
        DRange<2>(i, j, k, l),
        mappings_()
    { }

    /// Copy constructor
    GridCell(const GridCell& gc)
        : DRange<2>(gc)
    {
      setMappings(gc.mappings_);
    }

    /// Destructor
    virtual ~GridCell()
    {}
    //@}

    /// assignment operator
    GridCell& operator = (const GridCell& rhs);

    /** @name Accesssor methods
    */
    //@{
    /// Set transform
    inline void setMappings(const MappingVector& m) { mappings_ = m; }
    /// Mutable get transform
    inline MappingVector& getMappings() { return mappings_; }
    /// Get transform (non-mutable)
    inline const MappingVector& getMappings() const { return mappings_; }
    //@}

  protected:
    /// We estimate a different mapping for each coordinate and
    /// therefore store a vector of transformations
    MappingVector mappings_;
  }
  ; // end of class GridCell

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_GRIDCELL_H
