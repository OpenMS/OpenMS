// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: SparseVector.h,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_DATASTRUCTURES_SVECTOR_H
#define OPENMS_DATASTRUCTURES_SVECTOR_H

#include <vector>

namespace OpenMS
{
  /**
  sparse Vector for internal use in Similarity Matrices from ordered Spectra 
  Comparison<br>
  since only Spectra with similar parent_mass get similarity > 0 all other 
  similarities dont need to be saved<br>
  the resulting matrix is a banded matrix with varying bandwith <br> <br>
  **/
  class SparseVector 
  {
    friend class DoubleProxy;
      
    typedef float valuetype;

    /**
    proxy class that relays everything but zeroes to the vector <br>
    */
    class DoubleProxy
    {
    public:
      /** @name constructor, assignment operator <br> */
      //@{
      DoubleProxy(SparseVector& vec,uint index);
      DoubleProxy& operator=(const DoubleProxy& rhs);
      DoubleProxy& operator=(SparseVector::valuetype val);
      //@}
      operator SparseVector::valuetype() const;
    private:
      SparseVector& vec_;
      uint index_;
    };

  public:
    
    /** @name constructors, destructor, assignment operator <br> */
    //@{
    SparseVector(uint,uint);
    SparseVector();
    SparseVector(const SparseVector& source);
    ~SparseVector();
    SparseVector& operator=(const SparseVector& source);
    //@}
    //return type is distance of inserted double from the next existing value (to monitor the growth of the SparseVector
    /** \deprecated */
    uint insert(uint,valuetype);
  
    const DoubleProxy operator[] (uint pos) const;
    DoubleProxy operator[](uint);
    
    /** @name internal/debug stuff <br> */
    //@{
    void growfront(uint);
    void growback(uint);
    uint firstentry() { return firstentry_;}
    uint nonzero_size() { return leftarray_.size() + rightarray_.size();}
    // fraction of saved nonzeroes
    double used() const;
    //@}
  private:
    uint firstentry_;
    uint middle_;
    std::vector<valuetype> rightarray_;
    std::vector<valuetype> leftarray_;
  };
}
#endif //OPENMS_DATASTRUCTURES_SPARSEVECTOR_H
