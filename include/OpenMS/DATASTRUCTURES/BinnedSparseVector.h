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
// $Id: BinnedSparseVector.h,v 1.2 2006/03/29 13:06:18 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_DATASTRUCTURES_BINNEDSPARSEVECTOR_H
#define OPENMS_DATASTRUCTURES_BINNEDSPARSEVECTOR_H

#include <map>

namespace OpenMS
{

  /**
   sparse Vector for use in BinnedRep <br>
   since the BinnedReps tend to be very sparse at low binsizes <br>
   this should use less space than a normal vector ad distance functions can just 
   ignore (hop()) zeroes, so it should be a little bit faster, too 
  */
  class BinnedSparseVector
  {
    friend class BinnedSparseVectorIterator;
    friend class BinnedSparseVectorConstIterator;
    friend class DoubleProxy;
    
    /**
     class DoubleProxy allows the BinnedSparseVector to differentiate between 
     writing and reading, so zeros can be ignored<br>
     see "more effective c++" section 30
     */
    class DoubleProxy
    {
    public:
      DoubleProxy(BinnedSparseVector& vec, uint index);
      DoubleProxy& operator=(const DoubleProxy& rhs);
      DoubleProxy& operator=(double val);
      operator double() const;
    private:
      BinnedSparseVector& vec_;
      int index_;
    };

    // forward decl
    class BinnedSparseVectorConstIterator;
    
    /**
    rudimentary iterator for BinnedSparseVector<br> 
    */
    class BinnedSparseVectorIterator
    {
      friend class BinnedSparseVector;
      friend class BinnedSparseVector::BinnedSparseVectorConstIterator;
    public:
      BinnedSparseVectorIterator(const BinnedSparseVectorIterator& source );
      ~BinnedSparseVectorIterator();

      BinnedSparseVectorIterator& operator++(); //präfix
      BinnedSparseVectorIterator operator++(int); //postfix
      
      // go to the next nonempty position
      BinnedSparseVectorIterator& hop();
      // find out at what position the iterator is
      // useful in combination with hop()
      uint position() const;
      
      DoubleProxy operator*();
      
      bool operator!=(const BinnedSparseVectorIterator& other);
      
    private:
      BinnedSparseVectorIterator();
      BinnedSparseVectorIterator(BinnedSparseVector& vector, int position);
      
      // the position in BinnedSparseVector
      uint position_;
     
      BinnedSparseVector& vector_;

      // the position in the underlying map of BinnedSparseVector
      std::map<uint,double>::const_iterator valit_;
    };
    
    /**
    rudimentary const_iterator for BinnedSparseVector<br> 
    */
    class BinnedSparseVectorConstIterator
    {
      friend class BinnedSparseVector;
    public:
      BinnedSparseVectorConstIterator(const BinnedSparseVectorConstIterator& source);
      BinnedSparseVectorConstIterator(const BinnedSparseVector::BinnedSparseVectorIterator& source);
      ~BinnedSparseVectorConstIterator();

      BinnedSparseVectorConstIterator& operator++();
      BinnedSparseVectorConstIterator operator++(int);

      // go to the next nonempty position
      BinnedSparseVectorConstIterator& hop();
      
      // find out at what position the iterator is
      // useful in combination with hop()
      uint position() const;
      
      double operator*();
      
      bool operator!=(const BinnedSparseVectorConstIterator& other);

    private:
      BinnedSparseVectorConstIterator();
      BinnedSparseVectorConstIterator(const BinnedSparseVector& vector, int position);
      
      // the position in BinnedSparseVector
      uint position_;
      
      const BinnedSparseVector& vector_;
      
      // the position in the underlying map of BinnedSparseVector
      std::map<uint,double>::const_iterator valit_;
    };

  public:
    
    typedef BinnedSparseVectorConstIterator const_iterator;
    typedef BinnedSparseVectorIterator iterator;
    
    /** @name constructor, destructor, assignment operator <br> */
    //@{
    BinnedSparseVector();
    BinnedSparseVector(int size );
    BinnedSparseVector(const BinnedSparseVector& source);
    
    // dtor
    ~BinnedSparseVector();

    BinnedSparseVector& operator=(const BinnedSparseVector& source);
    //@}
    // size of the represented vector
    uint size() const;

    // number of nonzero elements
    // i.e. the space actually used
    uint nonzero_size() const;
    
    // DoubleProxy handles the conversion to int and 
    // the writing ( if != 0 )
    const DoubleProxy operator[] (uint pos) const;
    DoubleProxy operator[] (uint pos);
    
    double at(uint pos) const ;
    void push_back(double value);

    void clear();
    void resize(uint newsize);

    iterator begin();
    iterator end();
    
    const_iterator begin() const;
    const_iterator end() const;
  private:
    std::map<uint, double> values_;
    uint size_;
  };
  
}

#endif //OPENMS_DATASTRUCTURES_BINNEDSPARSEVECTOR_H

