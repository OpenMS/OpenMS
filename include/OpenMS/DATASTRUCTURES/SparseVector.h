// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_DATASTRUCTURES_SPARSEVECTOR_H
#define OPENMS_DATASTRUCTURES_SPARSEVECTOR_H

#include <vector>

namespace OpenMS
{
  /** @brief sparse vector implementation which does noch save zero-elements
	
  	sparse Vector for internal use in Similarity Matrices from ordered Spectra 
  	Comparison<br>
  	since only Spectra with similar parent_mass get similarity > 0 all other 
  	similarities dont need to be saved<br>
  	the resulting matrix is a banded matrix with varying bandwith <br> <br>

		@ingroup Datastructures 
  */
  class SparseVector 
  {
		// @name Friends
		// @{
    friend class DoubleProxy;
		// @}
		
    // @name Typedefs
		// @{
    typedef float valuetype;
		// @}

    /**
		   @brief proxy class that relays everything but zeroes to the vector <br>
    */
    class DoubleProxy
    {
    public:
      // @name Constructors and Destructors
      // @{
			/// detailed constructor
      DoubleProxy(SparseVector& vec, uint index);
			// @}
			
			// @name Operators
			// @{
			/// assignment operator
      DoubleProxy& operator=(const DoubleProxy& rhs);

			/// assignment operator
      DoubleProxy& operator=(SparseVector::valuetype val);
      
			/// return valuetype (see vector stl docs)
      operator SparseVector::valuetype() const;
			// @}
			
    private:
			
			///
      SparseVector& vec_;

			///
      uint index_;
    };

  public:
    
    // @name Constructors and Destructors
		// @{
		/// detailed constructor
    SparseVector(uint,uint);

		/// default constructor
    SparseVector();

		/// copy constructor
    SparseVector(const SparseVector& source);

		/// destructor
    virtual ~SparseVector();
		// @}

		// @name Operators
		// @{
		/// assignment operator
    SparseVector& operator=(const SparseVector& source);

		///
		const DoubleProxy operator[] (uint pos) const;

		///
		DoubleProxy operator[](uint);
    //@}

		// @name Accessors
		// @{
    //return type is distance of inserted double from the next existing value (to monitor the growth of the SparseVector
    uint insert(uint, valuetype);
		// @}
  
    // @name Internal
    // @{
		///
    void growfront(uint);

		///
    void growback(uint);

		///
    uint firstentry() { return firstentry_; }
		
		///
    uint nonzero_size() { return leftarray_.size() + rightarray_.size(); }
		
    /// fraction of saved nonzeroes
    double used() const;
    // @}
		
  private:

    uint firstentry_;
    uint middle_;
    std::vector<valuetype> rightarray_;
    std::vector<valuetype> leftarray_;
  };
}
#endif //OPENMS_DATASTRUCTURES_SPARSEVECTOR_H
