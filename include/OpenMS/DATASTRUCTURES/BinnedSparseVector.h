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
#ifndef OPENMS_DATASTRUCTURES_BINNEDSPARSEVECTOR_H
#define OPENMS_DATASTRUCTURES_BINNEDSPARSEVECTOR_H

#include <map>

namespace OpenMS
{

  /** @brief binned sparse vector implementation, which does not contain zero-elements
	
   		sparse Vector for use in BinnedRep <br>
   		since the BinnedReps tend to be very sparse at low binsizes <br>
   		this should use less space than a normal vector ad distance functions can just 
   		ignore (hop()) zeroes, so it should be a little bit faster, too 

			@ingroup Datastructures
  */
  class BinnedSparseVector
  {
		// @name Friends
		// @{
    friend class BinnedSparseVectorIterator;
    friend class BinnedSparseVectorConstIterator;
    friend class DoubleProxy;
		// @}
    
    /**
			@brief class DoubleProxy allows the BinnedSparseVector to differentiate between writing and reading, so zeros can be ignored
			
			See "more effective c++" section 30
    */
    class DoubleProxy
    {
		
	    public:

				// @name Constructors and Destructors
				// @{
				///
  	    DoubleProxy(BinnedSparseVector& vec, uint index);
				// @}

				// @name Operators
				// @{
				/// 
   	 		DoubleProxy& operator=(const DoubleProxy& rhs);
				
				///
   	  	DoubleProxy& operator=(double val);
				
				///
      	operator double() const;
				// @}
				
    	private:
			
				///
      	BinnedSparseVector& vec_;
				
				///
      	int index_;
    };

    // forward declaration
    class BinnedSparseVectorConstIterator;
    
    /// iterator for BinnedSparseVector
    class BinnedSparseVectorIterator
    {
			// @name Friends
			// @{
      friend class BinnedSparseVector;
      friend class BinnedSparseVector::BinnedSparseVectorConstIterator;
			// @}
			
    public:

			// @name Constructors and Destructors
			// @{
			/// copy constructor
      BinnedSparseVectorIterator(const BinnedSparseVectorIterator& source);
			
			/// destructor
      virtual ~BinnedSparseVectorIterator();
			// @}

			// @name Operators
			// @{
			///
      BinnedSparseVectorIterator& operator ++ (); // prefix

			///
      BinnedSparseVectorIterator operator ++ (int); // postfix

			///
			DoubleProxy operator * ();			
			// @}
      
		
			// @name Accessors
			// @{
      /// go to the next nonempty position
      BinnedSparseVectorIterator& hop();
			
      /// find out at what position the iterator is; useful in combination with hop()
      uint position() const;
			// @}
      
			// @name Predicates
			// @{
			/// inequality operator
      bool operator!=(const BinnedSparseVectorIterator& other);
			// @}
      
    private:

			/// default constructor
      BinnedSparseVectorIterator();

			/// 
      BinnedSparseVectorIterator(BinnedSparseVector& vector, int position);
      
      // the position in BinnedSparseVector
      uint position_;
     
		 	///
      BinnedSparseVector& vector_;

      // the position in the underlying map of BinnedSparseVector
      std::map<uint,double>::const_iterator valit_;
    };
    
    /// const_iterator for BinnedSparseVector
    class BinnedSparseVectorConstIterator
    {
			// @name Friends
			// @{
      friend class BinnedSparseVector;
			// @}
			
    public:
	
			// @name Constructors and Destructors
			// @{
			/// copy constructor
      BinnedSparseVectorConstIterator(const BinnedSparseVectorConstIterator& source);
			
			/// copy constructor from BinnedSparseVector::BinnedSparseVectorIterator
      BinnedSparseVectorConstIterator(const BinnedSparseVector::BinnedSparseVectorIterator& source);
			
			/// destructor
      virtual ~BinnedSparseVectorConstIterator();
			// @}

			// @name Operators
			// @{
			/// postincrement operator
      BinnedSparseVectorConstIterator& operator ++ ();

			/// immidiate increment operator
      BinnedSparseVectorConstIterator operator ++ (int);

			/// derefence operator
			double operator * ();
			// @}

			// @name Accessors
			// @{
      /// go to the next nonempty position
      BinnedSparseVectorConstIterator& hop();
      
      /// find out at what position the iterator is, useful in combination with hop()
      uint position() const;
			// @} 
			
			// @name Predicates
			// @{
			/// inequality operator
      bool operator != (const BinnedSparseVectorConstIterator& other);
			// @}

    private:

			/// default constructor
      BinnedSparseVectorConstIterator();

			/// detailed constructor
      BinnedSparseVectorConstIterator(const BinnedSparseVector& vector, int position);
      
      // the position in BinnedSparseVector
      uint position_;
      
			/// 
      const BinnedSparseVector& vector_;
      
      // the position in the underlying map of BinnedSparseVector
      std::map<uint, double>::const_iterator valit_;
    };

  public:
    
		// @name Typedefs
		// @{
    typedef BinnedSparseVectorConstIterator const_iterator;
		typedef BinnedSparseVectorConstIterator ConstIterator;
    typedef BinnedSparseVectorIterator iterator;
		typedef BinnedSparseVectorIterator Iterator;
		// @}
    
    /// @name Constructor and Desctructor
    // @{
		/// default constructor
    BinnedSparseVector();

		/// detailed constructor
    BinnedSparseVector(int size);

		/// copy constructor
    BinnedSparseVector(const BinnedSparseVector& source);
    
    /// destructor
    virtual ~BinnedSparseVector();
		// @}

		// @name Operators
		// @{
		/// assignment operator
    BinnedSparseVector& operator = (const BinnedSparseVector& source);

		/// DoubleProxy handles the conversion to int and ,the writing ( if != 0 )
		const DoubleProxy operator[] (uint pos) const;

		///
		DoubleProxy operator[] (uint pos);
    // @}

	
		// @name Accessors
		// @{
    /// size of the represented vector
    uint size() const;
		
    /// number of nonzero elements, i.e. the space actually used
    uint nonzero_size() const;
    
		/// at (see stl vector docs)
    double at(uint pos) const;
		
		/// push_back (see stl vector docs)
    void push_back(double value);

		/// removes all elements
    void clear();
		
		/// resizes the the vector to @param newsize
    void resize(uint newsize);
		// @}

		/// @name Iterators
		//@{
		/// begin iterator
    iterator begin();
		
		/// end iterator
    iterator end();
    
		/// const begin iterator 
    const_iterator begin() const;
		
		/// const end iterator
    const_iterator end() const;
		//@} 
		
  private:
	
		///
    std::map<uint, double> values_;

		///
    uint size_;
  };
  
}

#endif //OPENMS_DATASTRUCTURES_BINNEDSPARSEVECTOR_H

