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

#include <OpenMS/CONCEPT/Types.h>

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
  class SparseVector
  {
		// @name Friends
		// @{
    friend class SparseVectorIterator;
    friend class SparseVectorConstIterator;
    friend class DoubleProxy;
		// @}
    
    /**
			@brief class DoubleProxy allows the SparseVector to differentiate between writing and reading, so zeros can be ignored
			
			See "more effective c++" section 30
    */
    class DoubleProxy
    {
		
	    public:

				// @name Constructors and Destructors
				// @{
				///
  	    DoubleProxy(SparseVector& vec, UInt index);
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
      	SparseVector& vec_;
				
				///
      	int index_;
    };

    // forward declaration
    class SparseVectorConstIterator;
    
    /// iterator for SparseVector
    class SparseVectorIterator
    {
			// @name Friends
			// @{
      friend class SparseVector;
      friend class SparseVector::SparseVectorConstIterator;
			// @}
			
    public:

			// @name Constructors and Destructors
			// @{
			/// copy constructor
      SparseVectorIterator(const SparseVectorIterator& source);
			
			/// destructor
      virtual ~SparseVectorIterator();
			// @}

			// @name Operators
			// @{
			///
      SparseVectorIterator& operator ++ (); // prefix

			///
      SparseVectorIterator operator ++ (int); // postfix

			///
			DoubleProxy operator * ();			
			// @}
      
		
			// @name Accessors
			// @{
      /// go to the next nonempty position
      SparseVectorIterator& hop();
			
      /// find out at what position the iterator is; useful in combination with hop()
      UInt position() const;
			// @}
      
			// @name Predicates
			// @{
			/// inequality operator
      bool operator!=(const SparseVectorIterator& other);
			// @}
      
    private:

			/// default constructor
      SparseVectorIterator();

			/// 
      SparseVectorIterator(SparseVector& vector, int position);
      
      // the position in SparseVector
      UInt position_;
     
		 	///
      SparseVector& vector_;

      // the position in the underlying map of SparseVector
      std::map<UInt,double>::const_iterator valit_;
    };
    
    /// const_iterator for SparseVector
    class SparseVectorConstIterator
    {
			// @name Friends
			// @{
      friend class SparseVector;
			// @}
			
    public:
	
			// @name Constructors and Destructors
			// @{
			/// copy constructor
      SparseVectorConstIterator(const SparseVectorConstIterator& source);
			
			/// copy constructor from SparseVector::SparseVectorIterator
      SparseVectorConstIterator(const SparseVector::SparseVectorIterator& source);
			
			/// destructor
      virtual ~SparseVectorConstIterator();
			// @}

			// @name Operators
			// @{
			/// postincrement operator
      SparseVectorConstIterator& operator ++ ();

			/// immidiate increment operator
      SparseVectorConstIterator operator ++ (int);

			/// derefence operator
			double operator * ();
			// @}

			// @name Accessors
			// @{
      /// go to the next nonempty position
      SparseVectorConstIterator& hop();
      
      /// find out at what position the iterator is, useful in combination with hop()
      UInt position() const;
			// @} 
			
			// @name Predicates
			// @{
			/// inequality operator
      bool operator != (const SparseVectorConstIterator& other);
			// @}

    private:

			/// default constructor
      SparseVectorConstIterator();

			/// detailed constructor
      SparseVectorConstIterator(const SparseVector& vector, int position);
      
      // the position in SparseVector
      UInt position_;
      
			/// 
      const SparseVector& vector_;
      
      // the position in the underlying map of SparseVector
      std::map<UInt, double>::const_iterator valit_;
    };

  public:
    
		// @name Typedefs
		// @{
    typedef SparseVectorConstIterator const_iterator;
		typedef SparseVectorConstIterator ConstIterator;
    typedef SparseVectorIterator iterator;
		typedef SparseVectorIterator Iterator;
		// @}
    
    /// @name Constructor and Desctructor
    // @{
		/// default constructor
    SparseVector();

		/// detailed constructor
    SparseVector(int size);

		/// copy constructor
    SparseVector(const SparseVector& source);
    
    /// destructor
    virtual ~SparseVector();
		// @}

		// @name Operators
		// @{
		/// assignment operator
    SparseVector& operator = (const SparseVector& source);

		/// DoubleProxy handles the conversion to int and ,the writing ( if != 0 )
		const DoubleProxy operator[] (UInt pos) const;

		///
		DoubleProxy operator[] (UInt pos);
    // @}

	
		// @name Accessors
		// @{
    /// size of the represented vector
    UInt size() const;
		
    /// number of nonzero elements, i.e. the space actually used
    UInt nonzero_size() const;
    
		/// at (see stl vector docs)
    double at(UInt pos) const;
		
		/// push_back (see stl vector docs)
    void push_back(double value);

		/// removes all elements
    void clear();
		
		/// resizes the the vector to @param newsize
    void resize(UInt newsize);
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
    std::map<UInt, double> values_;

		///
    UInt size_;
  };
  
}

#endif //OPENMS_DATASTRUCTURES_SPARSEVECTOR_H

