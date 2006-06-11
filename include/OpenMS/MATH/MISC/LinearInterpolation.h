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
// $Id: LinearInterpolation.h,v 1.5 2006/03/28 16:19:59 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_MISC_LINEARINTERPOLATION_H
#define OPENMS_MATH_MISC_LINEARINTERPOLATION_H

#ifndef STD_VECTOR
#define STD_VECTOR
#include <vector>
#endif

namespace OpenMS
{  
  /**
  	@brief Smooth (interpolate) an isotope distribution or peak shape.
  	
  	You may ask for value() and derivative().
  	Interpolation is simply linear.  Original values beyond the boundaries are implicitly taken as zero.

  	@todo add to Math namespace
  	
  	@ingroup Math
  */
  template < typename Key = double, typename Value = Key >
  class LinearInterpolation
  {

  public:

		///\name Typedefs
		//@{
    typedef Value value_type;
    typedef Key key_type;
    typedef std::vector < value_type > container_type;

		typedef value_type      ValueType;
		typedef key_type        KeyType;
		typedef container_type  ContainerType;
		//@}

	 public:

    /**@brief Constructor. The first argument specifies the random access
			 container (valid keys are 0,1,...,size()-1) from which interpolated
			 values will be sampled.  The second argument is the scale which is
			 applied to the arguments of value() and derivative() before looking up
			 the interpolated values in the container.  The third arguments is the
			 offset, which is subtracted before everything else.
		 */
    LinearInterpolation ( key_type _scale = 1., key_type _offset = 0. )
      : scale_  ( _scale ),
				offset_ ( _offset ),
        data_   ( )
    {}

		/// Copy constructor.
		LinearInterpolation ( LinearInterpolation const & _arg )
      : scale_  ( _arg.scale_ ),
				offset_ ( _arg.offset_ ),
        data_   ( _arg.data_ )
    {}

		/// Destructor.
		~LinearInterpolation () {}

		
		/// The transformation from "outside" to "inside" coordinates.
		inline key_type key2index ( key_type pos ) const throw()
		{
			pos -= offset_;
			pos /= scale_;
			return pos;
		}

		/// The transformation from "inside" to "outside" coordinates.
		inline key_type index2key ( key_type pos ) const throw()
		{
			pos *= scale_;
			pos += offset_;
			return pos;
		}


    /// Returns the interpolated value.
    value_type
    value ( key_type _pos ) const throw()
    {

			// apply the key transformation
			key_type const pos = key2index(_pos);

			int const size_ = data_.size();
      int const left = int(pos); // this rounds towards zero

      if ( pos <= 0 )
        if ( left != 0 )
          return 0;
        else  // that is: -1 < pos <= 0
          return data_[ 0 ] * ( 1. + pos ) ;
      
      if ( left >= size_ - 1 )
        if ( left != size_ - 1 )
          return 0;
        else
          return data_[ left ] * ( size_ - pos );

      key_type factor = pos - key_type(left);

      return // weighted average
        data_[ left + 1 ] * factor +
        data_[ left ] * ( 1. - factor );
    }


    /// Returns the interpolated derivative.
    value_type
    derivative ( key_type _pos ) const throw()
    {

			// apply the key transformation
			key_type const pos = key2index(_pos);
			
			int const size_ = data_.size();
      int const left = int(pos+0.5); // rounds towards zero
      
      if ( left < 0 ) // quite small
        return 0;
      else 
        if ( left == 0 ) // at the border
          if ( pos >= -0.5 ) // that is: -0.5 <= pos < +0.5
            return 
              ( data_[1] - data_[0] ) * ( pos + 0.5 ) +
              ( data_[0] ) * ( 0.5 - pos );
          else // that is: -1.5 <= pos < -0.5
            return
              ( data_[0] ) * ( pos + 1.5 );

      // "else" case: to the right of the left margin


      key_type factor = key_type(left) - pos + 0.5;

      if ( left > size_ ) // quite large
        return 0;
      else
        if ( left < size_ - 1 ) // to the left of the right margin
          return // weighted average of derivatives for adjacent intervals
            ( data_[left] - data_[left-1] ) * factor +
            ( data_[left+1] - data_[left] ) * ( 1. - factor );
        else // somewhat at the border
          if ( left == size_ - 1 ) // at the border, first case
            return 
              ( data_[left] - data_[left-1] ) * factor +
              ( - data_[left] ) * ( 1. - factor );
      // else // that is: left == size_

      // We pull the last remaining case out of the "if" tree to avoid a
      // compiler warning ...

      return // at the border, second case
        ( - data_[left-1] ) * factor;
      
    }
		

		///\name Accessors
		//@{
		/// Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".
		inline key_type & getScale () throw() { return scale_; }
		/// Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".
		inline key_type const & getScale () const throw(){ return scale_; }
		/// Accessor.  "Scale" is the difference (in "outside" units) between consecutive entries in "Data".
		inline void setScale ( key_type const & _scale ) throw() { getScale() = _scale; }

		/// Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data[0]".
		inline key_type & getOffset () throw() { return offset_; }
		/// Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data[0]".
		inline key_type const & getOffset () const throw() { return offset_; }
		/// Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data[0]".
		inline void setOffset ( key_type const & _offset ) throw() { getOffset() = _offset; }

		/// Accessor.
		inline container_type & getData () throw() { return data_; }
		/// Accessor.
		inline container_type const & getData () const throw() { return data_; }
		/// Accessor.
		inline void setData ( container_type const & _data ) throw() { getData() = _data; }
		//@}


		/// Returns \c true if getData() is empty.
		inline bool empty () const throw() { return data_.empty(); }


		/// Lower boundary of the support, in "outside" coordinates.
		inline key_type supportMin() const throw()
		{ return index2key ( empty() ? 0 : -1 ); }

		/// Upper boundary of the support, in "outside" coordinates.
		inline key_type supportMax() const throw()
		{ return index2key ( data_.size() ); }

  private:
    
		key_type scale_;
		key_type offset_;
    container_type data_;
    
  };
  
} // namespace OpenMS

#endif // OPENMS_MATH_MISC_LINEARINTERPOLATION_H
