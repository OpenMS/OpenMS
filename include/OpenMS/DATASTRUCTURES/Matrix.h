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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_MATRIX_H
#define OPENMS_DATASTRUCTURES_MATRIX_H

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cmath> // pow()
#include <iomanip>
#include <vector>

namespace OpenMS
{

  /**
     @brief A two-dimensional matrix.  Similar to std::vector, but uses a binary
     operator(,) for element access.

		 Think of it as a random access container.  You can also generate gray
     scale images.  This is not designed to be used for linear algebra.

		 The following member functions of the base class std::vector<ValueType>
		 can also be used:

		 - begin
		 - end
		 - rbegin
		 - rend
		 - front
		 - back
		 - assign
		 - empty
		 - size
		 - capacity
		 - max_size
		 .

		 (It seems that Doxygen does not parse pure access declarations, so we
		 list them here.)

		 @ingroup Datastructures
  */
  template <typename Value>
	class Matrix : protected std::vector < Value >
	{
	 protected:
		typedef std::vector < Value > Base;

	 public:

		///@name STL compliance type definitions
		//@{
		typedef Base container_type;

		typedef typename Base::difference_type difference_type;
		typedef typename Base::size_type size_type;

		typedef typename Base::const_iterator const_iterator;
		typedef typename Base::const_reverse_iterator const_reverse_iterator;
		typedef typename Base::iterator iterator;
		typedef typename Base::reverse_iterator reverse_iterator;

		typedef typename Base::const_reference const_reference;
		typedef typename Base::pointer pointer;
		typedef typename Base::reference reference;
		typedef typename Base::value_type value_type;

		typedef typename Base::allocator_type allocator_type;
		//@}

		///@name OpenMS compliance type definitions
		//@{
		typedef Base ContainerType;
		typedef difference_type DifferenceType;
		typedef size_type SizeType;

		typedef const_iterator ConstIterator;
		typedef const_reverse_iterator ConstReverseIterator;
		typedef iterator Iterator;
		typedef reverse_iterator ReverseIterator;

		typedef const_reference ConstReference;
		typedef pointer Pointer;
		typedef reference Reference;
		typedef value_type ValueType;

		typedef allocator_type AllocatorType;
		//@}

		///@name Constructors, assignment, and destructor
		//@{
		Matrix ()
			: Base(),
				rows_(0),
				cols_(0)
		{}

		Matrix (SizeType rows, SizeType cols, ValueType value = ValueType())
			: Base(rows*cols,value),
				rows_(rows),
				cols_(cols)
		{}

		Matrix (const Matrix & source)
			: Base(source),
				rows_(source.rows_),
				cols_(source.cols_)
		{}

		Matrix & operator = (const Matrix & rhs)
		{
			Base::operator=(rhs);
			rows_ = rhs.rows_;
			cols_ = rhs.cols_;
			return *this;
		}

		~Matrix() {}
		//@}

		///@name Accessors
		//@{
		const_reference operator() (size_type const i, size_type const j) const
		{
			return getValue(i,j);
		}

		reference operator() (size_type const i, size_type const j)
		{
			return getValue(i,j);
		}

    const_reference getValue(size_type const i, size_type const j) const
		{
			return Base::operator[](index(i,j));
		}

    reference getValue(size_type const i, size_type const j)
		{
			return Base::operator[](index(i,j));
		}

    void setValue(size_type const i, size_type const j, value_type value)
		{
			Base::operator[](index(i,j)) = value;
		}

		//@}


		/**
			 @name Pure access declarations

			 These make begin(), end() etc. from container_type accessible.

			 @todo Fix check_test so that it will parse Matrix.h !  Workaround: To
			 see if Matrix_test.C is really complete, comment these access
			 declarations out.
		*/
		//@{
	 public:

		Base::begin;
		Base::end;
		Base::rbegin;
		Base::rend;

		Base::front;
		Base::back;
		Base::assign;

		Base::empty;
		Base::size;

		Base::capacity;
		Base::max_size;

		//@}

		void clear()
		{
			Base::clear();
			rows_ = 0;
			cols_ = 0;
		}

		void resize(size_type i, size_type j, value_type value = value_type())
		{
			rows_ = i;
			cols_ = j;
			Base::resize(rows_*cols_, value);
		}

		void resize(std::pair<Size,Size> const & size_pair, value_type value = value_type())
		{
			rows_ = size_pair.first;
			cols_ = size_pair.second;
			Base::resize(rows_*cols_, value);
		}

		/// Number of rows
		SizeType rows() const throw()
		{
			return rows_;
		}

		/// Number of columns
		SizeType cols() const throw()
		{
			return cols_;
		}

		std::pair<Size,Size> sizePair() const
		{
			return std::pair<Size,Size>(rows_,cols_);
		}

		/**@brief Calculate the index into the underlying vector from row and
			 column.  Note that Matrix uses the (row,column) lexicographic ordering
			 for indexing.
		*/
		SizeType const index(SizeType row, SizeType col) const
		{
#ifdef OPENMS_DEBUG
			if ( row >= rows_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,row, rows_);
			if ( col >= cols_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,col, cols_);
#endif
			return row * cols_ + col;
		}

		/**@brief Calculate the row and column from an index into the underlying
			 vector.  Note that Matrix uses the (row,column) lexicographic ordering
			 for indexing.
		*/
		std::pair<Size,Size> const indexPair(Size index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif
			return std::pair<SizeType,SizeType>(index/cols_,index%cols_);
		}

		/**@brief Calculate the column from an index into the underlying vector.
			 Note that Matrix uses the (row,column) lexicographic ordering for
			 indexing.
		*/
		SizeType colIndex(SizeType index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif
			return index%cols_;
		}

		/**@brief Calculate the row from an index into the underlying vector.
			 Note that Matrix uses the (row,column) lexicographic ordering for
			 indexing.
		*/
		SizeType rowIndex(SizeType index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif

			return index/cols_;
		}

		/**@brief Output content of Matrix as a grayscale image in plain PGM
			 format.


			 The abbreviation (and file extension) PGM stands for "portable
			 graymap".  PGM is one of the "portable pixmap" formats.  PGM was
			 designed to serve as a least-common-denominator for converting pixmap,
			 graymap, or bitmap files between different platforms.  It is a rather
			 stupid image format and very inefficient, but it has the advantage of
			 simplicity.  The "plain" variant of PGM even uses ASCII encoded decimal
			 numbers (as opposed to the raw format, which uses binary).  See
			 http://netpbm.sourceforge.net/doc/pgm.html or Wikipedia for further
			 information about PGM.

			 @param os The stream to write to.

			 @param maxval The maximal intensity of a pixel.  Values outside
			 [0:maxval] are truncated.  Note that 65535 ( == (1<<16)-1 ) is the
			 maximum for PGM (16 Bit), and values greater than 255 are not supported
			 by some older pieces of software.  The default is 0., which means to
			 set maxval to the maximum matrix entry, but of course never smaller
			 than 1.

			 @param scale A factor which is applied before rounding to gray levels.
			 The default is 0., which means automatic scaling such that the full
			 dynamic range is used.

			 @param reverse_video If @c true, the output is in <i>reverse video</i>.
			 Actually, this might be what you want in most cases - large matrix
			 entries show up dark.

			 @param gamma Parameter for the "gamma correction" (as known from image
			 processing) to be applied to the pixel intensities.  The gamma
			 correction replaces <i>x</i> with
			 <i>x<sup>gamma</sup>*maxval<sup>1-gamma</sup></i>.  (This is not
			 exactly the sRGB standard, which has a range with linear slope for
			 small intensities.)

			 @param comment A comment which will be embedded in the output. It can
			 consist of several lines separated by '\\n'.  By default, <i>no</i>
			 comment is written (not even an empty one).

			 Negative matrix entries are <i>never</i> mapped into the dynamic range.

			 Negative values for maxval or scale work like their absolute values.

			 Final remark: The ValueType must be "numeric" (or at least convertible
			 to double and int, comparable, etc.) for all this to work.

		*/
		std::ostream& writePGM ( std::ostream& os,
														 int maxval = 0,
														 double scale = 0,
														 double gamma = 1,
														 bool reverse_video = false,
														 std::string const& comment = String::EMPTY
													 ) const
		{

			if ( maxval < 0 )
			{
				maxval *= -1;
			}
			if ( scale < 0 )
			{
				scale *= -1;
			}

			// set automatic maxval and automatic scale, if requested
			if ( !maxval )
			{ // no maxval supplied, set default
				maxval = int( std::max( ValueType(1) , *std::max_element( begin(), end() ) ) );
				if ( !scale )
				{ // automatic scale
					scale = 1.0;
				}
			}
			else
			{ // maxval supplied
				if ( !scale )
				{ // automatic scale
					ValueType max_entry = *std::max_element( begin(), end() );
					if ( max_entry <= ValueType(0) )
					{ // All entries non-positive!  These are outside the dynamic range, so we can dump zeroes as well...
						scale = 0;
					}
					else
					{
						scale =  double(maxval) / double( max_entry);
					}
				}
			}

			// write PGM header
			os <<
				"P2\n"
				"# columns rows\n"
				 << cols() << ' ' << rows() << "\n"
				"# maxval\n"
				 << maxval << '\n' <<
				"# scaling factor:    " << scale << "\n"
				"# gamma correction:  " << gamma << ( gamma == 1. ? " (none)\n" : "\n") <<
				"# reverse video:     " << ( reverse_video ? "on\n" : "off\n" )
				;
			// Write out the comment, with "# " before each line.
			if ( ! comment.empty() )
			{
				os << "#----------\n";
				String quoted_comment(comment);
				std::vector<String> pieces;
				if ( quoted_comment.split('\n',pieces) )
				{
					quoted_comment.implode(pieces.begin(),pieces.end(),"\n# ");
					os << "# " << quoted_comment << '\n';
				}
				else
				{
					os << "# " << comment << '\n';
				}
				os << "#----------\n";
			}

			// write data
			const_iterator iter = begin();
			// Number of columns (per line) in the output file.  Note that we must
			// not exceed 70 chars per line, according to the specification of the
			// PGM format. Since 65535 takes six chars, ten seems a good choice.
			const unsigned int cols_output = 10;

			// We'll have similar stuff for normal and reverse video.  For ease of
			// maintenance I use a macro here to keep both versions in sync.
			// ((( Template programming, the HARD way ;-) )))
#define PGM_DATA_LOOP \
			for ( size_type i = 0; i < rows(); ++i )															\
			{																																			\
				if ( cols() > cols_output )																					\
				{																																		\
					os << "# row " << i << '\n';																			\
				}																																		\
				for ( Size j = 0, count = 0; ; )																		\
				{																																		\
					/* output gray value, rounded to nearest int */										\
					int gray = int( *iter * scale + 0.5 );														\
					if ( gray < 0 )																										\
					{																																	\
						os << PGM_REVERSE_VIDEO(0,maxval);															\
					}																																	\
					else if ( gray > maxval )																					\
					{																																	\
						os << PGM_REVERSE_VIDEO(maxval,0);															\
					}																																	\
					else																															\
					{																																	\
						os << PGM_GAMMA_CORRECTED(PGM_REVERSE_VIDEO(gray,maxval-gray));	\
					}																																	\
					/* loop increment */																							\
					++j; ++iter;																											\
					/* loop condition */																							\
					if ( j < cols()	)																									\
					{																																	\
						/* output whitespace */																					\
						if ( ++count == cols_output )																		\
						{																																\
							os << '\n';																										\
							count = 0;																										\
						}																																\
						else																														\
						{																																\
							os << ' ';																										\
						}																																\
					}																																	\
					else break;																												\
				}																																		\
				os << '\n';																													\
			}

			// Now here comes the "real" code.
			if ( gamma == 1. )
			{ // no gamma correction
#define PGM_GAMMA_CORRECTED(x) (x)
				if ( reverse_video )
				{
#define PGM_REVERSE_VIDEO(a,b) (b)
					PGM_DATA_LOOP;
#undef PGM_REVERSE_VIDEO
				}
				else
				{
#define PGM_REVERSE_VIDEO(a,b) (a)
					PGM_DATA_LOOP;
#undef PGM_REVERSE_VIDEO
				}
#undef PGM_GAMMA_CORRECTED
			}
			else
			{ // apply gamma correction
				double const gamma_reciprocal = 1./gamma;
				double const rescaling_factor = std::pow(maxval,1.-gamma_reciprocal);
#define PGM_GAMMA_CORRECTED(x) (int(std::pow( (x) ,gamma_reciprocal)*rescaling_factor + 0.5))
				if ( reverse_video )
				{
#define PGM_REVERSE_VIDEO(a,b) (b)
					PGM_DATA_LOOP;
#undef PGM_REVERSE_VIDEO
				}
				else
				{
#define PGM_REVERSE_VIDEO(a,b) (a)
					PGM_DATA_LOOP;
#undef PGM_REVERSE_VIDEO
				}
#undef PGM_GAMMA_CORRECTED
			}
#undef PGM_DATA_LOOP

			return os;
		}

		/**@brief Equality comparator.

		If matrices have different row or colmn numbers, throws a precondition exception.
		*/
		bool operator == ( Matrix const & rhs ) const
			throw (Exception::Precondition)
		{
			OPENMS_PRECONDITION(cols_ == rhs.cols_,
													"Matrices have different row sizes.");
			OPENMS_PRECONDITION(rows_ == rhs.rows_,
													"Matrices have different column sizes.");
			return static_cast < typename Matrix<Value>::Base const &>(*this) == static_cast < typename Matrix<Value>::Base const &>(rhs);
		}

		/**@brief Less-than comparator.  Comparison is done lexicographically: first by row, then by column.

		If matrices have different row or column numbers, throws a precondition exception.
		*/
		bool operator < (Matrix const & rhs) const
			throw (Exception::Precondition)
		{
			OPENMS_PRECONDITION(cols_ == rhs.cols_,
													"Matrices have different row sizes.");
			OPENMS_PRECONDITION(rows_ == rhs.rows_,
													"Matrices have different column sizes.");
			return static_cast < typename Matrix<Value>::Base const &>(*this) < static_cast < typename Matrix<Value>::Base const &>(rhs);
		}

	 protected:

		///@name Data members
		//@{
		/// Number of rows (height of a column)
		SizeType rows_;
		/// Number of columns (width of a row)
		SizeType cols_;
		//@}

	}; // class Matrix

	/**@brief Print the contents to a stream.

	@relatesalso Matrix
	*/
	template <typename Value>
	std::ostream& operator << (std::ostream& os, const Matrix<Value>& matrix)
	{
		typedef typename Matrix<Value>::size_type size_type;
	  for ( size_type i = 0; i < matrix.rows(); ++i ) {
	    for ( size_type j = 0; j < matrix.cols(); ++j ) {
				os << std::setprecision(6) << std::setw(6) << matrix(i,j) << ' ';
	    }
	    os << std::endl;
	  }
		return os;
	}


} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_MATRIX_H
