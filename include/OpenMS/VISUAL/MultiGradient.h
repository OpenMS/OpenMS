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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_MULTIGRADIENT_H
#define OPENMS_VISUAL_MULTIGRADIENT_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

//QT
#include <qcolor.h>

//STL
#include <map>
#include <vector>

namespace OpenMS {

/**
	@brief A gradient of multiple colors and arbitrary distances between colors.
	
	Positions associated with numbers range from 0 to 100.
	There is always a color associated with position 0 and 100.
	Stretching the gradient to a specified range, and precalculation and 
	caching is also possible.

	@ingroup Visual
*/
class MultiGradient
{
	public:

	/// Interploation mode. 
	enum InterpolationMode {
		IM_LINEAR, 	///< IM_LINEAR returns the linear interploation (default). 
		IM_STAIRS   ///< IM_STAIRS returns the color of the next lower position
		};
	
	///Constructor
	MultiGradient();
	///Destructor
	~MultiGradient();

	/// sets or replaces the color at position @p position
	void insert (SignedInt position, const QColor& color);
	/// removes the color at position @p position. Colors at positions 0 and 100 cannot be removed.
	bool remove (SignedInt position);
	/// returns if a value for position @p position exists
	bool exists (SignedInt position);
	/// returns the position of the @p index -th point
	UnsignedInt position(UnsignedInt index) throw (Exception::IndexUnderflow,Exception::IndexOverflow);
	/// returns the color of the @p index -th point
	const QColor& color(UnsignedInt index) throw (Exception::IndexUnderflow,Exception::IndexOverflow);


	/** 
		@brief Returns the color as @p position.
		
		If the @p position is higher or lower than the range [0,100] the highest, 
		respectively the lowest, color is returned.
	*/
	QColor interpolatedColorAt(double position) const;
	/** 
		@brief returns the color as @p position with the gradient stretched between @p min and @p max.
		 
		If the @p position is higher or lower than the range [min,max] the highest, 
		respectively the lowest, color is returned.
	*/
	QColor interpolatedColorAt(double position, double min, double max) const;

	/// activates the precalculation of values (only approximate results are given)
	void activatePrecalculationMode(double min, double max, UnsignedInt steps);
	/// deactivates the precalculation of values ( and deletes the precalculated values)
	void deactivatePrecalculationMode();
	/** 
		@brief Returns a precalculated color. 
		
		If the @p position is higher or lower than the the range
		specified in activatePrecalculationMode(...) the highest, respectively the lowest, color is returned.
		If precalcualtion mode is not activated, the exception is thrown.
	*/
	const QColor& precalculatedColorAt(double position) const throw (Exception::OutOfSpecifiedRange);

	///return the number of color points
	UnsignedInt size() const;

	/// sets the interploation mode (default or stairs). Default is linear
	void setInterpolationMode(UnsignedInt mode);
	/// returns the interpolaion mode
	UnsignedInt getInterpolationMode() const;

	///convert to string representation
	std::string toString() const;
	/**
		@brief Sets the gradient by string representation.
		
		The string represenation of a gradient starts with the interpolation mode: "Linear" or "Stairs" and the separator "|".
		It is followed by an arbitrary number of integer-color-pairs. 
		
		Such a pair consists of an integer (0-100) followed by a comma and
		a "#". Then follows a color in RGB notation "#RRGGBB" and finally a semicolon.
		
		Examples are:
		<UL>
			<LI> "Linear|0,#ffff00;100,#000000"
			<LI> "Stairs|0,#ffff00;11,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000"
		</UL>
	*/
	void fromString(const std::string& gradient);

	protected:
	/// Map of index and color
	std::map<UnsignedInt,QColor> pos_col_;
	/// Current interpolation mode
	UnsignedInt interpolation_mode_;
	/// Precalculated colors
	std::vector<QColor> pre_;
	/// Minimum of the precalculated color range
	double pre_min_;
	/// Size of the precalculated color range
	double pre_size_;
	/// Steps of the precalculated color range
	SignedInt pre_steps_;

};

}
#endif // OPENMS_VISUAL_MULTIGRADIENT_H

