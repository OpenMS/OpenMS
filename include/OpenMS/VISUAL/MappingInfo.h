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

#ifndef OPENMS_VISUAL_MAPPINGINFO_H
#define OPENMS_VISUAL_MAPPINGINFO_H

namespace OpenMS 
{
	class Param;
	
	/**
		@brief Describes the mapping of data to the screen.
		
		A class that describes the mapping (which dimension is x/y) and orientation (ascending/descending) 
		of the two dimensions x and y on the screen. Orientations are always seen from the lower left corner.
		
		Default is: m/z to x axis , orientation for x axis ascending, orientation for y axis ascending
	
		@ingroup Visual
	*/
  class MappingInfo
	{
		public:
			/// Default constructor
			inline MappingInfo():
				x_ascending_(true), 
				y_ascending_(true), 
				mz_to_x_(true) 
			{
				
			}
	
			/// set mapping of m/z to X-Axis
			inline void setMzToXAxis()  
			{	
				mz_to_x_ = true;
			}
			/// set mapping of m/z to Y-Axis
			inline void setMzToYAxis()	
			{	
				mz_to_x_ = false;
			}
	
			/// set orientation of X Axis (left -> right)
			inline void setXAxisAsc() 
			{	
				x_ascending_ = true;
			}
			/// set orientation of X Axis (left <- right)
			inline void setXAxisDesc()	
			{
				x_ascending_ = false;
			}
			/// set orientation of Y Axis (top <- bottom)
			inline void setYAxisAsc()	
			{	
				y_ascending_ = true;
			}
			/// set orientation of Y Axis (top -> bottom)
			inline void setYAxisDesc()	
			{	
				y_ascending_ = false;
			}
			
			/// mapping predicate
			inline bool isMzToXAxis() const  
			{  
				return mz_to_x_; 
			}
			/// orientation X axis predicate
			inline bool isXAxisAsc() const 
			{ 
				return x_ascending_; 
			}
			/// orientation Y axis predicate
			inline bool isYAxisAsc() const 
			{  
				return y_ascending_;
			}
			
			///returns the current settings 
			Param getParam();
			///sets the current settings
			void setParam(const Param& p);

		private:
			bool x_ascending_;
			bool y_ascending_;
			bool mz_to_x_;
	};

} //namespace

#endif
