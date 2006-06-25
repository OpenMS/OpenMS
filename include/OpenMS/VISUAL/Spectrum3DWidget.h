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
// $Id: Spectrum3DWidget.h,v 1.13 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM3DWIDGET_H

#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/config.h>

namespace OpenMS
{	
	class Spectrum3DCanvas;
	/**
		@brief Widget for 3D-visualization of map data
		
		
		
		@ingroup spectrum_widgets
	*/
	class Spectrum3DWidget:public SpectrumWidget
	{
		Q_OBJECT		

	public:	
		/**
		 *@brief Constructor
		 *
		 *Spectrum3DWidget Constructor
		 *@param parent The parent Widget
		 *@param name The Widget's name
		 *@param f Widget flags
		 *
		 */
		Spectrum3DWidget(QWidget* parent = 0, const char* name = "Spectrum3DWidget", WFlags f = 0);
		/**
		 *
		 *@brief Destructor
		 *
		 *Destroys the Widget and all assosiated data
		 */
		virtual ~Spectrum3DWidget();
		
		enum DotModes 
			{
				DOT_BLACK = 0,            ///< use black only
				DOT_GRADIENT = 1          ///< use gradient
			};
		
    enum ShadeModes 
      {
				SHADE_FLAT = 0,            
				SHADE_SMOOTH = 1         
      };
		enum IntScale 
      {
				INT_LINEAR = 0,            
				INT_LOG = 1         
      };
		/**
		 *@brief returns the Canvas Widget
		 *Returns the Canvas Widget
		 *@return the Canvas Widget
		 */
		Spectrum3DCanvas* canvas() const;
		
		/**
		 *@brief Creates a preferences dialog Page
		 *
		 *Creates a preferences dialog page
		 *
		 *@param parent the parent widget for the dialog page
		 */
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);  
		// Docu in base class
		virtual void recalculateAxes();
		///
		virtual void invalidate_();
		///	
		virtual void intensityModeChange_();
		///
		virtual Math::Histogram<UnsignedInt, float> createIntensityDistribution_();   
		
		void setMainPreferences(const Param& prefs);
		
	};
	
}//namespace

#endif
