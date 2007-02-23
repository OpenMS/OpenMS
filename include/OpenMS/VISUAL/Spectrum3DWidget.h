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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM3DWIDGET_H

#include <OpenMS/VISUAL/SpectrumWidget.h>

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
			///Constructor
			Spectrum3DWidget(QWidget* parent = 0);
			
			/// Destructor
			virtual ~Spectrum3DWidget();
			
			// Docu in base class
			Spectrum3DCanvas* canvas();
			
			// Docu in base class
			virtual PreferencesDialogPage* createPreferences(QWidget* parent);  
			// Docu in base class
			virtual void recalculateAxes_();
			// Docu in base class
			virtual Math::Histogram<UnsignedInt, float> createIntensityDistribution_();   
			// Docu in base class
			virtual QImage getImage(UnsignedInt width, UnsignedInt height);
	
			bool isLegendShown() ;
		
		public slots:
			void showLegend(int);

		protected:
			bool legend_shown_;

	};
	
}//namespace

#endif
