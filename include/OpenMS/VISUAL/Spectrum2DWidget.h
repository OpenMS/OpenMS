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

#ifndef OPENMS_VISUAL_SPECTRUM2DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM2DWIDGET_H

// OpenMS
#include <OpenMS/VISUAL/SpectrumWidget.h>

namespace OpenMS
{
	class Spectrum2DCanvas;
	
	/**
		@brief Widget for 2D-visualization of map data
		
		This class is the "model" and "controller" part of
		the two-dimensional spectrum view. The "view" part
		is Spectrum2DCanvas. This widget has the view and
		the two axis-widgets as children and controls them.
		It also provides an interface to these widgets.
		
		If you want to use the 2D spectrum widget, you need
		to use this class. The view features several view
		modes:
		
		- Dots: display peaks as small filled circles.
		- Contour lines: show an interpolated height map
		  by grouping peaks together.
		- Color map: show an interpolated height map as
		  a colored gradient background.
		
		The user can zoom, translate and select peaks. A
		zoom stack is provided for going back to an earlier
		view.
		
		@ingroup spectrum_widgets
	*/
	class Spectrum2DWidget : public SpectrumWidget
	{
		Q_OBJECT
	public:

		/// Default constructor
		Spectrum2DWidget(QWidget* parent = 0);
		/// Destructor
		virtual ~Spectrum2DWidget();
		
		// Docu in base class
		Spectrum2DCanvas* canvas();
		
		// Docu in base class
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);

	public slots:	
		// Docu in base class
		virtual void recalculateAxes_();
		
	signals:
		/**
			@brief Signal emitted whenever the visible area changes.
			
			@param area The new visible area.
		*/
		void visibleAreaChanged(DRange<2> area);

	protected:
		// Docu in base class
		virtual Math::Histogram<UnsignedInt,float> createIntensityDistribution_();
	};
}

#endif
