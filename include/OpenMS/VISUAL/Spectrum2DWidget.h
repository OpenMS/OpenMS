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
// $Id: Spectrum2DWidget.h,v 1.24 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
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
		/**	@name Type definitions */
		//@{
		///View modes for 2D dots.
		enum DotModes 
		{
			DOT_BLACK = 0,            ///< use black only
			DOT_GRADIENT = 1          ///< use gradient
		};
		//@}

		/**
			@brief Constructor
			
			Spectrum2DWidget constructor. See QWidget for details.
			
			@param parent The parent widget.
			@param name The widget's name.
			@param f Widget flags.
		*/
		Spectrum2DWidget(QWidget* parent = 0, const char* name = "Spectrum2DWidget", WFlags f = 0);
		
		/**
			@brief Destructor
			
			Destroys the Widget and all associated data.
		*/
		virtual ~Spectrum2DWidget();
		
		/**
			@brief Returns the canvas widget
			
			Returns the canvas widget which provides the data view.
			
			@return the canvas widget
		*/
		Spectrum2DCanvas* canvas() const;
		
		/**
			@brief Creates a preferences dialog page
			
			Creates a preferences dialog page for configuring this view.
			
			@param parent the parent widget for the dialog page
		*/
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);

	public slots:	
		// Docu in base class
		virtual void recalculateAxes();
		
	signals:
		/**
			@brief Signal emitted whenever the visible area changes.
			
			@param area The new visible area.
		*/
		void visibleAreaChanged(DRange<2> area);

	protected:
		virtual void intensityModeChange_();
		// Docu in base class
		virtual Math::Histogram<UnsignedInt,float> createIntensityDistribution_();
		
	private:
		/// shows the context menu at position p
		void showContextMenu_(QPoint p);
	};
}

#endif
