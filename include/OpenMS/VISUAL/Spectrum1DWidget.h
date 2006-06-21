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
// $Id: Spectrum1DWidget.h,v 1.21 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM1DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM1DWIDGET_H

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/config.h>

class QAction;

namespace OpenMS
{
	class Spectrum1DCanvas;

	/**
		@brief Widget for visualization of spectrum
		
		@ingroup spectrum_widgets
	*/
	class Spectrum1DWidget : public SpectrumWidget
	{
		Q_OBJECT
		
	public:
		/**
			@brief Default constructor
			
			@param parent The parent widget.
			@param name The widget name.
			@param f Qt::WidgetFlags that are passed on.
		*/
		Spectrum1DWidget(QWidget* parent = 0, const char* name = "Spectrum1DWidget", WFlags f = 0);
		
		///Destructor
		virtual ~Spectrum1DWidget();
	
		Spectrum1DCanvas* canvas() const;
		
		///PreferencesManager
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);
		
		inline Spectrum1DCanvas::LabelMode getLabelMode() 
		{ 
			return label_mode_; 
		}
		
		void minimizeToChart();

		void setMainPreferences(const Param& prefs);

		bool getSnapToMax();
		void setSnapToMax(bool b);

	public slots:
		void switchAxis(bool b);
		void setDrawMode(QAction*); //< Sets draw mode to one of the supported types
		void drawModePeaks();
		void drawModeLines();
		void intensityAxisAbsolute();
		void intensityAxisRelative();
		void setVisibleArea(double, double);	//< Sets visible area to [position1, position2] and emits visibleAreaChanged
		void mouseMoveEvent( QMouseEvent *e);
	
	signals:
		void visibleAreaChanged(double, double); //< Gets emitted whenever the visible area changes.		

	protected:
		// Docu in base class
		virtual void intensityModificationChange_();
		// Docu in base class
		virtual Math::Histogram<UnsignedInt,float> createIntensityDistribution_();
		
		// Docu in base class
		void legendModificationChange_();
		
		/// Label mode: percentage or absolut
		Spectrum1DCanvas::LabelMode label_mode_;
		
		/// Wrappers to retrieve and set label modes in a mapping-safe way. Programmers should use these whenever possible.
		void setIntensityAxisRelative_(); //< sets correct label_mode_ depending on mapping_info_ and previous label_mode_.
		void setIntensityAxisAbsolute_(); //< sets correct label_mode_ depending on mapping_info_ and previous label_mode_.
		bool isIntensityAxisAbsolute_() const;  //< returns true if the intensity axis is absolute, false if relative (percent)
		
		// Docu in base class
		virtual void recalculateAxes();
	};
} // namespace OpenMS

#endif
