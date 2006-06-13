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

// QT
#include <qimage.h>

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/SpectrumWidget.h>
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
		
		// Docu in SpectrumWidget
		QImage getImage(UnsignedInt width, UnsignedInt height, UnsignedInt flags=0); 
		
		///PreferencesManager
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);
		
		bool isAbsoluteIntensity() const;
		
		inline UnsignedInt getLabelMode() { return label_mode_; }
		
		void minimizeToChart();

		void setMainPreferences(const Param& prefs);

		bool getSnapToMax();
		void setSnapToMax(bool b);
		
	protected:
		virtual void intensityModificationChange_();
		virtual Math::Histogram<UnsignedInt,float> createIntensityDistribution_();
		
		void legendModificationChange_();
		
		/// state variables
		UnsignedInt label_mode_;
		bool log_switched_;
	
	public slots:
		void switchAxis(bool b);
		void setDrawMode(QAction*); //< Sets draw mode to one of the supported types
		void drawModePeaks();
		void drawModeLines();
		void intensityAxisAbsolute();
		void intensityAxisRelative();
		void setZoomFactor(double);  //< Sets zoom_factor_ to non default value.
		void setVisibleArea(double, double);	//< Sets visible area to [position1, position2] and emits visibleAreaChanged
		void mouseMoveEvent( QMouseEvent *e);
		void clearHighlighting(); // Clear canvas highlighting and statusbar message when out of canvas
	
	signals:
		void visibleAreaChanged(double, double); //< Gets emitted whenever the visible area changes.
		
	protected:
		void setLabelMode_(UnsignedInt label_mode);  //<set label mode of widget axis, side effect: set modes of axes as well
		/// Wrappers to retrieve and set label modes in a mapping-safe way. Programmers should use these whenever possible.
		void setIntensityAxisRelative_(); //< sets correct label_mode_ depending on mapping_info_ and previous label_mode_.
		void setIntensityAxisAbsolute_(); //< sets correct label_mode_ depending on mapping_info_ and previous label_mode_.
		bool isIntensityAxisAbsolute_() const;  //< returns true if the intensity axis is absolute, false if relative (percent)

		virtual void recalculateAxes();
	};
} // namespace OpenMS

#endif
