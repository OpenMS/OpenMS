// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Publicf
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

#ifndef OPENMS_VISUAL_SPECTRUM1DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM1DCANVAS_H

#include <OpenMS/config.h>

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
// #include <OpenMS/VISUAL/RubberBand.h>

//QT
#include <qpainter.h>

class QAction;

namespace OpenMS
{

	namespace Internal
	{
		class Spectrum1DCanvasPDP;
	}
	
	/**
		@brief Canvas for visualization of spectrum
		
		@todo Add measurement (Marc)
		
		@bug Logarithmic intensity mode is broken
		
		@ingroup spectrum_widgets
	*/
	class Spectrum1DCanvas : public SpectrumCanvas
	{
		Q_OBJECT
		
	public:
		/// Label modes (percentage or absolut) of x axis and y axis
		enum LabelMode
		{
			LM_XABSOLUTE_YABSOLUTE,
			LM_XPERCENT_YABSOLUTE,
			LM_XABSOLUTE_YPERCENT,
			LM_XPERCENT_YPERCENT
		};
		
		/**
			@brief Default constructor
			
			@param parent The parent widget.
			@param name The widget name.
			@param f Qt::WidgetFlags that are passed on.
		*/
		Spectrum1DCanvas(QWidget* parent = 0, const char* name = "Spectrum1DCanvas", WFlags f = 0);	
		
		/// Destructor
		virtual ~Spectrum1DCanvas();
	
		///Enumerate all avaiable paint styles
		enum DrawModes 
		{
			DM_PEAKS,						//< draw data as peak
			DM_CONNECTEDLINES		//< draw as connected lines
		};

		/// Function to mark a peak with an icon
		inline void setPeakIcon(unsigned int index, unsigned int icon) 
		{ 
			currentPeakData_()[0].getContainer()[index].setMetaValue(UnsignedInt(4),SignedInt(icon)); 
		}

		/**
			@brief Returns selected peaks.
			
			The actual selected peaks are framed by the first and last peak
		*/
		std::vector<SpectrumIteratorType> getSelectedPeaks();
		
		/// Returns the draw mode of the current layer
		DrawModes getDrawMode() const;

		/// Sets draw mode of the current layer
		void setDrawMode(DrawModes mode);
		
		// Docu in base class
		virtual void setMainPreferences(const Param& prefs);
		
		// Docu in base class
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);
		
	public slots:
		// Docu in base class
		void activateLayer(int layer_index);
		// Docu in base class
		void removeLayer(int layer_index);
		// Docu in base class
		SignedInt finishAdding(float low_intensity_cutoff = 0);
		
		/**
			@brief Sets the visible area.
			
			Sets the visible area to a new value. Note that it does not emit visibleAreaChanged()
			@param range the new visible area
		*/
		void setVisibleArea(DRange<2> range); //Do not change this to AreaType the signal needs QT needs the exact type...
		// Docu in base class
		virtual void horizontalScrollBarChange(int value);
		// Docu in base class
		virtual void verticalScrollBarChange(int value);
	
	protected:
		/// Draws the icon defined in the meta info of @p peak at the position @p p
		void drawIcon_(const PeakType& peak, const QPoint& p);
		
		/**
			@brief Changes visible area interval
			
			This method is for convenience only. It calls changeVisibleArea_(const AreaType&) .
		*/
		void changeVisibleArea_(double lo, double hi, bool add_to_stack = false);  
		
		/// Calls dataToWidget_(const PointType&, QPoint& point) but takes snap_factor_ and percentage_factor_ into account.
		void dataToWidget_(const PeakType& peak, QPoint& point);
		
/*
		/// RubberBand for zooming
		RubberBand rubber_band_;
*/
    /**
    	@brief Reblits the drawing buffer onto the screen.

    	Refreshes the screen. This function should be called
    	when the internal buffer is still current and the screen
    	representation was damaged by overdrawing. If the internal
    	buffer is outdated, invalidate_() should be called.
    */
    void refresh_();

		// Docu in base class
		virtual void invalidate_();
		
		/**
			@brief Sets the visible area
			
			Changes the visible area, adjustes the zoom stack and notifies interested clients about the change. 
			If parts of the area are outside of the data area, the new area will be adjusted.
			
			@param new_area The new visible area.
			@param add_to_stack If the new area is to add to the zoom_stack_
		*/
		virtual void changeVisibleArea_(const AreaType& new_area, bool add_to_stack = false);
		// Docu in base class
		virtual void recalculateSnapFactor_();
		// Docu in base class
		virtual void updateScrollbars_();
		
		/// Array of selected peak iterators
		std::vector<SpectrumIteratorType> selected_peaks_;
		
		/// The (one and only) painter. 
		QPainter painter_;

		/// Draw modes (for each spectrum)
		std::vector<DrawModes> draw_modes_;
		/// Iterators on first visible peak (for each spectrum)
		std::vector<SpectrumIteratorType> visible_begin_;
		/// Iterators after the last visible peak (for each spectrum)
		std::vector<SpectrumIteratorType> visible_end_;    
		/// Iterator on peak next to mouse position
		SpectrumIteratorType nearest_peak_;
		/// Find peak next to the given position
		SpectrumIteratorType findPeakAtPosition_(QPoint);  
		
		/// pen for drawing of normal peaks
		QPen norm_pen_;
		/// pen for drawing of highlighted peaks
		QPen high_pen_;
		/// pen for drawing of icons
		QPen icon_pen_; 

		/// Draws peaks for layer @p index
		void drawPeaks_(UnsignedInt index);
		/// Draws connectedLines for layer @p index
		void drawConnectedLines_(UnsignedInt index);

		/// QT Event
		void paintEvent( QPaintEvent * );
		/// QT Event
		void mousePressEvent( QMouseEvent *);
		/// QT Event
		void mouseDoubleClickEvent( QMouseEvent *);
		/// QT Event
		void mouseReleaseEvent( QMouseEvent *);
		/// QT Event
		void mouseMoveEvent( QMouseEvent *);


	};
} // namespace OpenMS

#endif
