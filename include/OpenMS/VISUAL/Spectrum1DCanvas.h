// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: Spectrum1DCanvas.h,v 1.36 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM1DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM1DCANVAS_H

#include <OpenMS/config.h>

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/RubberBand.h>

//QT
class QAction;

namespace OpenMS
{

	namespace Internal
	{
		class Spectrum1DCanvasPDP;
	}
	
	/**
		@brief Canvas for visualization of spectrum
		
		
		
		@ingroup spectrum_widgets
	*/
	class Spectrum1DCanvas : public SpectrumCanvas
	{
		Q_OBJECT
		
		friend class Internal::Spectrum1DCanvasPDP;
		friend class Spectrum1DWidget;
		
	public:
		/// Icons for marking peaks
		enum Icons
		{
			IT_NOICON,
			IT_CIRCLE,
			IT_TRIANGLE,
			IT_ASTERIX, 
			IT_SQUARE
		};
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
		virtual ~Spectrum1DCanvas();
	
		///Enumerate all avaiable paint styles
		enum DrawModes 
		{
			DM_PEAKS,						//< draw data as peak
			DM_CONNECTEDLINES		//< draw as connected lines
		};
	
		// some high level functions that use low level widget functions
		// zoom and translate operations (with animation, visualisation)
		void zoomOut(double x);
		void zoomIn(double x);
		void translate(double x, int animatedSteps=1);

		// function to mark a peak with an icon
		inline void setPeakIcon(unsigned int index, unsigned int icon) 
		{ 
			currentDataSet()[0].getContainer()[index].setMetaValue(UnsignedInt(4),SignedInt(icon)); 
		}

		/**
			@brief Returns selected peaks.
			
			The actual selected peaks are framed by the first and last peak
		*/
		std::vector<SpectrumIteratorType> getSelectedPeaks();
		
		inline int getDrawMode() const
		{ 
			return draw_modes_[current_data_]; 
		}
		
		// Docu in base class
		virtual void setMainPreferences(const Param& prefs);
		
		///PreferencesManager
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);
		
		/// Returns the snap_to_max_mode_ flag
		bool getSnapToMax();
		/// Sets the snap_to_max_mode_ flag
		void setSnapToMax(bool b);
		
		/**
			@brief returns the snap_factor_.
			
			@see snap_factor_
		*/
		double getSnapFactor();

		bool isAbsoluteIntensity() const;	

	public slots:
		
		void setDrawMode(QAction*); //< Sets draw mode to one of the supported types
		void drawModePeaks();
		void drawModeLines();
		void intensityAxisAbsolute();
		void intensityAxisRelative();
		
		// Docu in base class
		void activateDataSet(int data_set);
		// Docu in base class
		void removeDataSet(int data_set);
		// Docu in base class
		SignedInt finishAdding();

		void setVisibleArea(double lo, double hi);
		void setVisibleArea(DRange<2> range); //Do not change this to AreaType the signal needs QT needs the exact type...
	
	protected:
		
		void drawIcon(const PeakType& peak, const QPoint& p);
		
		/**
			@brief Changes visible area interval
			
			This method is for convenience only. It calls changeVisibleArea_(const AreaType&) .
		*/
		void changeVisibleArea_(double lo, double hi);  
		
		/// Calls chartToWidget_(const PointType&) but takes snap_factor_ and layer_factor_ into account.
		QPoint chartToWidget_(const PeakType& peak);
		
		/**
			@brief Updates visible_begin_ and visible_end_ .
			
			In snap-to-max-intensity mode it also updates the maximum intensity in the visible area.
		*/
		void updateVisibleAreaBounds_();
		
		// Docu in base class
		virtual void intensityModificationChange_();
		
		/// RubberBand for zooming
		RubberBand rubber_band_;

		// Docu in base class
		virtual void invalidate_();
		
		// Docu in base class
		void changeVisibleArea_(const AreaType& new_area);

		/// Array of selected peak iterators
		std::vector<SpectrumIteratorType> selected_peaks_;
		
		/// The (one and only) painter. 
		QPainter painter_;	

		/// Flag that indicates if intensity is absolute or relative
		bool absolute_intensity_;
		/// Scaling factor for relative scale with multiple layers
		double layer_factor_;
		/// Flag for 'snap to maximum intensity mode'.
		bool snap_to_max_mode_;
		/// Itensity multiplication factor for 'snap to maximum intensity mode'.
		double snap_factor_;

		/// start position of mouse actions
		QPoint action_start_pos_;
		/// current position of mouse actions
		QPoint action_current_pos_;

		/// Draw modes (for each spectrum)
		std::vector<DrawModes> draw_modes_;
		/// Iterators on first visible peak (for each spectrum)
		std::vector<SpectrumIteratorType> visible_begin_;
		/// Iterators after the last visible peak (for each spectrum)
		std::vector<SpectrumIteratorType> visible_end_;    
		/// Iterator on peak next to mouse position
		SpectrumIteratorType nearest_peak_;
		/// Find peak next to the given position
		SpectrumIteratorType findPeakAtPosition(QPoint);  
		
		/// pen for drawing of normal peaks
		QPen norm_pen_;
		/// pen for drawing of highlighted peaks
		QPen high_pen_;
		/// pen for drawing of icons
		QPen icon_pen_; 

		/// Draws peaks for dataset @p index
		void drawPeaks_(UnsignedInt index);
		/// Draws connectedLines for dataset @p index
		void drawConnectedLines_(UnsignedInt index);

		/// QT Event
		void viewportPaintEvent( QPaintEvent * );
		/// QT Event
		void contentsMousePressEvent( QMouseEvent *);
		/// QT Event
		void contentsMouseDoubleClickEvent( QMouseEvent *);
		/// QT Event
		void contentsMouseReleaseEvent( QMouseEvent *);
		/// QT Event
		void contentsMouseMoveEvent( QMouseEvent *);


	};
} // namespace OpenMS

#endif
