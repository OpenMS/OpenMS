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
		
		//icons
		static const UnsignedInt IT_NOICON = 0x00;
		static const UnsignedInt IT_CIRCLE = 0x01;
		static const UnsignedInt IT_TRIANGLE = 0x02;
		static const UnsignedInt IT_ASTERIX = 0x03;
		static const UnsignedInt IT_SQUARE = 0x04;

		// label modes of x axis and y axis
		static const UnsignedInt LM_XABSOLUTE_YABSOLUTE = 0x00;
		static const UnsignedInt LM_XPERCENT_YABSOLUTE = 0x01;
		static const UnsignedInt LM_XABSOLUTE_YPERCENT = 0x02;
		static const UnsignedInt LM_XPERCENT_YPERCENT = 0x03;
		
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
		void zoomOut(double x,int animatedSteps=1);
		void zoomIn(double x, int animatedSteps=1);
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
		
		inline double getZoomFactor() const 
		{ 
			return zoom_factor_; 
		}
		
		// Docu in base class
		virtual void setMainPreferences(const Param& prefs);
		
		///PreferencesManager
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);

		bool getSnapToMax();
		void setSnapToMax(bool b);
		/**
			@brief returns the snap_factor_.
			
			@see snap_factor_
		*/
		double getSnapFactor();

		bool isAbsoluteIntensity() const;	

		void clearHighlighting();

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

		void setZoomFactor(double);  //< Sets zoom_factor_ to non default value.
		void setVisibleArea(double lo, double hi);
		void setVisibleArea(DRange<2> range); //Do not change this to AreaType the signal needs QT needs the exact type...
	
	protected:
		// Docu in base class
		virtual const AreaType& getDataArea_();
		
		void drawIcon(const PeakType& peak, const QPoint& p);
		
		/**
			@brief Changes visible area interval
			
			This method is for convenience only. It calls changeVisibleArea_(const AreaType&) .
		*/
		void changeVisibleArea_(double lo, double hi);  
		
		/// Calls chartToWidget_(const PointType&) but takes snap_factor_ and layer_factor_ into account.
		QPoint chartToWidget_(const PeakType& peak);
		
		void setBounds_();
		
		// Docu in base class
		QPoint chartToWidget_(const PointType& pos);
		
		// Docu in base class
		virtual void intensityModificationChange_();

		void legendModificationChange_();

		/// scale data depending on log/linear and absolute/percent
		void scaleData_(bool is_log);
		void scaleAllData_(bool is_log);

		RubberBand rubber_band_;

		// Docu in base class
		virtual void invalidate_();
		
		// Docu in base class
		void changeVisibleArea_(const AreaType& new_area);

		/// Structure that (one day) will hold all style information (like pen color, font name, ...).
		struct Style
		{
			double border;
		} style_ ;

		/**
			*	Helper variables (get set in adjustBuffer_) that simplify layout.
			*/
		double grid_row_width_, grid_collumn_height_;

		std::vector<SpectrumIteratorType> selected_peaks_;

		QPainter painter_;	//< the (one and only) painter. 
		
		// state variables
		double zoom_factor_;
		std::vector<int> draw_modes_;
		double margin_;
		
		///for storing the old maximum when the intensities are transformed
		double old_max_intensity_;
		bool absolute_intensity_;

		//scaling factor for relative scale with multiple layers
		double layer_factor_;
		//layer with highest maximum
		UnsignedInt max_layer_;

		bool snap_to_max_mode_;
		
		///Itensity multiplication factor for 'snap to maximum intensity mode'.
		double snap_factor_;

		// selected diagram action mode and helper variables
		QPoint action_start_pos_;
		QPoint action_current_pos_;

		/// Area that encloses all peaks of all datasets
		AreaType data_area_;
		
		std::vector<SpectrumIteratorType> visible_begin_;  //< iterator on first visible peak
		std::vector<SpectrumIteratorType> visible_end_;    //< iterator on one after the last visible peak
		SpectrumIteratorType nearest_peak_; //< iterator on peak next to mouse position
		SpectrumIteratorType findPeakAtPosition(QPoint);  //< find peak next to position
		
		// data set drawing
		QPen norm_pen_; // pen for drawing of normal peaks
		QPen high_pen_; // pen for drawing of highlighted peaks
		QPen icon_pen_; // pen for drawing of icons

		QString zoom_status_;  // Status message in zoom-mode
		bool is_highlighted_;   // peak is highlighted

		void drawPoints_(UnsignedInt index);
		void drawPeaks_(UnsignedInt index);
		void drawConnectedLines_(UnsignedInt index);

		/// EVENTS
		/// Most of the drawing is done during paint event
		void viewportPaintEvent( QPaintEvent * );
		/// Mouse events.
		void contentsMousePressEvent( QMouseEvent *);
		void contentsMouseDoubleClickEvent( QMouseEvent *);
		void contentsMouseReleaseEvent( QMouseEvent *);
		void contentsMouseMoveEvent( QMouseEvent *);


	};
} // namespace OpenMS

#endif
