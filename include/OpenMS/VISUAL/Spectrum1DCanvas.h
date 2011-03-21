// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM1DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM1DCANVAS_H

// STL
#include <vector>
#include <utility>

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>

//QT
class QAction;

namespace OpenMS
{
	/**
		@brief Canvas for visualization of one or several spectra.
		
		@image html Spectrum1DCanvas.png
		
		The example image shows %Spectrum1DCanvas displaying a raw data layer and a peak data layer. 
		
		@todo Use spectrum StringDataArray with name 'label' for peak annotations (Hiwi, Johannes)
		
		@htmlinclude OpenMS_Spectrum1DCanvas.parameters
				
		@ingroup SpectrumWidgets
	*/
	class OPENMS_GUI_DLLAPI Spectrum1DCanvas 
		: public SpectrumCanvas
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
			
			/// Default constructor
			Spectrum1DCanvas(const Param& preferences, QWidget* parent = 0);
			/// Destructor
			virtual ~Spectrum1DCanvas();
		
			///Enumerate all avaiable paint styles
			enum DrawModes 
			{
				DM_PEAKS,						//< draw data as peak
				DM_CONNECTEDLINES		//< draw as connected lines
			};

			/// Returns the draw mode of the current layer
			DrawModes getDrawMode() const;
	
			/// Sets draw mode of the current layer
			void setDrawMode(DrawModes mode);
			
			// Docu in base class
			virtual void showCurrentLayerPreferences();

			// Docu in base class
			virtual void saveCurrentLayer(bool visible);
	
			/// Returns whether flipped layers exist or not
			bool flippedLayersExist();
			
			/// Flips the layer with @p index up/downwards
			void flipLayer(Size index);
			
			/// Returns whether this widget is currently in mirror mode
			bool mirrorModeActive();
			
			/// Sets whether this widget is currently in mirror mode
			void setMirrorModeActive(bool b);
		
			/// For convenience - calls dataToWidget(float, float, ...)
			void dataToWidget(const PeakType& peak, QPoint& point, bool flipped = false, bool percentage = true);
			
			/// Calls SpectrumCanvas::dataToWidget_(), takes mirror mode into account
			void dataToWidget(float x, float y, QPoint& point, bool flipped = false, bool percentage = false);
			
			/// For convenience - calls widgetToData(float, float, ...)
			PointType widgetToData(const QPoint& pos, bool percentage = false);
			
			/// Calls SpectrumCanvas::widgetToData_(), takes mirror mode into account
			PointType widgetToData(float x, float y, bool percentage = false);
			
      /// Draws all annotation items of @p layer_index on @p painter
      void drawAnnotations(Size layer_index, QPainter& painter);
			
			/// Performs an alignment of the layers with @p layer_index_1 and @p layer_index_2
			void performAlignment(Size layer_index_1, Size layer_index_2, const Param& param);
			
			/// Resets alignment_
			void resetAlignment();
			
			/// Draws the alignment on @p painter
			void drawAlignment(QPainter& painter);
			
			/// Returns the number of aligned pairs of peaks
			Size getAlignmentSize();
			
			/// Returns the score of the alignment
			DoubleReal getAlignmentScore();
			
			/// Sets current spectrum index of current layer to @p index
			void activateSpectrum(Size index, bool repaint=true);

      /// is the widget shown vertically? (for projections)
      void setSwappedAxis(bool swapped);

      /// Set's the Qt PenStyle of the active layer
      void setCurrentLayerPeakPenStyle(Qt::PenStyle ps);

      /// Actual painting takes place here
      void paint(QPainter* paint_device, QPaintEvent* e);
  signals:
      /// Requests to display all spectra in 2D plot
      void showCurrentPeaksAs2D();
      /// Requests to display all spectra in 3D plot
      void showCurrentPeaksAs3D();

		public slots:
			// Docu in base class
			void activateLayer(Size layer_index);
			// Docu in base class
			void removeLayer(Size layer_index);
      //docu in base class
      virtual void updateLayer(Size i);

			/**
				@brief Sets the visible area.
				
				Sets the visible area to a new value. Note that it does not emit visibleAreaChanged()
				@param range the new visible area
			*/
			void setVisibleArea(DRange<2> range); //Do not change this to AreaType the signal needs QT needs the exact type...
			// Docu in base class
			virtual void horizontalScrollBarChange(int value);

  	protected slots:
  		
			/// Reacts on changed layer paramters
			void currentLayerParamtersChanged_();

		protected:
			// Docu in base class
			bool finishAdding_();

      /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
	    void drawCoordinates_(QPainter& painter, const PeakIndex& peak);
	    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
	    void drawDeltas_(QPainter& painter, const PeakIndex& start, const PeakIndex& end);
			
			/**
				@brief Changes visible area interval
				
				This method is for convenience only. It calls changeVisibleArea_(const AreaType&, bool, bool) .
			*/
			void changeVisibleArea_(double lo, double hi, bool repaint = true, bool add_to_stack = false);  
			
			/// Draws a highlighted peak; if draw_elongation is true, the elongation line is drawn (for measuring)
			void drawHighlightedPeak_(Size layer_index, const PeakIndex& peak, QPainter& painter, bool draw_elongation = false);
			
			/// Draws a dashed line using the highlighted peak color parameter
			void drawDashedLine_(const QPoint& from, const QPoint& to, QPainter& painter);
			
      /// Recalculates the current scale factor based on the specified layer (= 1.0 if intensity mode != IM_PERCENTAGE)
      void updatePercentageFactor_(Size layer_index);

			/**
				@brief Sets the visible area
				
				Changes the visible area, adjustes the zoom stack and notifies interested clients about the change. 
				If parts of the area are outside of the data area, the new area will be adjusted.
				
				@param new_area The new visible area.
				@param repaint if repainting of the widget should ne triggered
				@param add_to_stack If the new area is to add to the zoom_stack_
			*/
			virtual void changeVisibleArea_(const AreaType& new_area, bool repaint = true, bool add_to_stack = false);
			// Docu in base class
			virtual void recalculateSnapFactor_();
			// Docu in base class
			virtual void updateScrollbars_();
			// Docu in base class
			virtual void intensityModeChange_();

			/// Draw modes (for each spectrum)
			std::vector<DrawModes> draw_modes_; 
      /// Draw style
      std::vector<Qt::PenStyle> peak_penstyle_;

      /// start point of "ruler" for measure mode
      QPoint measurement_start_point_;
      /// Indicates whether this widget is currently in mirror mode
			bool mirror_mode_;

			/// Indicates whether annotation items are just being moved on the canvas
			bool moving_annotations_;

      /// Indicates whether an alignment is currently visualized
      bool show_alignment_;
      /// Layer index of the first alignment layer
      Size alignment_layer_1_;
      /// Layer index of the second alignment layer
      Size alignment_layer_2_;
      /// Stores the alignment as MZ values of pairs of aligned peaks in both spectra
      std::vector<std::pair<DoubleReal, DoubleReal > > aligned_peaks_mz_delta_;
      /// Stores the peak indizes of pairs of aligned peaks in both spectra
      std::vector<std::pair<Size, Size> > aligned_peaks_indices_;

      /// Stores the score of the last alignment
			DoubleReal alignment_score_;
      /// is this widget showing data with swapped m/z and RT axis? (for drawCoordinates_ only)
      bool is_swapped_;

			/// Find peak next to the given position
			PeakIndex findPeakAtPosition_(QPoint);

      /// Shows dialog and calls addLabelAnnotation_
      void addUserLabelAnnotation_(const QPoint& screen_position);
      /// Adds an annotation item at the given screen position
      void addLabelAnnotation_(const QPoint& screen_position, QString label_text);
      /// Shows dialog and calls addPeakAnnotation_
      void addUserPeakAnnotation_(PeakIndex near_peak);
      /// Add an annotation item for the given peak
      void addPeakAnnotation_(PeakIndex peak_index, QString text);

			/// Ensure that all annotations are within data range
			void ensureAnnotationsWithinDataRange_();
	
	    /** @name Reimplemented QT events */
	    //@{
			void paintEvent(QPaintEvent* e);
			void mousePressEvent(QMouseEvent* e);
			void mouseReleaseEvent(QMouseEvent* e);
			void mouseMoveEvent(QMouseEvent* e);
			void keyPressEvent(QKeyEvent* e);
			void contextMenuEvent(QContextMenuEvent* e);
	    //@}
			
      ///Go forward in zoom history
	    virtual void zoomForward_();
      //docu in base class
			virtual void translateLeft_();
			//docu in base class
			virtual void translateRight_();
			//docu in base class
			virtual void paintGridLines_(QPainter& painter);
	};
} // namespace OpenMS

#endif
