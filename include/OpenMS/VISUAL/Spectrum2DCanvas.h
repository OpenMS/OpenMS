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


#ifndef OPENMS_VISUAL_SPECTRUM2DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM2DCANVAS_H

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/MultiGradient.h>

// QT
class QPainter;
class QMouseEvent;
class QWheelEvent;

namespace OpenMS
{
  /**
  	@brief Canvas for 2D-visualization of map data

  	This widget displays a 2D representation of a set
  	of peaks. There are 3 independent view modes:

  	- Dots: display peaks as small filled circles.
  	- Contour lines: show an interpolated height map
  	  by grouping peaks together.
  	- Color map: show an interpolated height map as
  	  a colored gradient background.

  	The user can zoom, translate and select peaks. A
  	zoom stack is provided for going back to an earlier
  	view.
		
		@ref Spectrum2DCanvas_Parameters are explained on a separate page.
		
  	@ingroup SpectrumWidgets
  */
  class Spectrum2DCanvas 
  	: public SpectrumCanvas
  {
      Q_OBJECT

    public:
      /// Default constructor
      Spectrum2DCanvas(const Param& preferences, QWidget* parent = 0);

      /// Destructor
      ~Spectrum2DCanvas();

      /**
      	@brief Draws the contents.

      	Device independent drawing function. Draws the contents on painter @p p.
      	This function follows the 

      	@param p The QPainter to draw the data to.
      	@param width The width of the canvas in pixels.
      	@param height The height of the canvas in pixels.
      */
      void print(QPainter& p, int width, int height);

			// Docu in base class
			virtual void showCurrentLayerPreferences();

    signals:
      /// Sets the data for the horizontal projection
      void showProjectionHorizontal(const MSExperiment<>&, Spectrum1DCanvas::DrawModes);
      /// Sets the data for the vertical projection
      void showProjectionVertical(const MSExperiment<>&, Spectrum1DCanvas::DrawModes);
      /// Shows the number of peaks and the intensity sum of the projection
      void showProjectionInfo(int, double, double);
      /// Requests to display the current peak data in 3D 
      void showCurrentPeaksAs3D();
			/// Requests to display the spectrum with index @p index in 1D
			void showSpectrumAs1D(int index);
		
    public slots:

      // Docu in base class
      void activateLayer(int layer_index);
      // Docu in base class
      void removeLayer(int layer_index);
      // Docu in base class
      Int finishAdding(float low_intensity_cutoff = 0);
      // Docu in base class
      virtual void horizontalScrollBarChange(int value);
      // Docu in base class
      virtual void verticalScrollBarChange(int value);
      /**
      	@brief Updates the projection data and emits some related signals.
      	
      	Emitted signals are showProjectionHorizontal(const MSExperiment<>&, Spectrum1DCanvas::DrawModes) and 
      	showProjectionVertical(const MSExperiment<>&, Spectrum1DCanvas::DrawModes).
      	
      	@see projection_mz_
      	@see projection_rt_
      */
      void showProjections();
      
    protected:
      /** @name Reimplemented QT events */
      //@{
      void mousePressEvent(QMouseEvent* e);
      void mouseReleaseEvent(QMouseEvent* e);
      void mouseMoveEvent(QMouseEvent* e);
      void wheelEvent(QWheelEvent* e);
			void mouseDoubleClickEvent(QMouseEvent* e);
			void paintEvent(QPaintEvent* e);
			void contextMenuEvent(QContextMenuEvent* e);
      //@}

      // Docu in base class
      virtual void updateScrollbars_();

      /**
      	@brief Paints individual peaks.

      	Paints the peaks as small ellipses. The peaks are colored according to the
      	selected dot gradient.
      	
      	@param layer_index The index of the layer.
      	@param p The QPainter to paint on.
      */
      void paintDots_(UInt layer_index, QPainter& p);

      /**
      	@brief Paints convex hulls (one for each mass trace) for a features of a layer.
      	
      	@param layer_index Int of the layer.
      	@param p The QPainter to paint on.
      */
      void paintTraceConvexHulls_(UInt layer_index, QPainter& p);

      /**
      	@brief Paints the convex hulls (one for each feature) for a features of a layer.
      	
      	@param layer_index Int of the layer.
      	@param p The QPainter to paint on.
      */
      void paintFeatureConvexHulls_(UInt layer_index, QPainter& p);

      /**
      	@brief Paints convex hulls (one for each mass trace) for a single feature.
      	
      	@param hulls Reference to convex hull vector.
      	@param p The QPainter to paint on.
      */
      void paintConvexHulls_(const std::vector<ConvexHull2D>& hulls, QPainter& p);

      /**
      	@brief Paints feature pair connections.
      	
      	@param layer_index Int of the layer.
      	@param p The QPainter to paint on.
      */
			void paintFeaturePairConnections_(UInt layer_index, QPainter& p);
			
      // Docu in base class
      virtual void intensityModeChange_();
      // Docu in base class
      virtual void intensityDistributionChange_();
      // DOcu in base class
      virtual void recalculateSnapFactor_();
			// Docu in base class
			virtual void currentLayerParamtersChanged_();
      /// recalculates the dot gradient of a layer
      void recalculateDotGradient_(UInt layer);

      /// m/z projection data
      MSExperiment<> projection_mz_;
      /// RT projection data
      MSExperiment<> projection_rt_;


      /// interpolation helper function
      float betweenFactor_(float v1, float v2, float val);
      /**
      	@brief Returns the color associated with @p val for the gradient @p gradient.
      	
      	Takes intensity modes into account.
      */
      const QColor& heightColor_(float val, const MultiGradient& gradient);

      /// Highlights a single peak
      void highlightPeak_(QPainter& p, const Feature* peak);

      /// Returns the nearest peak to position @p pos
      const Feature* findNearestPeak_(const QPoint& pos);

      /// the nearest peak/feature to the mouse cursor (DFeature to be able to store the convex hull too)
      const Feature* selected_peak_;
      /// start peak/feature of measuring mode
      const Feature* measurement_start_;
      /// end peak/feature of measuring mode
      const Feature* measurement_stop_;
      /// temporary peak/feature for findNearestPeak_
      Feature tmp_peak_;

  };
}

#endif
