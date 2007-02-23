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
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/VISUAL/MultiGradient.h>

// QT
class QPainter;
class QMouseEvent;
class QWheelEvent;

namespace OpenMS
{
  namespace Internal
  {
    class Spectrum2DCanvasPDP;
  }

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
  	
  	@todo Add reduction of 2D data (Marc)
  	
  	@ingroup spectrum_widgets
  */
  class Spectrum2DCanvas : public SpectrumCanvas
  {
      Q_OBJECT

    public:
      /**	@name Type definitions */
      //@{
      ///Dimensions of the peak/feature data
      enum DimensionId
      {
        MZ = DimensionDescription < LCMS_Tag >::MZ,
        RT = DimensionDescription < LCMS_Tag >::RT
    	};

      ///View modes for 2D dots.
      enum DotModes
      {
        DOT_BLACK = 0,            ///< use black only
        DOT_GRADIENT = 1          ///< use gradient
    	};

      //@}

      /// Default constructor
      Spectrum2DCanvas(QWidget* parent = 0);

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

      /**
      	@brief Sets the mode for 2D dots.

      	Sets the view mode for the dots. Note that this only affects the view
      	if peaks are actually shown as dots.

      	@param mode The new dot mode.
      */
      void setDotMode(SignedInt mode);

      /**
      	@brief Retunrs the mode for 2D dots.

      	Returns the currently set dot mode. If it has not been set yet, it is
      	read from the configuration file.

      	@returns The current dot mode.
      */
      SignedInt getDotMode();

      /**
      	@brief Sets the 2D dot gradient.

      	Sets the color gradient for peaks shown as dots. Peaks are colored
      	according to the gradient and their height, if the dot mode is set to
      	DOT_GRADIENT.

      	@param gradient A string containing the gradient description.
      */
      void setDotGradient(const std::string& gradient);

      /**
      	@brief Set 2D Surface gradient

      	Sets the color gradient used to paint the averaged peak intensities
      	in the background.

      	@param gradient a string representation of the gradient
      */
      void setSurfaceGradient(const std::string& gradient);

      // Docu in base class
      virtual PreferencesDialogPage* createPreferences(QWidget* parent);

      // Docu in base class
      void setMainPreferences(const Param& prefs);

      // Docu in base class
      virtual void repaintAll();

    signals:
      /// Sets the data for the horizontal projection
      void showProjectionHorizontal(const MSExperiment<>&);
      /// Sets the data for the vertical projection
      void showProjectionVertical(const MSExperiment<>&);

    public slots:

      /// Sets if countours are shown
      void showContours(bool on);
      /// Sets if colored surface is shown
      void showSurface(bool on);
      /// Sets if dots are shown
      void showPoints(bool on);

      /// Returns if countours are shown
      bool contoursAreShown();
      /// Returns if colored surface is shown
      bool surfaceIsShown();
      /// Returns if dots are shown
      bool dotsAreShown();

      // Docu in base class
      void activateLayer(int layer_index);
      // Docu in base class
      void removeLayer(int layer_index);
      // Docu in base class
      SignedInt finishAdding(float low_intensity_cutoff = 0);
      // Docu in base class
      virtual void horizontalScrollBarChange(int value);
      // Docu in base class
      virtual void verticalScrollBarChange(int value);

    protected:
      /** @name Reimplemented QT events */
      //@{
      void mousePressEvent(QMouseEvent* e);
      void mouseReleaseEvent(QMouseEvent* e);
      void mouseMoveEvent(QMouseEvent* e);
      void wheelEvent(QWheelEvent* e);
			void mouseDoubleClickEvent(QMouseEvent* e);
			void paintEvent(QPaintEvent* e);
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
      void paintDots_(UnsignedInt layer_index, QPainter& p);

      /**
      	@brief Paints data as a height map.

      	Paints the peak data as interpolated contour lines.
      	The data is shown as a height map such that higher
      	areas are enclosed by more lines than lower areas.
      	
      	@param layer_index The index of the layer.
      	@param p The QPainter to paint on.
      */
      void paintContours_(UnsignedInt layer_index, QPainter& p);

      /**
      	@brief Paints data as a colored surface gradient.

      	Paints the peak data as an interpolated surface gradient.
      	The data is shown according to the gradien which can be
      	set with the setSurfaceGradient() member function.
      	
      	@param layer_index The index of the layer.
      	@param p The QPainter to paint on.
      */
      void paintSurface_(UnsignedInt layer_index, QPainter& p);

      /**
      	@brief Paints all feature convex hulls for a feature layer.
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
      void paintConvexHulls_(UnsignedInt layer_index, QPainter& p);

      /**
      	@brief Paints feature convex hulls.
      	
      	@param hulls Reference to convex hull vector.
      	@param p The QPainter to paint on.
      */
      void paintConvexHulls_(const DFeature<2>::ConvexHullVector& hulls, QPainter& p);

      /**
      	@brief Paints feature pair connections.
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
			void paintFeaturePairConnections_(UnsignedInt layer_index, QPainter& p);
			
      // Docu in base class
      virtual void intensityModeChange_();
      // Docu in base class
      virtual void intensityDistributionChange_();
      // DOcu in base class
      virtual void recalculateSnapFactor_();
      /// recalculates the surface gradient inerpolation values. Use after Intensites or gradient changed
      void recalculateSurfaceGradient_();
      /// recalculates the dot gradient inerpolation values. Use after Intensites or gradient changed
      void recalculateDotGradient_();

      /**
      	@brief Updates the projection data and emits some related signals.
      	
      	Emitted signals are showProjectionHorizontal(const MSExperiment<>&) and 
      	showProjectionVertical(const MSExperiment<>&).
      	
      	@see projection_mz_
      	@see projection_rt_
      */
      void createProjections_(const AreaType& area, bool shift_pressed, bool ctrl_pressed);

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

      /// Performs the marching squares calculations for a layer and stores the matrix in marching_squares_matrices_
      void calculateMarchingSquareMatrix_(UnsignedInt layer_index);
      /// Returns the marching square cell with the smallest data coordinates
      AreaType getOriginCell_(UnsignedInt layer_index);

      /// Highlights a single peak
      void highlightPeak_(QPainter& p, DFeature<2>* peak);

      /// Returns the nearest peak to position @p pos
      DFeature<2>* findNearestPeak_(const QPoint& pos);

      /// marching squares matrices for the layers
      std::vector< std::vector< std::vector<float> > > marching_squares_matrices_;
      /// Contains the highes value in the marching squares matrix foreach layer
      std::vector<float> max_values_;

      /// Flags whether or not to show the height map for each layer
      std::vector<bool> show_contours_;
      /// Flags whether or not to show the surface gradient for each layer
      std::vector<bool> show_surface_;
      /// Flags whether or not to show individual peaks for each layer
      std::vector<bool> show_dots_;

      /// the nearest peak/feature to the mouse cursor (DFeature to be able to store the convex hull too)
      DFeature<2>* selected_peak_;
      /// start peak/feature of measuring mode
      DFeature<2>* measurement_start_;
      /// end peak/feature of measuring mode
      DFeature<2>* measurement_stop_;
      /// temporary peak/feature for findNearestPeak_
      DFeature<2> tmp_peak_;

      /// Gradient for dots
      MultiGradient dot_gradient_;
      /// Gradient for surface
      MultiGradient surface_gradient_;

  };
}

#endif
