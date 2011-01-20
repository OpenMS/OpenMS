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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_SPECTRUM2DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM2DCANVAS_H

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/KERNEL/PeakIndex.h>

// QT
class QPainter;
class QMouseEvent;

namespace OpenMS
{
    /**
  	@brief Canvas for 2D-visualization of peak map, feature map and consensus map data

  	This widget displays a 2D representation of a set of peaks, features or consensus elements.

        @image html Spectrum2DCanvas.png

        The example image shows %Spectrum2DCanvas displaying a peak layer and a feature layer.

        @htmlinclude OpenMS_Spectrum2DCanvas.parameters

        @improvement Add RT interpolation mode for high zoom in 2D View (Hiwi)

        @improvement Snap also to min intensity (Hiwi)
		
  	@ingroup SpectrumWidgets
  */
    class OPENMS_DLLAPI Spectrum2DCanvas
  	: public SpectrumCanvas
    {
        Q_OBJECT

    public:
        /// Default constructor
        Spectrum2DCanvas(const Param& preferences, QWidget* parent = 0);

        /// Destructor
        ~Spectrum2DCanvas();

        // Docu in base class
        virtual void showCurrentLayerPreferences();

        // Docu in base class
        virtual void saveCurrentLayer(bool visible);

        /// Merges the features in @p map into the features layer @p i
        void mergeIntoLayer(Size i, FeatureMapSharedPtrType map);

        /// Merges the consensus features in @p map into the features layer @p i
        void mergeIntoLayer(Size i, ConsensusMapSharedPtrType map);

        /// Merges the peptide identifications in @p peptides into the peptide layer @p i
        void mergeIntoLayer(Size i, std::vector<PeptideIdentification>& peptides);

    signals:
        /// Sets the data for the horizontal projection
        void showProjectionHorizontal(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes);
        /// Sets the data for the vertical projection
        void showProjectionVertical(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes);
        /// Shows the number of peaks and the intensity sum of the projection
        void showProjectionInfo(int, double, double);
        /// Signal emitted when the projections are to be shown/hidden
        void toggleProjections();
        /// Requests to display the spectrum with index @p index in 1D
        void showSpectrumAs1D(int index);
        /// Requests to display all spectra in 3D plot
        void showCurrentPeaksAs3D();

    public slots:
        // Docu in base class
        void activateLayer(Size layer_index);
        // Docu in base class
        void removeLayer(Size layer_index);
        //docu in base class
        virtual void updateLayer(Size i);
        // Docu in base class
        virtual void horizontalScrollBarChange(int value);
        // Docu in base class
        virtual void verticalScrollBarChange(int value);
        /**
      	@brief Updates the projection data and emits some related signals.
      	
        Emitted signals are showProjectionHorizontal(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes) and
        showProjectionVertical(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes).
      	
      	@see projection_mz_
      	@see projection_rt_
      */
        void updateProjections();

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

        /** @name Reimplemented QT events */
        //@{
        void mousePressEvent(QMouseEvent* e);
        void mouseReleaseEvent(QMouseEvent* e);
        void mouseMoveEvent(QMouseEvent* e);
        void paintEvent(QPaintEvent* e);
        void contextMenuEvent(QContextMenuEvent* e);
        void keyPressEvent(QKeyEvent* e);
        void keyReleaseEvent(QKeyEvent* e);
        void mouseDoubleClickEvent(QMouseEvent* e);
        //@}

        // Docu in base class
        virtual void updateScrollbars_();

        /**
      	@brief Paints individual peaks.

        Calls different painting methods depending on the layer type and the density of displayed peaks
      	
      	@param layer_index The index of the layer.
      	@param p The QPainter to paint on.
      */    
        void paintDots_(Size layer_index, QPainter& p);

        void paintAllIntensities_(Size layer_index, DoubleReal average_spacing_mz, DoubleReal average_spacing_rt, QPainter& painter);
        /**
        @brief Paints maximum intensity of individual peaks.

        Paints the peaks as small ellipses. The peaks are colored according to the
        selected dot gradient.

        @param layer_index The index of the layer.
        @param p The QPainter to paint on.
      */
        void paintMaximumIntensities_(Size layer_index, Size rt_pixel_count, Size mz_pixel_count, QPainter& p);

        /**
        @brief Paints the precursor peaks.

        @param layer_index The index of the layer.
        @param p The QPainter to paint on.
       */
        void paintPrecursorPeaks_(Size layer_index, QPainter& painter);

        /**
        @brief Paints feature data.

        @param layer_index The index of the layer.
        @param p The QPainter to paint on.
       */
        void paintFeatureData_(Size layer_index, QPainter& p);

        /**
      	@brief Paints convex hulls (one for each mass trace) of a features layer.
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
        void paintTraceConvexHulls_(Size layer_index, QPainter& p);

        /**
      	@brief Paints the convex hulls (one for each feature) of a features layer.
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
        void paintFeatureConvexHulls_(Size layer_index, QPainter& p);

        /**
      	@brief Paints peptide identifications (for idXML and unassigned peptides in featureXML).
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
        void paintIdentifications_(Size layer_index, QPainter& p);

        /**
      	@brief Paints the consensus elements of a consensus features layer.
      	
      	@param layer_index Index of the layer.
      	@param p The QPainter to paint on.
      */
        void paintConsensusElements_(Size layer_index, QPainter& p);

        /**
      	@brief Paints one consensus element of a consensus features layer.
      	
      	@param layer_index Index of the layer.
      	@param cf Reference to the feature to be painted.
      	@param p The QPainter to paint on.
      	@param use_buffer Flag to switch between painting on the buffer and screen.
      */		
        void paintConsensusElement_(Size layer_index, const ConsensusFeature& cf, QPainter& p, bool use_buffer);

        /**
      	@brief checks if any element of a consensus feature is currently visible.
      	
      	@param layer_index Index of the layer.
      	@param ce The ConsensusFeature that needs checking
      */
        bool isConsensusFeatureVisible_(const ConsensusFeature& ce, Size layer_index);

        /**
      	@brief Paints convex hulls (one for each mass trace) for a single feature.
      	
      	@param hulls Reference to convex hull vector.
      	@param p The QPainter to paint on.
      */
        void paintConvexHulls_(const std::vector<ConvexHull2D>& hulls, QPainter& p);

        // Docu in base class
        virtual void intensityModeChange_();
        // DOcu in base class
        virtual void recalculateSnapFactor_();
        /// recalculates the dot gradient of a layer
        void recalculateDotGradient_(Size layer);

        /// m/z projection data
        ExperimentType projection_mz_;
        /// RT projection data
        ExperimentType projection_rt_;

        /**
      	@brief Returns the color associated with @p val for the gradient @p gradient.
      	
      	Takes intensity modes into account.
      */
        inline QRgb heightColor_(Real val, const MultiGradient& gradient, DoubleReal snap_factor)
        {
            switch (intensity_mode_)
            {
            case IM_NONE:
                return gradient.precalculatedColorAt(val).rgb();
                break;
            case IM_PERCENTAGE:
                return gradient.precalculatedColorAt(val*percentage_factor_).rgb();
                break;
            case IM_SNAP:
                return gradient.precalculatedColorAt(val*snap_factor).rgb();
                break;
            case IM_LOG:
                return gradient.precalculatedColorAt(std::log(val + 1)).rgb();
            default:
                throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
            }
        }

        /**
          @brief Convert chart to widget coordinates

          Translates chart coordinates to widget coordinates.
          @param x the chart coordinate x
          @param y the chart coordinate y
          @param point returned widget coordinates
        */
        inline void dataToWidget_(float x, float y, QPoint& point)
        {
          if (!isMzToXAxis())
          {
            point.setX( int((y - visible_area_.minY()) / visible_area_.height() * width()));
            point.setY(height() - int((x - visible_area_.minX()) / visible_area_.width() * height()));
          }
          else
          {
            point.setX( int((x - visible_area_.minX()) / visible_area_.width() * width()));
            point.setY( height() - int((y-visible_area_.minY())/visible_area_.height()*height()));
          }
        }

        /// Highlights a single peak and prints coordinates to screen
        void highlightPeak_(QPainter& p, const PeakIndex& peak);

        /// Returns the nearest peak to position @p pos
        PeakIndex findNearestPeak_(const QPoint& pos);

        /// Paints a peak icon for feature and consensus feature peaks
        void paintIcon_(const QPoint& pos, const QRgb& color, const String& icon, Size s, QPainter& p) const;

        /// the nearest peak/feature to the mouse cursor
        PeakIndex selected_peak_;
        /// start peak/feature of measuring mode
        PeakIndex measurement_start_;

        //docu in base class
        virtual void translateLeft_();
        //docu in base class
        virtual void translateRight_();
        //docu in base class
        virtual void translateForward_();
        //docu in base class
        virtual void translateBackward_();

        /// Finishes context menu after customization to peaks, features or consensus features
        void finishContextMenu_(QMenu* context_menu, QMenu* settings_menu);
    };
}

#endif
