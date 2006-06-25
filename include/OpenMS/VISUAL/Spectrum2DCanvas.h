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
// $Id: Spectrum2DCanvas.h,v 1.39 2006/06/08 15:51:32 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_SPECTRUM2DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM2DCANVAS_H

// OpenMS
#include <OpenMS/DATASTRUCTURES/QuadTree.h>
#include <OpenMS/KERNEL/DSpectrum.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/VISUAL/MultiGradient.h>

// QT
class QPainter;

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
		
		@todo Marching squares: wird nicht feiner beim Zoomen, log scale (Marc)
		
		@ingroup spectrum_widgets
	*/
	class Spectrum2DCanvas : public SpectrumCanvas
	{
		Q_OBJECT
		
		friend class Spectrum2DWidget;
		
		friend class Internal::Spectrum2DCanvasPDP;
		
		
	public:
		/**	@name Type definitions */
		//@{
    ///
		typedef DimensionDescription < DimensionDescriptionTagLCMS > DimensionDescription;
		///
		enum DimensionId { MZ = DimensionDescription::MZ, RT = DimensionDescription::RT };		
		///View modes for 2D dots.
		enum DotModes 
		{
			DOT_BLACK = 0,            ///< use black only
			DOT_GRADIENT = 1          ///< use gradient
		};
		//@}
		
		/**
			@brief Constructor
		
			Spectrum2DCanvas constructor. See QWidget for details.
		
			@param parent The parent widget.
			@param name The widget's name.
		*/
		Spectrum2DCanvas(QWidget* parent = 0, const char* name = "Spectrum2DCanvas");
		
		/**
			@brief Destructor
		
			Destroys the Widget and all associated data.
		*/
		~Spectrum2DCanvas();
		
		/**
			@brief Draws the contents.
		
			Device independent drawing function. Draws the contents on painter @p p.
			This function follows the 
		
			@param p The QPainter to draw the data to.
			@param width The width of the canvas in pixels.
			@param height The height of the canvas in pixels.
		*/
		void print(QPainter* p, int width, int height);
		
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
		
		
		/**
			@brief Scale dots according to intesity
		
			Sets whether or not a peak representing dot should be resized
			according to the peak's intensity.
		
			@param on if @c true, the dots are scaled, if @c false all dots are of the same size
		*/
		void setIntensityScaledDots(bool on);
		
		/**
			@brief Returns whether or not dots are scaled
		
			Returns whether or not dots are scaled according to their peak's
			intensity.
		
			@return @c true, if dots are scaled, @c false otherwise
		*/
		bool isIntensityScaledDots() { return intensity_scaled_dots_; }
		
		/**
			@brief Creates a preferences dialog page
		
			Creates a preferences dialog page for configuring this view.
			This is used by the PreferencesManager. Reimplement to provide
			an own preferences dialog page.
		
			@param parent the parent widget for the dialog page
		*/
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);
		
		// Docu in base class
		void setMainPreferences(const Param& prefs);
		
	signals:
		void selectedHorz(const DSpectrum<1>&);
		void selectedVert(const DSpectrum<1>&);
	
	public slots:
	
		void showContours(bool on);
		void showColors(bool on);
		void showPoints(bool on);

		bool getShowContours();
		bool getShowColors();
		bool getShowPoints();

		// Docu in base class
		void activateDataSet(int data_set);
		// Docu in base class
		void removeDataSet(int data_set);
		// Docu in base class
		SignedInt finishAdding();
	
	protected:
		//* @name Mouse events */
		//@{	
		virtual void mousePressEvent(QMouseEvent* e);
		virtual void mouseReleaseEvent(QMouseEvent* e);
		virtual void mouseMoveEvent(QMouseEvent* e);
		virtual void wheelEvent(QWheelEvent* e);
		//@}
		
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
			@brief Paints individual peaks.
		
			Paints the peaks as small ellipses. Peaks are
			drawn in sorted order such that lower peaks are
			drawn first. This ensures that high peaks are always
			visible. The peaks are colored according to the
			selected dot gradient.
			
			@param data_set The index of the dataset.
			@param p The QPainter to paint on.
		*/
		void paintPoints_(UnsignedInt data_set, QPainter* p);
		
		/**
			@brief Paints data as a height map.
		
			Paints the peak data as interpolated contour lines.
			The data is shown as a height map such that higher
			areas are enclosed by more lines than lower areas.
			
			@param data_set The index of the dataset.
			@param p The QPainter to paint on.
		*/
		void paintContourLines_(UnsignedInt data_set, QPainter* p);
		
		/**
			@brief Paints data as a colored surface gradient.
		
			Paints the peak data as an interpolated surface gradient.
			The data is shown according to the gradien which can be
			set with the setSurfaceGradient() member function.
			
			@param data_set The index of the dataset.
			@param p The QPainter to paint on.
		*/
		void paintColorMap_(UnsignedInt data_set, QPainter* p);
		
		// Docu in base class
		virtual void intensityModeChange_();
		// Docu in base class
		virtual void intensityDistributionChange_();
		
		/// recalculates the surface gradient inerpolation values. Use after Intensites or gradient changed
		void recalculateSurfaceGradient_();
		/// recalculates the dot gradient inerpolation values. Use after Intensites or gradient changed
		void recalculateDotGradient_();
		
		void createHorzScan_(float min, float max);
		void createVertScan_(float min, float max);
		
		typedef QuadTree<KernelTraits, PeakType > QuadTreeType_;
		
		// zooms around position pos with factor.
		void zoom_(const PointType& pos, float factor);
		// zooms in around position pos with a fixed factor.
		void zoomIn_(const PointType& pos);
		// zooms out around position pos with a fixed factor.
		void zoomOut_(const PointType& pos);
		
		// interpolation helper function
		float betweenFactor_(float v1, float v2, float val);
		/// Returns the color associated with @p val for the surface gradient. Takes Log mode into accout
		const QColor& heightColor_(float val);
		/// Performs the marching squares calculations for a dataset and stores the matrix in marching_squares_matrices_
		void calculateMarchingSquareMatrix_(UnsignedInt data_set);
		/// Returns the marching square cell with the smallest data coordinates
		AreaType getOriginCell_(UnsignedInt data_set);
		
		/// Highlights peak under cursor and start/stop peak for measurement
		void highlightPeaks_();
		/// Highlights a single peak
		void highlightPeak_(QPainter* p, DPeak<2>* peak);
		
		/// Returns the nearest peak to position @p pos
		DPeak<2>* findNearestPeak_(QPoint pos);
		
		/**
			@brief This quad tree stores the peaks which are actually shown. 
			
			It's a pointer since the constructor of a QuadTree needs the bounding area
		  of all points ever inserted, and this area is not known in the Spectrum2DCanvas constructor.
		*/
		std::vector<QuadTreeType_*> trees_;
		
		/// marching squares matrices for the datasets
		std::vector< std::vector< std::vector<float> > > marching_squares_matrices_;
		/// Contains the highes value in the marching squares matrix foreach dataset
		std::vector<float> max_values_;
		
		/// Flags whether or not to show the height map for each dataset
		std::vector<bool> show_contours_;
		/// Flags whether or not to show the surface gradient for each dataset
		std::vector<bool> show_colors_;
		/// Flags whether or not to show individual peaks for each dataset
		std::vector<bool> show_points_;
		/// Flag whether or not to show points scaled for each dataset
		bool intensity_scaled_dots_;
		
		// the nearest peak to the mouse cursor.
		DPeak<2>* nearest_peak_;
		/// start peak of measuring mode
		DPeak<2>* measurement_start_;
		/// end peak of measuring mode
		DPeak<2>* measurement_stop_;
		// temporary peak that is constructed out of the 1D Peak and the RT (for findNearestPeak_)
		DPeak<2> tmp_peak_;
		
		// Gradient for dots
		MultiGradient dot_gradient_;
		// Gradient for surface
		MultiGradient surface_gradient_;
		
	};
}

#endif
