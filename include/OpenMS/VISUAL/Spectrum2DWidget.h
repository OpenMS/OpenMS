// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_VISUAL_SPECTRUM2DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM2DWIDGET_H

// OpenMS
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>

class QGroupBox;
class QLabel;

namespace OpenMS
{
	class Spectrum2DCanvas;
	
	/**
		@brief Widget for 2D-visualization of peak map and feature map data
		
		The widget composes of two scoll bars, two AxisWidget and a Spectrum2DCanvas as central widget.
		
		@image html Spectrum2DWidget.png
		
		The example image shows %Spectrum2DWidget displaying a peak layer and a feature layer. 
		
		@ingroup SpectrumWidgets
	*/
	class Spectrum2DWidget 
		: public SpectrumWidget
	{
		Q_OBJECT
	public:

		/// Default constructor
		Spectrum2DWidget(const Param& preferences, QWidget* parent = 0);
		/// Destructor
		virtual ~Spectrum2DWidget();
		
		/// This method is overwritten to make the class specific members accessable
		inline Spectrum2DCanvas* canvas()
		{
			return static_cast<Spectrum2DCanvas*>(canvas_);
		}

		/// const reference to the horizontal projection
		const Spectrum1DWidget* getHorizontalProjection() const;
		/// const reference to the vertical projection
		const Spectrum1DWidget* getVerticalProjection() const;

	public slots:	
		// Docu in base class
		virtual void recalculateAxes_();
		/// Hides the projections
		void hideProjections();
		// Docu in base class
    virtual void showGoToDialog();
    
	signals:
		/**
			@brief Signal emitted whenever the visible area changes.
			
			@param area The new visible area.
		*/
		void visibleAreaChanged(DRange<2> area);
		/// Requests to display the spectrum with index @p index in 1D
		void showSpectrumAs1D(int index);

	protected:
		// Docu in base class
		virtual Math::Histogram<UInt,float> createIntensityDistribution_();
		
		/// Vertical projection widget
		Spectrum1DWidget* projection_vert_;
		/// Horizontal projection widget
		Spectrum1DWidget* projection_horz_;
		/// Group box that shows information about the projections
		QGroupBox* projection_box_;
		/// Number of peaks of the projection
		QLabel* projection_peaks_;
		/// Intensity sum of the projection
		QLabel* projection_sum_;
		/// Intensity maximum of the projection
		QLabel* projection_max_;

	private slots:
		/// shows the horizontal projection with the given data and draw mode
		void horizontalProjection(const MSExperiment<>&, Spectrum1DCanvas::DrawModes);
		/// shows the vertical projection with the given data and draw mode
		void verticalProjection(const MSExperiment<>&, Spectrum1DCanvas::DrawModes);
		/// shows projections information
		void projectionInfo(int peaks, double intensity, double max);

	};
}

#endif
