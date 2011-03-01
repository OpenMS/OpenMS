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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM1DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM1DWIDGET_H

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

class QAction;
class QSpacerItem;

namespace OpenMS
{
	class Spectrum1DCanvas;

	/**
		@brief Widget for visualization of several spectra
		
		The widget composes of a scoll bar, an AxisWidget and a Spectrum1DCanvas as central widget.
		
		@image html Spectrum1DWidget.png
		
		The example image shows %Spectrum1DWidget displaying a raw data layer and a peak data layer. 
		
		@ingroup SpectrumWidgets
	*/
	class OPENMS_DLLAPI Spectrum1DWidget 
		: public SpectrumWidget
	{
		Q_OBJECT
		
	public:
		/// Default constructor
		Spectrum1DWidget(const Param& preferences, QWidget* parent = 0);
		///Destructor
		virtual ~Spectrum1DWidget();
		
		/// This method is overwritten to make the class specific members accessable
		inline Spectrum1DCanvas* canvas()
		{
			return static_cast<Spectrum1DCanvas*>(canvas_);
		}
		
		// Docu in base class
		virtual void hideAxes();
		
		// Docu in base class
		virtual void showLegend(bool show);
		
		/// Switches to mirror view, displays another y-axis for the second spectrum
		void toggleMirrorView(bool mirror);
		
		/// Performs an alignment of the layers with @p layer_index_1 and @p layer_index_2
		void performAlignment(Size layer_index_1, Size layer_index_2, const Param& param);
		
		/// Resets the alignment
		void resetAlignment();
		
    // Docu in base class
    virtual void saveAsImage();

	signals:
		/// Is emitted whenever the visible area changes.		
		void visibleAreaChanged(double, double); 

    /// Requests to display the whole spectrum in 2D view
    void showCurrentPeaksAs2D();

    /// Requests to display the whole spectrum in 3D view
    void showCurrentPeaksAs3D();

	public slots:
		// Docu in base class
    virtual void showGoToDialog();

	protected:
		// Docu in base class
		virtual Math::Histogram<> createIntensityDistribution_() const;
		// Docu in base class
		virtual Math::Histogram<> createMetaDistribution_(const String& name) const;
		// Docu in base class
		virtual void recalculateAxes_();
		
		/// The second y-axis for the mirror view
		AxisWidget* flipped_y_axis_;
		
		/// Spacer between the two y-axes in mirror mode (needed when visualizing an alignment)
		QSpacerItem* spacer_;

	};
} // namespace OpenMS

#endif
