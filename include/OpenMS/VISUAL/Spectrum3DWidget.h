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
// $Authors: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM3DWIDGET_H

#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>

namespace OpenMS
{	
	class Spectrum3DCanvas;
	/**
		@brief Widget for 3D-visualization of map data
		
		@image html Spectrum3DWidget.png
		
		@ingroup SpectrumWidgets
	*/
	class OPENMS_DLLAPI Spectrum3DWidget
		: public SpectrumWidget
	{
		Q_OBJECT		

		public:	
			///Constructor
			Spectrum3DWidget(const Param& preferences, QWidget* parent = 0);
			
			/// Destructor
			virtual ~Spectrum3DWidget();
			
			/// This method is overwritten to make the class specific members accessable
			inline Spectrum3DCanvas* canvas()
			{
				return static_cast<Spectrum3DCanvas*>(canvas_);
			}
		
			// Docu in base class
			virtual void recalculateAxes_();
			// Docu in base class
			virtual Math::Histogram<> createIntensityDistribution_() const;
			// Docu in base class
			virtual Math::Histogram<> createMetaDistribution_(const String& name) const;
			
			//docu in base class
			bool isLegendShown() const;
			//docu in base class
			virtual void showLegend(bool show);

    signals:
      /// Requests to display all spectra in 2D plot
      void showCurrentPeaksAs2D();

		public slots:
			// Docu in base class
	    virtual void showGoToDialog();
	};
	
}//namespace

#endif
