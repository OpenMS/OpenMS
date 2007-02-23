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

#ifndef OPENMS_VISUAL_SPECTRUM1DWIDGET_H
#define OPENMS_VISUAL_SPECTRUM1DWIDGET_H

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

class QAction;

namespace OpenMS
{
	class Spectrum1DCanvas;

	/**
		@brief Widget for visualization of spectrum
		
		@ingroup spectrum_widgets
	*/
	class Spectrum1DWidget : public SpectrumWidget
	{
		Q_OBJECT
		
	public:
		/// Default constructor
		Spectrum1DWidget(QWidget* parent = 0);
		///Destructor
		virtual ~Spectrum1DWidget();
		
		// Docu in base class
		Spectrum1DCanvas* canvas();
		
		// Docu in base class
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);
	
	signals:
		/// Is emitted whenever the visible area changes.		
		void visibleAreaChanged(double, double); 

	protected:
		// Docu in base class
		virtual Math::Histogram<UnsignedInt,float> createIntensityDistribution_();
		// Docu in base class
		virtual void recalculateAxes_();
	};
} // namespace OpenMS

#endif
