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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM2DWINDOW_H
#define OPENMS_VISUAL_SPECTRUM2DWINDOW_H

#include <OpenMS/config.h>

// OpenMS
#include <OpenMS/VISUAL/SpectrumWindow.h>
#include <OpenMS/KERNEL/DSpectrum.h>

class QPopupMenu;
class QGridLayout;

namespace OpenMS
{
	class AxisWidget;
	class Spectrum1DWidget;
	class Spectrum2DWidget;
	
	/**
		@brief Window for 2D-visualization of map data
		
		
		
		@ingroup spectrum_widgets
	*/
	class Spectrum2DWindow : public SpectrumWindow
	{
		Q_OBJECT
	public:
		/// Constructor
		Spectrum2DWindow(QWidget* parent=0, const char* name="Spectrum2DWindow", WFlags f=0);
		/// Destructor
		~Spectrum2DWindow();
		
		// Docu in base class
		Spectrum2DWidget* widget();

		/// returns the mode for 2D dots   
		SignedInt getDotMode();
    
		// Docu in base class
		virtual PreferencesDialogPage* createPreferences(QWidget* parent);

	public slots:
		/// Shows or hides the projects
		void show1DProjections(bool on);
		/// Changes the visibilty of the projections
		void changeShow1DProjections();
		// Docu in base class
    virtual void showGoToDialog();    

	protected slots:
		virtual void createContextMenu_();
	protected:
		QGridLayout* grid_;
		// Widget Data is drawn on
		Spectrum1DWidget* tic_;
		Spectrum1DWidget* projection_;
	
	private slots:
		void horizontalSpectrum(const DSpectrum<1>&);
		void verticalSpectrum(const DSpectrum<1>&);
	};
}

#endif

