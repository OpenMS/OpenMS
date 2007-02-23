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

#ifndef OPENMS_VISUAL_SPECTRUMWINDOW_H
#define OPENMS_VISUAL_SPECTRUMWINDOW_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/PreferencesManager.h>

//STL
#include <string>

//QT
#include <QtGui/QMainWindow>

namespace OpenMS 
{
	class SpectrumWidget;
	
	/**
		@brief Base class for MDI windows
		
		This class is the base class for the different MDI window
		types in the TOPPView application. For each type of spectrum
		view (such as 1D view, 2D view etc.), there should be a
		corresponding class derived from this class.
		
		To integrate a new spectrum view (i.e. classes derived from
		SpectrumWidget and SpectrumCanvas) into the TOPPView application,
		another class must be derived from this class which holds an
		instance of the SpectrumWidget class as a child widget (see
		setWidget_() and widget()).
		
		@ingroup spectrum_widgets
	*/
	class SpectrumWindow 
		: public QMainWindow, 
			public PreferencesManager
	{
		Q_OBJECT
		public:
			
			/// Constructor
			SpectrumWindow(QWidget* parent=0);
			/// Destructor
			virtual ~SpectrumWindow();
			
			///connect the signals/slots of window and widget (status messages, mode changes)
			void connectWidgetSignals(SpectrumWidget* sw);

			// Docu in base class
			virtual PreferencesDialogPage* createPreferences(QWidget* parent)=0;
			
			/// Returns a pointer to the child widget
			SpectrumWidget* widget();

			///Set the main Param object
			void setMainPreferences(const Param& prefs);
			
			int window_id;
			
		signals:
			/// Display a status message. See TOPPViewBase::showStatusMessage .
			void sendStatusMessage(std::string,OpenMS::UnsignedInt);
			/// Display coordinates (mz, rt, intensity)
			void sendCursorStatus(double mz=-1.0, double intens=-1.0, double rt=-1.0);
		  ///signals that draw or display mode changed (e.g. used to update the tool bar)
		  void modesChanged(QWidget*);
		  /// Shows the main preferences dialog
		  void openPreferences();
		  /// Message about the destruction of this widget
		  void aboutToBeDestroyed(int window_id);
			
		public slots:
			/// Displays a status message. See TOPPViewBase::showStatusMessage .
			void showStatusMessage(std::string,OpenMS::UnsignedInt);
			/// Displays coordinates (mz, rt, intensity)
			void showCursorStatus(double mz, double intens, double rt);
			/// Emits the modesChanged signal
			void modesChangedSlot(QWidget*);
			/// Displays a GoTo dialog
      virtual void showGoToDialog() = 0;		      

		protected:
			/// Sets the pointer to the child widget and the back pointer
			void setWidget_(SpectrumWidget* widget);
      /// Pointer to the child widget
			SpectrumWidget* widget_;
	};
}
#endif
