// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPViewOpenDialog.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

class QAbstractButton;

namespace OpenMS 
{
	class Param;
	class String;
	/**
		@brief Dataset opening options for TOPPView
		
		@ingroup TOPPView_elements
	*/
	class OPENMS_GUI_DLLAPI TOPPViewOpenDialog
		: public QDialog,
  		public Ui::TOPPViewOpenDialogTemplate
	{
		Q_OBJECT
		
		public:			
			/// Constructor
			TOPPViewOpenDialog(const String& data_name, bool as_window, bool as_2d, bool cutoff, QWidget* parent = 0 );
			/// Destructor
			virtual ~TOPPViewOpenDialog();
			
			/// Returns true, if 2D mode is to be used for maps
			bool viewMapAs2D() const;
			/// Returns true, if 1D mode is to be used for maps
			bool viewMapAs1D() const;
			/// Returns of the low intensity peaks should be hidden
			bool isCutoffEnabled() const;
			/// Returns true, if the data should be opened in a new window
			bool openAsNewWindow() const;
			///Returns the index of the selected merge layer. If the option is not selected -1 is returned.
			Int getMergeLayer() const;
			
			/// Disables view dimension section and sets the selected option
			void disableDimension(bool as_2d);
			/// Disables cutoff section and sets the selected option
			void disableCutoff(bool cutoff_on);
			/// Disables opening location section and sets the selected option
			void disableLocation(bool window);
			/**
				@brief Sets the possible merge layers (index and name) and activates the the option
				
				It is deactivated by default and can be deactivated manually by passing an empty list.
			*/
			void setMergeLayers(const Map<Size,String>& layers);
			
		protected slots:
			///slot that disables 2D/3D options, when as layer is selected
			void updateViewMode_(QAbstractButton* button);
		
		protected:
			///Stores if this option is disabled, to avoid activating it in updateViewMode_()
			bool map_as_2d_disabled_;
	};

}
#endif // OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H

