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


#ifndef OPENMS_VISUAL_DIALOGS_DATAFILTERDIALOG_H
#define OPENMS_VISUAL_DIALOGS_DATAFILTERDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_DataFilterDialog.h>
#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>


namespace OpenMS 
{
	/**
		@brief Dialog for creating and changing a DataFilter
		
	*/
	class OPENMS_GUI_DLLAPI DataFilterDialog
		: public QDialog,
  		public Ui::DataFilterDialogTemplate
	{
		Q_OBJECT
		
		public:
			/// constructor
			DataFilterDialog(DataFilters::DataFilter& filter, QWidget* parent);
		
		protected slots:
			/// Checks if the settings are valid and writes them to filter_ if so
			void check_();
			/// Is called when field_ changes and enables/disables the meta data functionality as needed
			void field_changed_(const QString&);
			/// Is called when op_ changes and disables the value field, if operation is "exists", else enables it
			void op_changed_(const QString&);
		
		protected:
			/// Reference to the filter that is modified
			DataFilters::DataFilter& filter_;
		
		private:
			///Not implemented
			DataFilterDialog();
};

}
#endif // OPENMS_VISUAL_DIALOGS_OPENDIALOG_H

