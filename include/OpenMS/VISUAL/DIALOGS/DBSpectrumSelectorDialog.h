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


#ifndef OPENMS_VISUAL_DIALOGS_DBSPECTRUMSELECTORDIALOG_H
#define OPENMS_VISUAL_DIALOGS_DBSPECTRUMSELECTORDIALOG_H

#include <vector>
#include <QtGui/QDialog>
#include <OpenMS/CONCEPT/Types.h>

class QLineEdit;
class QTableWidget;

namespace OpenMS
{
	class DBConnection;
	
	/**
		@brief Dialog that allow selecting a spectrum from a DB.
		
		@ingroup Dialogs
	*/
	class DBSpectrumSelectorDialog 
		: public QDialog
	{
		Q_OBJECT
		public:
			/**
				@brief Constructor
				
				The spectrum ids to load are inserted into the @p result vector.
				
				An external DB connection is used by handing over @p adapter.
			*/
			DBSpectrumSelectorDialog(DBConnection& adapter, std::vector<UnsignedInt>& result,QWidget* parent=0);
			/// Destructor
			~DBSpectrumSelectorDialog();

		private slots:
			/// Slot for accepting the selection
			void ok();
			/// Slot for refreshing the shown spectra
			void loadSpectra();

		protected:
			/// DB connection
			DBConnection& adapter_;
			/// reference to the result vector
			std::vector<UnsignedInt>& result_;
			/// pointer to the search string lineedit
			QLineEdit* search_string_;
			/// pointer to the table for displaying the overview
			QTableWidget* table_;
	};

} //namespace

#endif //OPENMS_VISUAL_DIALOGS_DBSPECTRUMSELECTORDIALOG_H

