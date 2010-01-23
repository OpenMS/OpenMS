// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_DEMODIALOG_H
#define OPENMS_VISUAL_DIALOGS_DEMODIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_DemoDialog.h>

#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <QtGui/QDialog>

namespace OpenMS
{
	/**
		@brief Dialog that is used to browse a list of HTML pages, e.g. for a demo or tutorial.
		
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI DemoDialog 
		: public QDialog,
  		public Ui::DemoDialogTemplate
	{
		Q_OBJECT
		
		public:
			
			/// Constructor
			DemoDialog(QWidget* parent);
			/// Destructor
			~DemoDialog();
			
			/// Sets the window title
			void setTitle(const String& title);
			/// Sets the ordered list pages (file names) and shows the first page
			void setPages(const StringList& pages);
			
		protected:
			
			/// Ordered list of pages
			StringList pages_;
			/// The index of the current page
			Size current_page_;
			/// window title
			String window_title_;
			
			/// Shows page @p i
			void show_(Size i);
			
		protected slots:
			
			/// show previous page, if there is one
			void previous_();
			/// show next page, if there is one
			void next_();
		
		private:
			
			/// not implemented
			DemoDialog();

	};

} //namespace

#endif //OPENMS_VISUAL_DIALOGS_DEMODIALOG_H

