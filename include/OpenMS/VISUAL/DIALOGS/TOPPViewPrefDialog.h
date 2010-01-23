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


#ifndef OPENMS_VISUAL_DIALOGS_TOPPVIEWPREFDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPVIEWPREFDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPViewPrefDialog.h>

namespace OpenMS 
{
	namespace Internal
	{
		/**
			@brief Preferences dialog for TOPPView
			
			@ingroup TOPPView_elements
		*/
		class OPENMS_DLLAPI TOPPViewPrefDialog
			: public QDialog,
	  		public Ui::TOPPViewPrefDialogTemplate
		{
			Q_OBJECT
			
			public:
				TOPPViewPrefDialog(QWidget * parent);
				
			protected slots:
				void browseDefaultPath_();
				void browseTempPath_();
		};
	}
}
#endif // OPENMS_VISUAL_DIALOGS_TOPPVIEWPREFDIALOG_H

