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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPViewPrefDialog.h>
#include <QtGui/QFileDialog>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
		TOPPViewPrefDialog::TOPPViewPrefDialog(QWidget * parent)
			: QDialog(parent)
		{
			setupUi(this);
			connect(findChild<QPushButton*>("browse_default"),SIGNAL(clicked()),this,SLOT(browseDefaultPath_()));
			connect(findChild<QPushButton*>("browse_temp"),SIGNAL(clicked()),this,SLOT(browseTempPath_()));
		}
	
		void TOPPViewPrefDialog::browseDefaultPath_()
		{
			QString path = QFileDialog::getExistingDirectory(this, "Choose a directory", findChild<QLineEdit*>("default_path")->text());
			if (path!="")
			{
				findChild<QLineEdit*>("default_path")->setText(path);
			}
		}
		
		void TOPPViewPrefDialog::browseTempPath_()
		{
			QString path = QFileDialog::getExistingDirectory(this, "Choose a directory", findChild<QLineEdit*>("temp_path")->text());
			if (path!="")
			{
				findChild<QLineEdit*>("temp_path")->setText(path);
			}
		}

	} //namespace Internal
}//namspace OpenMS



