// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>

#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QCompleter>
#include <QtGui/QDirModel>
#include <QtCore/QFileInfo>

#include <iostream>

namespace OpenMS
{
	TOPPASInputFileDialog::TOPPASInputFileDialog(const QString& file_name)
	{
		setupUi(this);
		
		line_edit->setText(file_name);
		// disable completer for windows, causes crashes
		#ifndef OPENMS_WINDOWSPLATFORM
		QCompleter* completer = new QCompleter(this);
		completer->setModel(new QDirModel(completer));
		line_edit->setCompleter(completer);
		#endif
		connect (browse_button,SIGNAL(clicked()),this,SLOT(showFileDialog()));
		connect (ok_button,SIGNAL(clicked()),this,SLOT(checkValidity_()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
	}
	
	void TOPPASInputFileDialog::showFileDialog()
	{
		QString file_name = QFileDialog::getOpenFileName(this, tr("Specify input file"), tr(""), tr(/*valid formats*/""));
		if (file_name != "")
		{
			line_edit->setText(file_name);
		}
	}
	
	QString TOPPASInputFileDialog::getFilename()
	{
		return line_edit->text();
	}
	
	void TOPPASInputFileDialog::checkValidity_()
	{
			// we ALLOW non-existing filenames (e.g. for FASTA files, which
			// are searched in other paths via OpenMS.ini:id_db_dir
		if (!fileNameValid(line_edit->text()))
		{
      QMessageBox::warning(0,"Invalid file name","Warning: filename does not exist!");
		}
		
		accept();
	}
	
	bool TOPPASInputFileDialog::fileNameValid(const QString& file_name)
	{
		QFileInfo fi(file_name);
		return (fi.exists() && fi.isReadable() && (!fi.isDir()));
	}
	
} // namespace
