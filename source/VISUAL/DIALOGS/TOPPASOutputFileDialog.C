// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFileDialog.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>

#include <iostream>

namespace OpenMS
{
	TOPPASOutputFileDialog::TOPPASOutputFileDialog(TOPPASOutputFileVertex* parent)
		: parent_(parent)
	{
		setupUi(this);
		
		line_edit->setText(parent->getFilename());
		
		connect (browse_button,SIGNAL(clicked()),this,SLOT(showFileDialog()));
		connect (ok_button,SIGNAL(clicked()),this,SLOT(checkValidity_()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
	}
	
	void TOPPASOutputFileDialog::showFileDialog()
	{
		QFileDialog fd;
		fd.setAcceptMode(QFileDialog::AcceptSave);
		fd.setFileMode(QFileDialog::AnyFile);
		if (File::exists(File::path(line_edit->text())))
		{
			fd.setDirectory(File::path(line_edit->text()).toQString());
		}
		//fd.setFilter("*.mzData;*.mzML;*.dta; .....");
		if (fd.exec() && !fd.selectedFiles().empty())
		{
			line_edit->setText(fd.selectedFiles().first());
		}
	}
	
	QString TOPPASOutputFileDialog::getFilename()
	{
		return line_edit->text();
	}
	
	void TOPPASOutputFileDialog::checkValidity_()
	{
		if (!(parent_->fileNameValid(line_edit->text())))
		{
			QMessageBox::warning(0,"Invalid file name","The specified file name is invalid!");
			return;
		}
		
		accept();
	}
	
} // namespace
