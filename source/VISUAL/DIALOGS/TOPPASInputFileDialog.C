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
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>

#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>

#include <iostream>

namespace OpenMS
{
	TOPPASInputFileDialog::TOPPASInputFileDialog(TOPPASInputFileVertex* parent)
		: parent_(parent)
	{
		setupUi(this);
		
		line_edit->setText(parent->getFilename());
		
		connect (browse_button,SIGNAL(clicked()),this,SLOT(showFileDialog()));
		connect (ok_button,SIGNAL(clicked()),this,SLOT(checkValidity_()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
	}
	
	void TOPPASInputFileDialog::showFileDialog()
	{
		QFileDialog fd;
		fd.setFileMode(QFileDialog::ExistingFile);
		//fd.setFilter("*.mzData;*.mzML;*.dta; .....");
		if (fd.exec() && !fd.selectedFiles().empty())
		{
			line_edit->setText(fd.selectedFiles().first());
		}
	}
	
	QString TOPPASInputFileDialog::getFilename()
	{
		return line_edit->text();
	}
	
	void TOPPASInputFileDialog::checkValidity_()
	{
		if (!(parent_->fileNameValid(line_edit->text())))
		{
			QMessageBox::warning(0,"Invalid file name","The specified file does not exist!");
			return;
		}
		
		accept();
	}
	
} // namespace
