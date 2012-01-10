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
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>

#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>
#include <QtGui/QCompleter>
#include <QtGui/QDirModel>
#include <QtCore/QDir>
#include <QtCore/QFileInfo>

#include <iostream>

namespace OpenMS
{
	TOPPASOutputFilesDialog::TOPPASOutputFilesDialog(const QString& dir_name, int num_jobs)
	{
		setupUi(this);
		if (dir_name != "")
		{
			line_edit->setText(dir_name);
		}
		else
		{
			line_edit->setText(QDir::currentPath());
		}
    if (num_jobs >= 1)
    {
      num_jobs_box->setValue(num_jobs);
    }
		QCompleter* completer = new QCompleter(this);
		QDirModel* dir_model = new QDirModel(completer);
		dir_model->setFilter(QDir::AllDirs);
		completer->setModel(dir_model);
		line_edit->setCompleter(completer);
		connect (browse_button,SIGNAL(clicked()),this,SLOT(showFileDialog()));
		connect (ok_button,SIGNAL(clicked()),this,SLOT(checkValidity_()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
	}
	
	void TOPPASOutputFilesDialog::showFileDialog()
	{
		QString dir = File::exists(File::path(line_edit->text())) ? File::path(line_edit->text()).toQString() : "";
		QString selected_dir = QFileDialog::getExistingDirectory(this, tr("Select output directory"), dir);
		if (selected_dir != "")
		{
			line_edit->setText(selected_dir);
		}
	}
	
	QString TOPPASOutputFilesDialog::getDirectory()
	{
		return line_edit->text();
	}

  int TOPPASOutputFilesDialog::getNumJobs()
  {
    return num_jobs_box->value();
  }
	
	void TOPPASOutputFilesDialog::checkValidity_()
	{
		if (!dirNameValid(line_edit->text()))
		{
			QMessageBox::warning(0,"Invalid directory","Either the specified path is no directory, or you have no permission to write there.");
			return;
		}
		
		accept();
	}
	
	bool TOPPASOutputFilesDialog::dirNameValid(const QString& dir_name)
	{
		QFileInfo fi(dir_name);
		QString file_name = dir_name;
		if (!file_name.endsWith(QDir::separator()))
		{
		 file_name += QDir::separator();
		}
		file_name += "test_file";
		return (fi.isDir() && File::writable(file_name));
	}
	
} // namespace
