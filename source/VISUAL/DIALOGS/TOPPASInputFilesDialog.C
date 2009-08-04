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
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>

#include <QtGui/QFileDialog>

#include <iostream>

namespace OpenMS
{
	TOPPASInputFilesDialog::TOPPASInputFilesDialog(TOPPASInputFileListVertex* parent)
		: parent_(parent)
	{
		setupUi(this);
		
		input_file_list->setSortingEnabled(true);
		input_file_list->addItems(parent->getFilenames());
		
		connect (ok_button,SIGNAL(clicked()),this,SLOT(accept()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
		connect (add_button,SIGNAL(clicked()),this,SLOT(showFileDialog()));
		connect (remove_button,SIGNAL(clicked()),this,SLOT(removeSelected()));
	}
	
	void TOPPASInputFilesDialog::showFileDialog()
	{
		QFileDialog fd;
		fd.setFileMode(QFileDialog::ExistingFiles);
		//fd.setFilter("*.mzData;*.mzML;*.dta; .....");
		if (fd.exec())
		{
			input_file_list->addItems(fd.selectedFiles());
		}
	}
	
	void TOPPASInputFilesDialog::removeSelected()
	{
		QList<QListWidgetItem*> selected_items = input_file_list->selectedItems();
		foreach (QListWidgetItem* item, selected_items)
		{
			input_file_list->takeItem(input_file_list->row(item));
		}
	}
	
	void TOPPASInputFilesDialog::getFilenames(QStringList& files)
	{
		files.clear();
		for (int i = 0; i < input_file_list->count(); ++i)
		{
			files.push_back(input_file_list->item(i)->text());
		}
	}
	
} // namespace
