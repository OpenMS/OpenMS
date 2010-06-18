// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>

#include <QtGui/QFileDialog>

#include <iostream>

namespace OpenMS
{
	TOPPASInputFilesDialog::TOPPASInputFilesDialog(const QStringList& list)
	{
		setupUi(this);
		
		//input_file_list->setSortingEnabled(true);
		input_file_list->addItems(list);
		
		connect (ok_button,SIGNAL(clicked()),this,SLOT(accept()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
		connect (add_button,SIGNAL(clicked()),this,SLOT(showFileDialog()));
		connect (remove_button,SIGNAL(clicked()),this,SLOT(removeSelected()));
		connect (remove_all_button,SIGNAL(clicked()),this,SLOT(removeAll()));
		connect (edit_button,SIGNAL(clicked()),this,SLOT(editCurrentItem()));
		connect (up_button,SIGNAL(clicked()),this,SLOT(moveCurrentItem()));
		connect (down_button,SIGNAL(clicked()),this,SLOT(moveCurrentItem()));
	}
	
	void TOPPASInputFilesDialog::showFileDialog()
	{
		QStringList file_names = QFileDialog::getOpenFileNames(this, tr("Select input file(s)"), tr(""), tr(/*valid filetypes*/""));
		if (!file_names.isEmpty())
		{
			input_file_list->addItems(file_names);
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
	
	void TOPPASInputFilesDialog::removeAll()
	{
		input_file_list->clear();
	}
	
	void TOPPASInputFilesDialog::getFilenames(QStringList& files)
	{
		files.clear();
		for (int i = 0; i < input_file_list->count(); ++i)
		{
			files.push_back(input_file_list->item(i)->text());
		}
	}
	
	void TOPPASInputFilesDialog::editCurrentItem()
	{
		QListWidgetItem* item = input_file_list->currentItem();
		if (!item)
		{
			return;
		}
		TOPPASInputFileDialog tifd(item->text());
		if (tifd.exec())
		{
			item->setText(tifd.getFilename());
		}
	}
	
	void TOPPASInputFilesDialog::moveCurrentItem()
	{
		if (input_file_list->count() < 2)
		{
			return;
		}
		int row = input_file_list->currentRow();
		if (row < 0)
		{
			return;
		}
		
		bool direction;
		if (QObject::sender() == up_button)
		{
			direction = true;
		}
		else if (QObject::sender() == down_button)
		{
			direction = false;
		}
		else
		{
			return;
		}
		
		if (direction == true) // move upwards
		{
			if (row == 0)
			{
				return;
			}
			QListWidgetItem* item = input_file_list->takeItem(row);
			input_file_list->insertItem(row-1, item);
			input_file_list->setCurrentItem(item);
		}
		else // move downwards
		{
			if (row == input_file_list->count()-1)
			{
				return;
			}
			QListWidgetItem* item = input_file_list->takeItem(row);
			input_file_list->insertItem(row+1, item);
			input_file_list->setCurrentItem(item);
		}
	}
	
} // namespace
