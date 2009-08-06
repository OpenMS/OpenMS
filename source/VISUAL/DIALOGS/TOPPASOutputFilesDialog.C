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
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFileDialog.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>

#include <QtGui/QFileDialog>

#include <iostream>

namespace OpenMS
{
	TOPPASOutputFilesDialog::TOPPASOutputFilesDialog(const QStringList& list)
	{
		setupUi(this);
		
		output_file_list->setSelectionMode(QAbstractItemView::SingleSelection);
		output_file_list->setSortingEnabled(false); // same order as output files of tool
		
		output_file_list->addItems(list);
		
		connect (ok_button,SIGNAL(clicked()),this,SLOT(checkValidity_()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
		connect (edit_button,SIGNAL(clicked()),this,SLOT(editCurrentItem()));
	}
	
	void TOPPASOutputFilesDialog::editCurrentItem()
	{
		QListWidgetItem* item = output_file_list->currentItem();
		if (!item)
		{
			return;
		}
		
		TOPPASOutputFileDialog tofd(item->text());
		if (tofd.exec())
		{
			item->setText(tofd.getFilename());
		}
	}
	
	
	void TOPPASOutputFilesDialog::getFilenames(QStringList& files)
	{
		files.clear();
		for (int i = 0; i < output_file_list->count(); ++i)
		{
			const QString& text = output_file_list->item(i)->text();
			if (text != "<edit filename>")
			{
				files.push_back(text);
			}
		}
	}
	
	void TOPPASOutputFilesDialog::checkValidity_()
	{
		// ...
		accept();
	}
	
} // namespace
