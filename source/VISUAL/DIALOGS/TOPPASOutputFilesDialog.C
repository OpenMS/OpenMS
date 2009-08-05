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
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>

#include <QtGui/QFileDialog>

#include <iostream>

namespace OpenMS
{
	TOPPASOutputFilesDialog::TOPPASOutputFilesDialog(TOPPASOutputFileListVertex* parent)
		: parent_(parent)
	{
		setupUi(this);
		
		output_file_list->setSelectionMode(QAbstractItemView::SingleSelection);
		output_file_list->setSortingEnabled(false); // same order as output files of tool
		
		// output vertex has exactly 1 in edge:
		TOPPASEdge* in_edge = *(parent->inEdgesBegin());
		TOPPASToolVertex* in_tool = qobject_cast<TOPPASToolVertex*>(in_edge->getSourceVertex());
		const QVector<QStringList>& files_vector = in_tool->getOutputFileNames();
		int param_index = in_edge->getSourceOutParam();
		if (param_index != -1)
		{
			const QStringList& files = files_vector[param_index];
			
			int specified_files_count = parent->getFilenames().size();
			int tmp_files_count = files.size();
			if (specified_files_count <= tmp_files_count)
			{
				output_file_list->addItems(parent->getFilenames());
				// if too few file names specified, fill the rest
				for (int i = specified_files_count; i < tmp_files_count; ++i)
				{
					output_file_list->addItem("<edit filename>");
				}
			}
			else
			{
				// too many file names specified, only show as many as needed
				const QStringList& save_names = parent->getFilenames();
				for (int i = 0; i < tmp_files_count; ++i)
				{
					output_file_list->addItem(save_names[i]);
				}
			}
		}
		
		connect (ok_button,SIGNAL(clicked()),this,SLOT(checkValidity_()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
		connect (edit_button,SIGNAL(clicked()),this,SLOT(showFileDialog()));
	}
	
	void TOPPASOutputFilesDialog::showFileDialog()
	{
		if (output_file_list->selectedItems().empty())
		{
			return;
		}
		
		QFileDialog fd;
		fd.setAcceptMode(QFileDialog::AcceptSave);
		fd.setFileMode(QFileDialog::AnyFile);
		if (fd.exec() && !fd.selectedFiles().empty())
		{
			(output_file_list->selectedItems().first())->setText(fd.selectedFiles().first());
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
