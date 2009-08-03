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

#include <QtGui/QFileDialog>

#include <iostream>

namespace OpenMS
{
	TOPPASOutputFilesDialog::TOPPASOutputFilesDialog(TOPPASOutputFileListVertex* parent)
	{
		setupUi(this);
		
		output_file_list->setSortingEnabled(false); // same order as output files of tool
		
		// output vertex has exactly 1 in edge:
		TOPPASEdge* in_edge = *(parent->inEdgesBegin());
		TOPPASToolVertex* in_tool = qobject_cast<TOPPASToolVertex*>(in_edge->getSourceVertex());
		const QVector<QStringList>& files_vector = in_tool->getOutputFileNames();
		int param_index = in_edge->getSourceOutParam();
		if (param_index != -1)
		{
			const QStringList& files = files_vector[param_index];
			output_file_list->addItems(files);
		}
		
		connect (ok_button,SIGNAL(clicked()),this,SLOT(checkValidity_()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
	}
	
// 	void TOPPASOutputFilesDialog::showFileDialog()
// 	{
// 		QFileDialog fd;
// 		fd.setFileMode(QFileDialog::ExistingFiles);
// 		//fd.setFilter("*.mzData;*.mzML;*.dta; .....");
// 		if (fd.exec())
// 		{
// 			input_file_list->addItems(fd.selectedFiles());
// 		}
// 	}
	
	
	void TOPPASOutputFilesDialog::getFilenames(QStringList& files)
	{
		files.clear();
		for (int i = 0; i < output_file_list->count(); ++i)
		{
			files.push_back(output_file_list->item(i)->text());
		}
	}
	
	void TOPPASOutputFilesDialog::checkValidity_()
	{
		// ...
		accept();
	}
	
} // namespace
