// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  this library is distributed in the hope that it will be useful,
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

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/DocumentIdentifierVisualizer.h>
#include <OpenMS/FORMAT/FileHandler.h>

//QT
#include <QtGui/QLineEdit>

using namespace std;

namespace OpenMS
{

	DocumentIdentifierVisualizer::DocumentIdentifierVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<DocumentIdentifier>()
	{
		addLabel_("Modify DocumentIdentifier information");
		addSeparator_();
		addLineEdit_(identifier_, "Identifier");
		addSeparator_();
		addLineEdit_(file_path_,"Loaded from file");
		addLineEdit_(file_type_,"File type");
		finishAdding_();
	}


	void DocumentIdentifierVisualizer::update_()
	{
	  identifier_->setText(temp_.getIdentifier().c_str());
	  file_path_->setText(temp_.getLoadedFilePath().c_str());
	  file_type_->setText(FileHandler::typeToName(temp_.getLoadedFileType()).c_str());
		file_path_->setReadOnly(true);
		file_type_->setReadOnly(true);
	}

	void DocumentIdentifierVisualizer::store()
	{
		ptr_->setIdentifier(identifier_->text());

		temp_ = (*ptr_);
	}

	void DocumentIdentifierVisualizer::undo_()
	{
		update_();
	}

}
