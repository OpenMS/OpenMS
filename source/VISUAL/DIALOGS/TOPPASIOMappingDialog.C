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
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/TOPPASInputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>

#include <iostream>
#include <string>
#include <sstream>

namespace OpenMS
{
	TOPPASIOMappingDialog::TOPPASIOMappingDialog(TOPPASEdge* parent)
	{
		edge_ = parent;
		setupUi(this);
		connect (ok_button,SIGNAL(clicked()),this,SLOT(accept()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
		
		fillComboBoxes_();
	}
	
	void TOPPASIOMappingDialog::fillComboBoxes_()
	{
		TOPPASVertex* source = edge_->getSourceVertex();
		TOPPASVertex* target = edge_->getTargetVertex();
		if (qobject_cast<TOPPASToolVertex*>(source))
		{
			qobject_cast<TOPPASToolVertex*>(source)->getRequiredOutputFiles(source_output_files_);
			source_label->setText(source->getName().toQString());
			if (source->getType() != "")
			{
				source_type_label->setText("(" + source->getType().toQString() + ")");
			}
			else
			{
				source_type_label->setVisible(false);
			}
			source_combo->addItem("<select>");
			foreach (TOPPASToolVertex::IOInfo info, source_output_files_)
			{
				String item_name;
				if (info.type == TOPPASToolVertex::IOInfo::IOT_FILE)
				{
					item_name = "File: ";
				}
				else
				{
					item_name = "List: ";
				}
				item_name += info.param_name + " ";
				std::ostringstream ss;
				ss << info.valid_types;
				item_name += ss.str();
				
				source_combo->addItem(item_name.toQString());
			}
		}
		else if (edge_->getEdgeType() == TOPPASEdge::ET_FILE_TO_TOOL)
		{
			source_label->setText("File");
			source_type_label->setVisible(false);
			source_combo->setVisible(false);
			source_parameter_label->setVisible(false);
		}
		else if (edge_->getEdgeType() == TOPPASEdge::ET_LIST_TO_TOOL)
		{
			source_label->setText("List");
			source_type_label->setVisible(false);
			source_combo->setVisible(false);
			source_parameter_label->setVisible(false);
		}
		
		if (qobject_cast<TOPPASToolVertex*>(target))
		{
			qobject_cast<TOPPASToolVertex*>(target)->getRequiredInputFiles(target_input_files_);
			target_label->setText(target->getName().toQString());
			if (target->getType() != "")
			{
				target_type_label->setText("(" + target->getType().toQString() + ")");
			}
			else
			{
				target_type_label->setVisible(false);
			}
			target_combo->addItem("<select>");
			foreach (TOPPASToolVertex::IOInfo info, target_input_files_)
			{
				String item_name;
				if (info.type == TOPPASToolVertex::IOInfo::IOT_FILE)
				{
					item_name = "File: ";
				}
				else
				{
					item_name = "List: ";
				}
				item_name += info.param_name + " ";
				std::ostringstream ss;
				ss << info.valid_types;
				item_name += ss.str();
				
				target_combo->addItem(item_name.toQString());
			}
		}
		else if (edge_->getEdgeType() == TOPPASEdge::ET_TOOL_TO_FILE)
		{
			target_label->setText("File");
			target_type_label->setVisible(false);
			target_combo->setVisible(false);
			target_parameter_label->setVisible(false);
		}
		else if (edge_->getEdgeType() == TOPPASEdge::ET_TOOL_TO_LIST)
		{
			target_label->setText("List");
			target_type_label->setVisible(false);
			target_combo->setVisible(false);
			target_parameter_label->setVisible(false);
		}
		resize(width(),0);
		//updateGeometry();
	}
	
} // namespace
