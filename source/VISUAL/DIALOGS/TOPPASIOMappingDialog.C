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

#include <QtGui/QMessageBox>

#include <iostream>
#include <sstream>

namespace OpenMS
{
	TOPPASIOMappingDialog::TOPPASIOMappingDialog(TOPPASEdge* parent)
	{
		edge_ = parent;
		setupUi(this);
		connect (ok_button,SIGNAL(clicked()),this,SLOT(checkValidity_()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
		
		fillComboBoxes_();
	}
	
	void TOPPASIOMappingDialog::fillComboBoxes_()
	{
		TOPPASVertex* source = edge_->getSourceVertex();
		TOPPASVertex* target = edge_->getTargetVertex();
		
		TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
		TOPPASToolVertex* target_tool = qobject_cast<TOPPASToolVertex*>(target);
		
		if (source_tool)
		{
			QVector<TOPPASToolVertex::IOInfo> source_output_files;
			source_tool->getOutputParameters(source_output_files);
			source_label->setText(source_tool->getName().toQString());
			if (source_tool->getType() != "")
			{
				source_type_label->setText("(" + source_tool->getType().toQString() + ")");
			}
			else
			{
				source_type_label->setVisible(false);
			}
			source_combo->addItem("<select>");
			foreach (TOPPASToolVertex::IOInfo info, source_output_files)
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
		
		if (target_tool)
		{
			QVector<TOPPASToolVertex::IOInfo> target_input_files;
			target_tool->getInputParameters(target_input_files);
			target_label->setText(target_tool->getName().toQString());
			if (target_tool->getType() != "")
			{
				target_type_label->setText("(" + target_tool->getType().toQString() + ")");
			}
			else
			{
				target_type_label->setVisible(false);
			}
			target_combo->addItem("<select>");
			foreach (TOPPASToolVertex::IOInfo info, target_input_files)
			{
				// check if parameter occupied by another edge already
				bool occupied = false;
				for (TOPPASVertex::EdgeIterator it = target->inEdgesBegin(); it != target->inEdgesEnd(); ++it)
				{
					int param_index = (*it)->getTargetInParam();
					if (param_index >= 0)
					{
						if (info.param_name == target_input_files[param_index].param_name)
						{
							occupied = true;
							break;
						}
					}
				}
				if (occupied)
				{
					continue;
				}

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
		
		int source_out = edge_->getSourceOutParam();
		int target_in = edge_->getTargetInParam();
		if (source_out != -1)
		{
			source_combo->setCurrentIndex(source_out + 1);
		}
		if (target_in != -1)
		{
			target_combo->setCurrentIndex(target_in + 1);
		}
		
		resize(width(),0);
	}
	
	void TOPPASIOMappingDialog::checkValidity_()
	{
		const QString& source_text = source_combo->currentText();
		const QString& target_text = target_combo->currentText();

		if (source_text == "<select>")
		{
			QMessageBox::warning(0,"Invalid selection","You must specify the source output parameter!");
			return;
		}
		if (target_text == "<select>")
		{
			QMessageBox::warning(0,"Invalid selection","You must specify the target input parameter!");
			return;
		}
		
		edge_->setSourceOutParam(source_combo->currentIndex()-1);
		edge_->setTargetInParam(target_combo->currentIndex()-1);
		
		edge_->updateColor();
		
		if (edge_->getEdgeStatus() == TOPPASEdge::ES_VALID ||
				edge_->getEdgeStatus() == TOPPASEdge::ES_NOT_READY_YET)
		{
			accept();
		}
	}
} // namespace
