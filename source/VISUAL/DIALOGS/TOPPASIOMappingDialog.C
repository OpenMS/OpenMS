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
	
	void TOPPASIOMappingDialog::checkValidity_()
	{
		TOPPASVertex* source = edge_->getSourceVertex();
		TOPPASVertex* target = edge_->getTargetVertex();
		const QString& source_text = source_combo->currentText();
		const QString& target_text = target_combo->currentText();

		if (source_text == "<select>" || target_text == "<select>")
		{
			QMessageBox::information(0,"Invalid selection","You must specify the output and input parameter!");
			return;
		}
		
		StringList source_param_types;
		StringList target_param_types;
		bool source_param_has_list_type = false;
		bool target_param_has_list_type = false;
		bool valid = false;
		if (qobject_cast<TOPPASToolVertex*>(source))
		{
			int source_param_index = source_combo->currentIndex()-1;
			TOPPASToolVertex::IOInfo& source_param = source_output_files_[source_param_index];
			source_param_types = source_param.valid_types;
			source_param_has_list_type = source_param.type == TOPPASToolVertex::IOInfo::IOT_LIST;
		}
		if (qobject_cast<TOPPASToolVertex*>(target))
		{
			int target_param_index = target_combo->currentIndex()-1;
			TOPPASToolVertex::IOInfo& target_param = target_input_files_[target_param_index];
			target_param_types = target_param.valid_types;
			target_param_has_list_type = target_param.type == TOPPASToolVertex::IOInfo::IOT_LIST;
		}
		if (edge_->getEdgeType() == TOPPASEdge::ET_FILE_TO_TOOL)
		{
			if (target_param_has_list_type)
			{
				QMessageBox::information(0,"Invalid selection","The selected target input parameter is of type 'list', but must be 'file'");
				return;
			}
			if (target_param_types.empty())
			{
				// no restrictions specified
				valid = true;
			}
			else
			{
				const String& file_name = String(qobject_cast<TOPPASInputFileVertex*>(source)->getFilename());
				String::SizeType extension_start_index = file_name.rfind(".");
				if (extension_start_index != String::npos)
				{
					const String& extension = file_name.substr(extension_start_index+1);
					for (StringList::iterator it = target_param_types.begin(); it != target_param_types.end(); ++it)
					{
						if (*it == extension)
						{
							valid = true;
							break;
						}
					}
				}
			}
		}
		else if (edge_->getEdgeType() == TOPPASEdge::ET_LIST_TO_TOOL)
		{
			if (!target_param_has_list_type)
			{
				QMessageBox::information(0,"Invalid selection","The selected target input parameter is of type 'file', but must be 'list'");
				return;
			}
			if (target_param_types.empty())
			{
				// no restrictions specified
				valid = true;
			}
			else
			{
				const QStringList& file_names = qobject_cast<TOPPASInputFileListVertex*>(source)->getFilenames();
				foreach (const QString& q_file_name, file_names)
				{
					const String& file_name = String(q_file_name);
					String::SizeType extension_start_index = file_name.rfind(".");
					if (extension_start_index != String::npos)
					{
						const String& extension = file_name.substr(extension_start_index+1);
						for (StringList::iterator it = target_param_types.begin(); it != target_param_types.end(); ++it)
						{
							if (*it == extension)
							{
								valid = true;
								break;
							}
						}
					}
				}
			}
		}
		else if (edge_->getEdgeType() == TOPPASEdge::ET_TOOL_TO_FILE)
		{
			if (source_param_has_list_type)
			{
				QMessageBox::information(0,"Invalid selection","The selected source output parameter is of type 'list', but must be 'file'");
				return;
			}
			
			valid = true;
		}
		else if (edge_->getEdgeType() == TOPPASEdge::ET_TOOL_TO_LIST)
		{
			if (!source_param_has_list_type)
			{
				QMessageBox::information(0,"Invalid selection","The selected source output parameter is of type 'file', but must be 'list'");
				return;
			}
			
			valid = true;
		}
		else if (edge_->getEdgeType() == TOPPASEdge::ET_TOOL_TO_TOOL)
		{
			// check
			valid = true;
		}
		else
		{
			QMessageBox::critical(0,"Invalid edge type","This should not have happened. Please contact the OpenMS mailing list and report this bug.");
			return;
		}
		
		if (valid)
		{
			accept();
		}
		else
		{
			QMessageBox::information(0,"Invalid selection","The file types/extensions of output and input do not match!");
		}
	}
	
} // namespace
