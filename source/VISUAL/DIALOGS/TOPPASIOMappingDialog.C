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

#include <QtGui/QTableWidgetItem>

#include <iostream>

namespace OpenMS
{
	TOPPASIOMappingDialog::TOPPASIOMappingDialog(TOPPASEdge* parent)
	{
		edge_ = parent;
		setupUi(this);
		connect (ok_button,SIGNAL(clicked()),this,SLOT(accept()));
		connect (cancel_button,SIGNAL(clicked()),this,SLOT(reject()));
		resizeEvent(0);
		
		fillTable_();
	}
	
	void TOPPASIOMappingDialog::resizeEvent(QResizeEvent* /*event*/)
	{
		if (table->columnCount() == 2)
		{
			int width = table->width() / 2;
			table->setColumnWidth(0, width);
			table->setColumnWidth(1, width);
		}
	}
	
	void TOPPASIOMappingDialog::determineEdgeType_()
	{
		bool source_vertex_is_a_tool = false;
		bool source_vertex_is_a_list = false;
		TOPPASVertex* source = edge_->getSourceVertex();
		TOPPASVertex* target = edge_->getTargetVertex();
		
		if (qobject_cast<TOPPASToolVertex*>(source))
		{
			source_vertex_is_a_tool = true;
			qobject_cast<TOPPASToolVertex*>(source)->getRequiredOutputFiles(source_output_files_);
		}
		else if (qobject_cast<TOPPASInputFileListVertex*>(source))
		{
			source_vertex_is_a_list = true;
			// fill that
		}
		if (source_vertex_is_a_tool)
		{
			if (qobject_cast<TOPPASToolVertex*>(target))
			{
				edge_type_ = ET_TOOL_TO_TOOL;
				qobject_cast<TOPPASToolVertex*>(target)->getRequiredInputFiles(target_input_files_);
			}
			else if (qobject_cast<TOPPASOutputFileListVertex*>(target))
			{
				edge_type_ = ET_TOOL_TO_LIST;
				// here too
			}
			else
			{
				edge_type_ = ET_TOOL_TO_FILE;
				// and so on
			}
		}
		else if (qobject_cast<TOPPASToolVertex*>(target))
		{
			qobject_cast<TOPPASToolVertex*>(target)->getRequiredInputFiles(target_input_files_);
			
			if (source_vertex_is_a_list)
			{
				edge_type_ = ET_LIST_TO_TOOL;
			}
			else
			{
				edge_type_ = ET_FILE_TO_TOOL;
			}
		}
		else
		{
			edge_type_ = ET_INVALID;
		}
	}
	
	void TOPPASIOMappingDialog::fillTable_()
	{
		determineEdgeType_();
		
		int counter;
		Size overall_size, i;
		switch (edge_type_)
		{
			case ET_FILE_TO_TOOL:
				//bla
				break;
			
			case ET_LIST_TO_TOOL:
				//bla
				break;
			
			case ET_TOOL_TO_FILE:
				//bla
				break;
			
			case ET_TOOL_TO_LIST:
				//bla
				break;
			
			case ET_TOOL_TO_TOOL:
				// future: store mapping of IOInfos / index (row) and vice versa [store this in edge itself?]
				overall_size = source_output_files_.size() > target_input_files_.size() ?
					source_output_files_.size() : target_input_files_.size();
				for (i = 0; i < overall_size; ++i)
				{
					table->insertRow(i);
				}
				counter = 0;
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
					item_name += info.param_name;
					
					QTableWidgetItem* item = new QTableWidgetItem(item_name.toQString());
					table->setItem(counter++,0,item);
				}
				counter = 0;
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
					item_name += info.param_name;
					
					QTableWidgetItem* item = new QTableWidgetItem(item_name.toQString());
					table->setItem(counter++,1,item);
				}
				break;
				
			case ET_INVALID:
				break;
			default:
				break;
				// error message
		}
	}
	
} // namespace
