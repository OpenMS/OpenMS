// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free DataProcessing Foundation; either
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
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/DataProcessingVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QListWidget>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

	DataProcessingVisualizer::DataProcessingVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<DataProcessing>()
	{
		addLabel_("Modify data processing information.");	
		addSeparator_();  
		
		addLineEdit_(completion_time_, "Completion time" );
		addListView_(actions_,"Processing actions");
		finishAdding_();
	}
	
	void DataProcessingVisualizer::update_()
	{
		//time
		completion_time_->setText(temp_.getCompletionTime().get().c_str()); 
		
		//actions
		actions_->clear();
		for (Size i=0; i<DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
		{
			QListWidgetItem* item = new QListWidgetItem(actions_);
			item->setText(QString::fromStdString(DataProcessing::NamesOfProcessingAction[i]));
			if (temp_.getProcessingActions().count(DataProcessing::ProcessingAction(i))==1)
			{
				item->setCheckState(Qt::Checked);
			}
			else
			{
				item->setCheckState(Qt::Unchecked);
			}
			if (isEditable())
			{
				item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
			}
			else
			{
				item->setFlags(Qt::ItemIsEnabled);
			}
			actions_->addItem(item);
		}
	}
	
	void DataProcessingVisualizer::store()
	{
		DateTime date;
		try
		{
			date.set(completion_time_->text());
			ptr_->setCompletionTime(date);
		}
		catch(exception& /*e*/)
		{
			if(date.isNull())
			{
				std::string status= "Format of date in DATAPROCESSING is not correct.";
				emit sendStatus(status);
			}
		}
		
		//actions
		ptr_->getProcessingActions().clear();
		for (UInt i=0; i<DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
		{
			if (actions_->item(i)->checkState()==Qt::Checked)
			{
				ptr_->getProcessingActions().insert(DataProcessing::ProcessingAction(i));
			}
		}
		
		temp_=(*ptr_);
	}
	
	void DataProcessingVisualizer::undo_()
	{
		update_();
	}

}
