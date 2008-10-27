// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/DataProcessingVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

	DataProcessingVisualizer::DataProcessingVisualizer(bool editable, QWidget *parent)
		: BaseVisualizer(editable, parent)
	{
		type_="DataProcessing";
	  
		addLabel("Modify data processing information.");	
		addSeparator();  
		
		addLineEdit(software_completion_time_, "Completion time" );
		
		finishAdding_();
		
	}
	
	void DataProcessingVisualizer::load(DataProcessing &s)
	{
	  //Pointer to current object to keep track of the actual object
		ptr_ = &s;
		
		//Copy of current object for restoring the original values
		tempprocessingmethod_=s;

	  String str;
	  s.getCompletionTime().get(str);
		software_completion_time_->setText(str.c_str()); 
	}
	
	void DataProcessingVisualizer::store_()
	{
		try
		{
			String m(software_completion_time_->text().toStdString());
			DateTime date;
			
			try
			{
				date.set(m);
				ptr_->setCompletionTime(date);
			}
			catch(exception& e)
			{
				if(date.isNull())
				{
					std::string status= "Format of date in DATAPROCESSING is not correct.";
					emit sendStatus(status);
				}
			}
			
			tempprocessingmethod_=(*ptr_);
			
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to store the new processing method data. "<<e.what()<<endl;
		}
		
	}
	
	void DataProcessingVisualizer::reject_()
	{
		try
		{
			load(tempprocessingmethod_);
		}
		catch(exception e)
		{
			cout<<"Error while trying to restore original processing method data. "<<e.what()<<endl;
		}
	}

}
