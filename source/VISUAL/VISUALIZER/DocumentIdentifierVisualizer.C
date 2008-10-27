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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/DocumentIdentifierVisualizer.h>

//QT
#include <QtGui/QLineEdit>

using namespace std;

namespace OpenMS
{

	DocumentIdentifierVisualizer::DocumentIdentifierVisualizer(bool editable, QWidget *parent) 
		: BaseVisualizer(editable, parent)
	{
		addLabel("Modify DocumentIdentifier information");		
		addSeparator();
		addLineEdit(identifier_, "Idenifier");
		finishAdding_();
	}
	
	
	void DocumentIdentifierVisualizer::load(DocumentIdentifier &h)
	{
	  ptr_ = &h;
		
		//Copy of current object for restoring the original values
		tempDocumentIdentifier_=h;
	  identifier_->setText(h.getIdentifier().c_str());
			
		identifier_->adjustSize();		
	}
	
	void DocumentIdentifierVisualizer::store_()
	{
		try
		{
			ptr_->setIdentifier(identifier_->text().toStdString());
			tempDocumentIdentifier_ = (*ptr_);
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to store the new document identifier. "<<e.what()<<endl;
		}
	}
	
	void DocumentIdentifierVisualizer::reject_()
	{
		try
		{
			load(tempDocumentIdentifier_);
		}
		catch(exception e)
		{
			cout<<"Error while trying to restore original document identifier.. "<<e.what()<<endl;
		} 
	}

}
