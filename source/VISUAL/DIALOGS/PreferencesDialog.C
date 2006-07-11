// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/DIALOGS/PreferencesDialog.h>
#include <OpenMS/VISUAL/ListStack.h>

#include <iostream>

#include <qdialog.h>
#include <qpushbutton.h>
#include <qlayout.h>
#include <qmessagebox.h>

using namespace OpenMS;
using namespace std;


PreferencesDialog::PreferencesDialog()
:QDialog(),
pages_()
{
	QGridLayout* layout;
	QPushButton* button;
	
	setCaption("Preferences");

	//layout
	layout = new QGridLayout(this,2,5);
	layout->setSpacing(4);
	layout->setMargin(6);
	
	//buttons
	button = new QPushButton("&OK",this);
	button->setSizePolicy(QSizePolicy::QSizePolicy::Fixed,QSizePolicy::Fixed);
	connect(button,SIGNAL(clicked()),SLOT(ok_()));
	layout->addWidget(button,1,1);	

	button = new QPushButton("&Cancel",this);
	button->setSizePolicy(QSizePolicy::QSizePolicy::Fixed,QSizePolicy::Fixed);
	connect(button,SIGNAL(clicked()),SLOT(cancel_()));
	layout->addWidget(button,1,2);
	

	button = new QPushButton("&Apply",this);
	button->setSizePolicy(QSizePolicy::QSizePolicy::Fixed,QSizePolicy::Fixed);
	connect(button,SIGNAL(clicked()),SLOT(apply_()));
	layout->addWidget(button,1,3);
	
	button = new QPushButton("&Help",this);
	button->setSizePolicy(QSizePolicy::QSizePolicy::Fixed,QSizePolicy::Fixed);
	connect(button,SIGNAL(clicked()),SLOT(help_()));
	layout->addWidget(button,1,4);	

	//liststack
	stack_ = new ListStack(this);
	layout->addMultiCellWidget(stack_,0,0,0,4);
}

PreferencesDialog::~PreferencesDialog()
{
	
}

void PreferencesDialog::addPage(std::string name, PreferencesDialogPage* page, PreferencesManager* creator, PreferencesManager* parent) 
{
	pages_.push_back(page);
	stack_->addWidget(name,page,creator,parent);
	stack_->expand();
}


void PreferencesDialog::ok_()
{
	for (vector<PreferencesDialogPage*>::iterator it =pages_.begin();it!=pages_.end();++it)
	{
		(*it)->save();
	}
	accept();
}

void PreferencesDialog::cancel_()
{
	for (vector<PreferencesDialogPage*>::iterator it =pages_.begin();it!=pages_.end();++it)
	{
		(*it)->load();
	}
	reject();
}


void PreferencesDialog::apply_()
{
	for (vector<PreferencesDialogPage*>::iterator it =pages_.begin();it!=pages_.end();++it)
	{
		(*it)->save();
	}
}


void PreferencesDialog::help_()
{
	QMessageBox::information( this, "Dialog page help", QString(dynamic_cast<PreferencesDialogPage*>(stack_->activeWidget())->getHelpText().c_str()) );
}





