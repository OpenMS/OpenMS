// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/DemoDialog.h>

using namespace std;

namespace OpenMS
{

	DemoDialog::DemoDialog(QWidget* parent)
		: QDialog(parent),
			pages_(),
			current_page_(0),
			window_title_("DemoDialog")
	{
		setupUi(this);
		
		//signals and slots
		connect(previous_button_,SIGNAL(clicked()),this,SLOT(previous_()));
		connect(next_button_,SIGNAL(clicked()),this,SLOT(next_()));
	}
	
	DemoDialog::~DemoDialog()
	{ 	
	}

	void DemoDialog::setTitle(const String& title)
	{
		window_title_ = title;
		setWindowTitle(window_title_.toQString());
	}
	
	void DemoDialog::setPages(const StringList& pages)
	{
		pages_ = pages;
		show_(0);
	}

	void  DemoDialog::previous_()
	{
		show_(current_page_-1);
	}
	
	void  DemoDialog::next_()
	{
		show_(current_page_+1);
	}


	void  DemoDialog::show_(Size i)
	{
		//update current page
		current_page_ = i;
		
		//enable/disable buttons
		if (i==0)
		{
			previous_button_->setEnabled(false);
		}
		else
		{
			previous_button_->setEnabled(true);
		}
		if (i==pages_.size()-1)
		{
			next_button_->setEnabled(false);
		}
		else
		{
			next_button_->setEnabled(true);
		}

		//show page
		text_browser_->setSource(QUrl::fromLocalFile(pages_[current_page_].toQString()));
			
		//update window title
		String current_title = window_title_ + " (" + (i+1) + "/" + pages_.size() + ")";
		setWindowTitle(current_title.toQString());
	}
	
}


