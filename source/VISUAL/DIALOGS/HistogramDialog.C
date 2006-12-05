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


#include <OpenMS/VISUAL/DIALOGS/HistogramDialog.h>

#include <sstream>

#include <qpushbutton.h>
#include <qlayout.h>
#include <qtabwidget.h>
#include <qlabel.h>
#include <qnamespace.h>
#include <qstring.h>
#include <qlineedit.h>
#include <qradiobutton.h>
#include <qvbuttongroup.h>


using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	HistogramDialog::HistogramDialog(const Histogram<UnsignedInt,float>& distribution, QWidget* parent, const char* name) :QDialog(parent,name)
	{ 	
		setCaption("Intensity Distribution");
		
		//layout
		QGridLayout* layout = new QGridLayout(this,3,2);
		layout->setMargin(6);
		layout->setSpacing(6);
		layout->setRowStretch(0,100);
	
		//ok
		QPushButton* ok_button_ = new QPushButton("&OK",this);
		ok_button_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);
		connect(ok_button_,SIGNAL(clicked()),this,SLOT(accept()));
		layout->addWidget(ok_button_,1,1);
		
		//cancel
		QPushButton* cancel_button_ = new QPushButton("&Cancel",this);
		cancel_button_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);
		connect(cancel_button_,SIGNAL(clicked()),this,SLOT(reject()));
		layout->addWidget(cancel_button_,1,2);
		
		//distribution
		mw_ = new HistogramWidget(distribution, this, "HistogramWidget");
		mw_->showSplitters(true);
		layout->addMultiCellWidget(mw_,0,0,0,2);
		
		//resize dialog
		adjustSize();
	}
	
	HistogramDialog::~HistogramDialog()
	{ 	
		
	}
	
	double HistogramDialog::getLeftSplitter()
	{ 	
		return mw_->getLeftSplitter();	
	}
	
	double HistogramDialog::getRightSplitter()
	{ 	
		return mw_->getRightSplitter();	
	}
	
	void HistogramDialog::setLeftSplitter(double position)
	{
		mw_->setLeftSplitter(position);
	}
	
	void HistogramDialog::setRightSplitter(double position)
	{
		mw_->setRightSplitter(position);
	}

} //namespace OpenMS
