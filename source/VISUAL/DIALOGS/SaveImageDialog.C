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
//

#include <OpenMS/VISUAL/DIALOGS/SaveImageDialog.h>

// Qt
#include <qlayout.h>
#include <qpushbutton.h>  
#include <qcombobox.h>
#include <qlabel.h>
#include <qimage.h>
#include <qapplication.h>
#include <qvalidator.h>

// STL
#include <iostream>
#include <math.h>

namespace OpenMS
{


	SaveImageDialog::SaveImageDialog( QWidget * parent, const char * name, bool modal, WFlags f ):
	QDialog(parent,name,modal,f)
	{
		size_ratio_=1;
		//create dialog and layout (grid)
		QGridLayout* grid=new QGridLayout(this,5,2,10,5);
		
		//add accept/cancel buttons (and their layout)
		QBoxLayout* box_layout = new QHBoxLayout();
		grid->addLayout(box_layout,5,1);
		
		QPushButton* button = new QPushButton(this);
		button->setText("Cancel");
		connect(button,SIGNAL(clicked()),this,SLOT(reject()));
		box_layout->addWidget(button);
		
		button = new QPushButton(this);
		button->setText("Accept");
		button->setDefault(true);
		connect(button,SIGNAL(clicked()),this,SLOT(checkSize()));
		box_layout->addWidget(button);
	
		//add picture format selector
		QLabel* label = new QLabel("Picture format:",this);
		grid->addWidget(label,0,0);
		format_ = new QComboBox(this,false);
		format_->insertStrList(QImage::outputFormats());
		grid->addWidget(format_,0,1,Qt::AlignLeft);
		//set format to PNG/JPEG if available
		int png=-1;
		int jpeg=-1;
		for (int i=0;i<format_->count();i++)
		{	
			if (format_->text(i)=="PNG")
			{
				png=i;
			}
			if (format_->text(i)=="JPEG")
			{
				jpeg=i;
			}			
		}
		if (png!=-1)
		{
			format_->setCurrentItem(png);
		}
		else if (jpeg!=-1)
		{
			format_->setCurrentItem(jpeg);
		}
		
		//add size boxes and label (and their layout)
		label = new QLabel("Size (WxH):",this);
		grid->addWidget(label,1,0);
		
		QValidator* v = new QIntValidator(1,10000,this);
		box_layout = new QHBoxLayout();
		grid->addLayout(box_layout,1,1);
		size_x_ = new QLineEdit(this);
		size_x_->setValidator(v);
		connect(size_x_,SIGNAL(textChanged(const QString&)),this,SLOT(xSizeChanged(const QString&)));
		box_layout->addWidget(size_x_);
		label = new QLabel("x",this);
		box_layout->addWidget(label);
		size_y_ = new QLineEdit(this);
		size_y_->setValidator(v);
		connect(size_y_,SIGNAL(textChanged(const QString&)),this,SLOT(ySizeChanged(const QString&)));
		box_layout->addWidget(size_y_);
		label = new QLabel("pixel",this);
		box_layout->addWidget(label);		
		
		size_proportions_ = new QCheckBox("keep proportions",this);
		size_proportions_->setChecked(true);
		connect(size_proportions_,SIGNAL(toggled(bool)),this,SLOT(proportionsActivated(bool)));
		grid->addWidget(size_proportions_,2,1);
	}
	
	void SaveImageDialog::setSize(int x, int y)
	{
		QString* temp =new QString();
		temp->setNum(x);
		size_x_->setText(*temp);
		temp->setNum(y);
		size_y_->setText(*temp);
		setSizeRatio_(float(x)/float(y));
	}
		
	void SaveImageDialog::setSizeRatio_(float r)
	{
		if (r==0.0)
		{
			size_ratio_=1.0;
		}
		else
		{
			size_ratio_=r;
		}
	}	
			
	void SaveImageDialog::xSizeChanged(const QString& s)
	{
		if(size_proportions_->isChecked() && size_x_==qApp->focusWidget())
		{
			QString* temp = new QString();
			temp->setNum((int)rint(s.toInt()/size_ratio_));
			size_y_->setText(*temp);
		}
	}	
	
	void SaveImageDialog::ySizeChanged(const QString& s)
	{
		if(size_proportions_->isChecked() && size_y_==qApp->focusWidget())
		{
			QString* temp = new QString();
			temp->setNum((int)rint(s.toInt()*size_ratio_));
			size_x_->setText(*temp);
		}
	}	
	
	void SaveImageDialog::proportionsActivated(bool state)
	{
		if (state==true){
			setSizeRatio_(QString(size_x_->text()).toFloat()/QString(size_y_->text()).toFloat());
			}
	}	
	
	void SaveImageDialog::checkSize()
	{
		int x =size_x_->text().toInt();
		int y =size_y_->text().toInt();
		if (x>0 && y>0){
			accept();
			}
	}
	
	int SaveImageDialog::getXSize()
	{
		return size_x_->text().toInt();
	}
	
	int SaveImageDialog::getYSize()
	{
		return size_y_->text().toInt();
	}
	
	QString SaveImageDialog::getFormat()
	{
		return format_->currentText();
	}

} //namespace

