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
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/DataTable.h>

//QT
#include <QtGui/QLayout>
#include <QtGui/QWidget>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit> 
#include <QtGui/QLabel>
#include <QtGui/QTextEdit>
#include <QtGui/QGridLayout>
#include <QtGui/QPushButton>
#include <QtGui/QHBoxLayout>
#include <QtGui/QListWidget>

namespace OpenMS
{
	
	DataTable::DataTable(bool editable, QWidget *parent) 
		: QWidget(parent),
			editable_(editable)
	{
	  //for the actual metadata
		mainlayout_= new QGridLayout(this);
		mainlayout_->setMargin(0);
		row_=0;
	}
	
	
	void DataTable::addLabel_(QString label, UInt row)
	{
		QLabel *label_item = new QLabel(label, this);
		mainlayout_->addWidget(label_item, row, 0);
	}
	
	void DataTable::addLabel(QString label)
	{
		QLabel *label_item = new QLabel(label, this);
		mainlayout_->addWidget(label_item, row_, 0,1,3);
		row_++;
	}
	
	void DataTable::addLineEdit(QLineEdit*& ptr, QString label)
	{
		ptr = new QLineEdit(this);
		ptr->setMinimumWidth(180);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1, 1, 2);
		ptr->setReadOnly(!isEditable());
		row_++;
	}
	
	void DataTable::addIntLineEdit(QLineEdit*& ptr, QString label)
	{
		ptr = new QLineEdit(this);
		ptr->setMinimumWidth(180);
		QIntValidator *vali = new QIntValidator(ptr);
		ptr->setValidator(vali);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1, 1, 2);
		ptr->setReadOnly(!isEditable());
		row_++;
		
	}
	
	void DataTable::addDoubleLineEdit(QLineEdit*& ptr, QString label)
	{
		ptr = new QLineEdit(this);
		ptr->setMinimumWidth(180);
		QDoubleValidator *vali = new QDoubleValidator(ptr);
		ptr->setValidator(vali);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1, 1, 2);
		ptr->setReadOnly(!isEditable());
		row_++;
		
	}
	void DataTable::addLineEditButton(QString label, QLineEdit*& ptr1, QPushButton*& ptr2, QString buttonlabel)
	{
		QLabel* label_item = new QLabel(label, this);
		ptr1 = new QLineEdit(this);
		ptr1->setMinimumWidth ( 180 );
		ptr2 = new QPushButton(buttonlabel, this);
		mainlayout_->addWidget(label_item, row_, 0);
		mainlayout_->addWidget(ptr1, row_, 1);
		mainlayout_->addWidget(ptr2, row_, 2);
		
		ptr1->setReadOnly(!isEditable());
		ptr2->setEnabled(isEditable());
		row_++;
	}
	
	void DataTable::addTextEdit(QTextEdit*& ptr, QString label)
	{
		ptr = new QTextEdit(this);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1, 1, 2);
		ptr->setReadOnly(!isEditable());
		row_++;
	}
	
	
	void DataTable::addComboBox(QComboBox*& ptr,  QString label)
	{
		ptr = new QComboBox(this);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1,1, 2);
		ptr->blockSignals(true);
		row_++;	
	}
	
	void DataTable::addBooleanComboBox(QComboBox*& ptr,   QString label)
	{
		ptr = new QComboBox(this);
		ptr->insertItem(0,"false");
		ptr->insertItem(1,"true");
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1,1, 2);
		ptr->blockSignals(true);
		row_++;
	}
	
	void DataTable::fillComboBox(QComboBox*& ptr,  const std::string* items, int agr)
	{
		for(int i=0; i < agr; ++i)
		{
			ptr->insertItem(i,QString(items[i].c_str()));
		}	
	}
	
	void DataTable::addButton(QPushButton*& ptr, QString label )
	{
		ptr = new QPushButton(label, this);
		QHBoxLayout* box = new QHBoxLayout();
		box->addStretch(1);
		box->addWidget(ptr);
		mainlayout_->addLayout(box, row_, 0, 1, 3);
	
		ptr->setEnabled(isEditable());
		row_++;	
	}
	
	void DataTable::addVSpacer()
	{
		mainlayout_->setRowStretch(row_,1);
		row_++;
	}
	
	
	
	void DataTable::add2Buttons(QPushButton*& ptr1, QString label1, QPushButton*& ptr2, QString label2 )
	{
		ptr1 = new QPushButton(label1, this);
		ptr2 = new QPushButton(label2, this);
		QHBoxLayout* box = new QHBoxLayout();
		box->addStretch(1);
		box->addWidget(ptr1);
		box->addWidget(ptr2);
		mainlayout_->addLayout(box, row_, 0, 1, 3);
		row_++;
	}
	
	
	void DataTable::addSeparator()
	{
		QLabel* pLabel = new QLabel(this);
		pLabel->setFrameShape(QFrame::HLine); 
		mainlayout_->addWidget(pLabel, row_, 0, 1, 3);
		row_++;
	}
			

	bool DataTable::isEditable() const
	{
		return editable_;
	}


	void DataTable::addListView(QListWidget*& ptr, QString label)
	{
		ptr = new QListWidget(this);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1, 1, 2);
		row_++;
	}

} //namespace
