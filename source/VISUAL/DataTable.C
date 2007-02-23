// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free software; you can redistribute it and/or
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
// $Maintainer: stefan_heess $
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

using namespace OpenMS;

//Constructor
DataTable::DataTable(bool editable, QWidget *parent) 
	: QWidget(parent),
		editable_(editable)
{
  //for the actual metadata
	mainlayout_= new QGridLayout(this);
  mainlayout_->setSpacing(6);
  mainlayout_->setMargin(11);
	row_=0;
	column_=0;
}


//protected member functions.
void DataTable::addLabel_(const QString &labelName, UnsignedInt row)
{
	QLabel *label = new QLabel(labelName, this);
	mainlayout_->addWidget(label, row, 0);
}

//public member functios
void DataTable::addLabel(const QString &labelName)
{
	addLabel_(labelName, row_);
	row_++;
}

void DataTable::addLineEdit(QLineEdit* &ptr, const QString &label)
{
	ptr = new QLineEdit(this);
	addLabel_(label, row_);
	mainlayout_->addWidget(ptr, row_, 1);
	ptr->setReadOnly(!isEditable());
	row_++;
	
}

void DataTable::addLineEditButton(const QString &labelname, QLineEdit* &ptr1, QPushButton* &ptr2, const QString &buttonlabel)
{
	QWidget* buttonField = new QWidget(this);
	QHBoxLayout* bfLayout = new QHBoxLayout(buttonField);
	bfLayout->setSpacing(6);
  bfLayout->setMargin(11); 
	QLabel* label = new QLabel(labelname, buttonField);
	ptr1 = new QLineEdit(buttonField);
	ptr2 = new QPushButton(buttonlabel, buttonField);
	bfLayout->addWidget(label);
	bfLayout->addStretch(1);
	bfLayout->addWidget(ptr1);
	bfLayout->addStretch(1);
	bfLayout->addWidget(ptr2);
	mainlayout_->addWidget(buttonField, row_, row_, 0, 2);
	ptr1->setReadOnly(!isEditable());
	ptr2->setEnabled(isEditable());
	row_++;

}


void DataTable::addTextEdit(QTextEdit* &ptr ,const QString &label)
{
	ptr = new QTextEdit(this);
	addLabel_(label, row_);
	mainlayout_->addWidget(ptr, row_, 1);
	ptr->setReadOnly(!isEditable());
	row_++;
}


void DataTable::addComboBox(QComboBox* &ptr , const QString &label)
{
	ptr = new QComboBox(this);
	addLabel_(label, row_);
	mainlayout_->addWidget(ptr, row_, 1);
	ptr->blockSignals(true);
	row_++;	
}


void DataTable::fillComboBox(QComboBox* &ptr , const std::string* items, int agr)
{
	for(int i=0; i < agr; ++i)
	{
		ptr->insertItem(i,QString(items[i].c_str()));
	}	
}

void DataTable::addButton(QPushButton* &ptr, const QString &label )
{
	QWidget* buttonField = new QWidget(this);
	QHBoxLayout* bfLayout = new QHBoxLayout(buttonField);
	bfLayout->setSpacing(6);
  bfLayout->setMargin(11); 
	ptr = new QPushButton(label, buttonField);
	bfLayout->addStretch(1);
	bfLayout->addWidget(ptr);
	mainlayout_->addWidget(buttonField, row_, row_, 0, 2);
	ptr->setEnabled(isEditable());
	row_++;
	
}

void DataTable::addVSpacer()
{
	QSpacerItem* vspacer = new QSpacerItem( 0, 0, QSizePolicy::Expanding, QSizePolicy::Expanding );
	mainlayout_->addItem(vspacer, row_, row_, 0, 2);;
	row_++;
}



void DataTable::addHorizontalButtons(QPushButton* &ptr1, const QString &label1, QPushButton* &ptr2, const QString &label2 )
{
  QWidget* buttonField = new QWidget(this);
	QHBoxLayout* bfLayout = new QHBoxLayout(buttonField);
	bfLayout->setSpacing(6);
  bfLayout->setMargin(11); 
	ptr1 = new QPushButton(label1, buttonField);
	ptr2 = new QPushButton(label2, buttonField);
	bfLayout->addStretch(1);
	bfLayout->addWidget(ptr1);
	bfLayout->addWidget(ptr2);
	mainlayout_->addWidget(buttonField, row_, row_, 0, 2);
	row_++;
}

void DataTable::addSeperator()
{
	QLabel* pLabel = new QLabel(this);
	mainlayout_->addWidget(pLabel, row_, row_,0, 2);
	row_++;
	
}
		
		
void DataTable::addEmptyLine()
{
	QLabel *label1 = new QLabel(this);
	QLabel *label2 = new QLabel(this);
	mainlayout_->addWidget(label1, row_, 0);
	mainlayout_->addWidget(label2, row_, 1);
	row_++;
	
}

bool DataTable::isEditable() const
{
	return editable_;
}

