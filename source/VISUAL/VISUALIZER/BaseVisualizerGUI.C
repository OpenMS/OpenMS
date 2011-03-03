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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

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
	BaseVisualizerGUI::BaseVisualizerGUI(bool editable, QWidget* parent)
		: QWidget(parent),
			undo_button_(0),
			mainlayout_(0),
			row_(0),
			editable_(editable)
	{
		mainlayout_= new QGridLayout(this);
		mainlayout_->setMargin(0);
	}

	bool BaseVisualizerGUI::isEditable() const
	{
		return editable_;
	}

	void BaseVisualizerGUI::finishAdding_()
	{
		if(isEditable())
		{	
			addSeparator_();
			addButton_(undo_button_, "Undo");
			connect(undo_button_, SIGNAL(clicked()), this, SLOT(undo_()) );
		}
		addVSpacer_();
	}

	void BaseVisualizerGUI::addLabel_(QString label, UInt row)
	{
		QLabel *label_item = new QLabel(label, this);
		mainlayout_->addWidget(label_item, row, 0);
	}
	
	void BaseVisualizerGUI::addLabel_(QString label)
	{
		QLabel *label_item = new QLabel(label, this);
		mainlayout_->addWidget(label_item, row_, 0,1,3);
		row_++;
	}
	
	void BaseVisualizerGUI::addLineEdit_(QLineEdit*& ptr, QString label)
	{
		ptr = new QLineEdit(this);
		ptr->setMinimumWidth(180);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1, 1, 2);
		ptr->setReadOnly(!isEditable());
		row_++;
	}
	
	void BaseVisualizerGUI::addIntLineEdit_(QLineEdit*& ptr, QString label)
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
	
	void BaseVisualizerGUI::addDoubleLineEdit_(QLineEdit*& ptr, QString label)
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
	void BaseVisualizerGUI::addLineEditButton_(QString label, QLineEdit*& ptr1, QPushButton*& ptr2, QString buttonlabel)
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
	
	void BaseVisualizerGUI::addTextEdit_(QTextEdit*& ptr, QString label)
	{
		ptr = new QTextEdit(this);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1, 1, 2);
		ptr->setReadOnly(!isEditable());
		row_++;
	}
	
	void BaseVisualizerGUI::addComboBox_(QComboBox*& ptr,  QString label)
	{
		ptr = new QComboBox(this);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1,1, 2);
		ptr->blockSignals(true);
		row_++;	
	}
	
	void BaseVisualizerGUI::addBooleanComboBox_(QComboBox*& ptr,   QString label)
	{
		ptr = new QComboBox(this);
		ptr->insertItem(0,"false");
		ptr->insertItem(1,"true");
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1,1, 2);
		ptr->blockSignals(true);
		row_++;
	}
	
	void BaseVisualizerGUI::fillComboBox_(QComboBox*& ptr,  const std::string* items, int agr)
	{
		for(int i=0; i < agr; ++i)
		{
			ptr->insertItem(i,QString(items[i].c_str()));
		}	
	}
	
	void BaseVisualizerGUI::addButton_(QPushButton*& ptr, QString label )
	{
		ptr = new QPushButton(label, this);
		QHBoxLayout* box = new QHBoxLayout();
		box->addStretch(1);
		box->addWidget(ptr);
		mainlayout_->addLayout(box, row_, 0, 1, 3);
	
		ptr->setEnabled(isEditable());
		row_++;	
	}
	
	void BaseVisualizerGUI::addVSpacer_()
	{
		mainlayout_->setRowStretch(row_,1);
		row_++;
	}
	
	void BaseVisualizerGUI::add2Buttons_(QPushButton*& ptr1, QString label1, QPushButton*& ptr2, QString label2 )
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
	
	void BaseVisualizerGUI::addSeparator_()
	{
		QLabel* pLabel = new QLabel(this);
		pLabel->setFrameShape(QFrame::HLine); 
		mainlayout_->addWidget(pLabel, row_, 0, 1, 3);
		row_++;
	}

	void BaseVisualizerGUI::addListView_(QListWidget*& ptr, QString label)
	{
		ptr = new QListWidget(this);
		addLabel_(label, row_);
		mainlayout_->addWidget(ptr, row_, 1, 1, 2);
		row_++;
	}

}
