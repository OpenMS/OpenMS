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

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/DataFilterDialog.h>

//Qt includes
#include <QtGui/QDoubleValidator>
#include <QtGui/QIntValidator>
#include <QtGui/QMessageBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

	DataFilterDialog::DataFilterDialog(DataFilters::DataFilter& filter, QWidget* parent)
		: QDialog(parent),
			filter_(filter)
	{
		setupUi(this);
		connect(ok_button_,SIGNAL(clicked()),this,SLOT(check_()));
		connect(field_,SIGNAL(activated(const QString&)),this,SLOT(field_changed_(const QString&)));
		connect(op_,SIGNAL(activated(const QString&)),this,SLOT(op_changed_(const QString&)));
		
		//set values for edit mode
		field_->setCurrentIndex((UInt)filter.field);
		op_->setCurrentIndex((UInt)filter.op);
		if(filter.field == DataFilters::META_DATA)
		{
			meta_name_field_->setText(filter.meta_name.toQString());
			// if the value stored in filter is numerical, get value from filter.value (a double)
			if(filter.value_is_numerical)
			{
				value_->setText(QString::number(filter.value));
			}
			else // get value from filter.value_string (a String)
			{
				value_->setText(filter.value_string.toQString());
			}
			meta_name_field_->setEnabled(true);
			meta_name_label_->setEnabled(true);
			if(filter.op == DataFilters::EXISTS)
			{
				value_->setEnabled(false);
				value_label_->setEnabled(false);
			}
		}
		else // for non meta data, the value is always numerical
		{
			value_->setText(QString::number(filter.value));
		}
		
		//focus the value if this is an edit operation
		if (filter!=DataFilters::DataFilter())
		{
			value_->selectAll();
			setTabOrder(value_, cancel_button_);
			setTabOrder(cancel_button_, ok_button_);
			setTabOrder(ok_button_, field_);
			setTabOrder(field_, meta_name_field_);
			setTabOrder(meta_name_field_, op_);
		}
	}
	
	void DataFilterDialog::field_changed_(const QString& field)
	{
		QString op(op_->currentText());
		if(field == "Meta data")
		{
			meta_name_field_->setEnabled(true);
			meta_name_label_->setEnabled(true);
		}
		else
		{
			meta_name_field_->setEnabled(false);
			meta_name_label_->setEnabled(false);
		}
	}
	
	void DataFilterDialog::op_changed_(const QString& op)
	{
		QString field(field_->currentText());
		if(op != "exists")
		{
			value_->setEnabled(true);
			value_label_->setEnabled(true);
		}
		else
		{
			value_->setEnabled(false);
			value_label_->setEnabled(false);
		}
	}
	
  void DataFilterDialog::check_()
  {
		QString field = field_->currentText();
		QString op = op_->currentText();
		QString value = value_->text();
		QString meta_name_field = meta_name_field_->text();
		bool not_numerical = true;
		int tmp;
		
		//meta data
		if (field == "Meta data")
		{
			QDoubleValidator dv(this);
			not_numerical = dv.validate(value,tmp) == QValidator::Invalid;
			
			if (meta_name_field.isEmpty())
			{
				QMessageBox::warning(this,"Insufficient arguments", "You must specify a meta name!");
				return;
			}
			if (op == "<=" || op == ">=")
			{
				if(not_numerical)
				{
					QMessageBox::warning(this,"Invalid value","<= and >= are defined for numerical values only!");
					return;
				}
			}
		}
		//intensity, quality, charge:
		else
		{
			if (op == "exists")
			{
				QMessageBox::warning(this,"Invalid operation","Operation \"exists\" is defined for meta data only!");
				return;
			}
			//double
			if (field=="Intensity" || field=="Quality")
			{
				QDoubleValidator v(this);
				if (v.validate(value,tmp)==QValidator::Invalid)
				{
					QMessageBox::warning(this,"Invalid value","A real value is required!");
					return;
				}
			}
			//int
			if (field=="Charge" || field=="Size")
			{
				QIntValidator v(this);
				if (v.validate(value,tmp)==QValidator::Invalid)
				{
					QMessageBox::warning(this,"Invalid value","An integer value is required!");
					return;
				}
			}
		}
		
		//write result
		if (field=="Intensity") filter_.field = DataFilters::INTENSITY;
		else if (field=="Quality") filter_.field = DataFilters::QUALITY;
		else if (field=="Charge") filter_.field = DataFilters::CHARGE;
		else if (field=="Size") filter_.field = DataFilters::SIZE;
		else if (field=="Meta data") 
		{
			filter_.field = DataFilters::META_DATA;
			filter_.meta_name = meta_name_field;
			if (not_numerical) // entered value is not numerical, store it in value_string (as String)
			{
				filter_.value_string = String(value);
				filter_.value_is_numerical = false;
			}
			else // value is numerical, store it in value (as double)
			{
				filter_.value = value.toDouble();
				filter_.value_is_numerical = true;
			}
		}
		
		if (op==">=") filter_.op = DataFilters::GREATER_EQUAL;
		else if (op=="=") filter_.op = DataFilters::EQUAL;
		else if (op=="<=") filter_.op = DataFilters::LESS_EQUAL;
		else if (op=="exists") filter_.op = DataFilters::EXISTS;

		if (field=="Intensity" || field=="Quality") filter_.value = value.toDouble();
		else if (field=="Charge" || field=="Size")  filter_.value = value.toInt();
		
		accept();	
  }

}
