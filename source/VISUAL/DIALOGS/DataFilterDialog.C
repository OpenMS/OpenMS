// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/DataFilterDialog.h>

//Qt includes
#include <QtGui/QDoubleValidator>
#include <QtGui/QIntValidator>
#include <QtGui/QMessageBox>

using namespace std;

namespace OpenMS
{

	DataFilterDialog::DataFilterDialog(DataFilters::DataFilter& filter, QWidget* parent)
		: QDialog(parent),
			filter_(filter)
	{
		setupUi(this);
		connect(ok_button_,SIGNAL(clicked()),this,SLOT(check_()));
		
		//set values for edit mode
		field_->setCurrentIndex((UInt)filter.field);
		op_->setCurrentIndex((UInt)filter.op);
		value_->setText(QString::number(filter.value));
	}
	
  void DataFilterDialog::check_()
  {
		QString field = field_->currentText();
		QString op = op_->currentText();
		QString value = value_->text();
		int tmp;
		
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
		if (field=="Charge")
		{
			QIntValidator v(this);
			if (v.validate(value,tmp)==QValidator::Invalid)
			{
				QMessageBox::warning(this,"Invalid value","An integer value is required!");
				return;
			}
		}
		
		//write result
		if (field=="Intensity") filter_.field = DataFilters::INTENSITY;
		else if (field=="Quality") filter_.field = DataFilters::QUALITY;
		else if (field=="Charge") filter_.field = DataFilters::CHARGE;

		if (op==">=") filter_.op = DataFilters::GREATER_EQUAL;
		else if (op=="=") filter_.op = DataFilters::EQUAL;
		else if (op=="<=") filter_.op = DataFilters::LESS_EQUAL;

		if (field=="Intensity" || field=="Quality") filter_.value = value.toDouble();
		else if (field=="Charge")  filter_.value = value.toInt();
		
		accept();	
  }

}



