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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/BaselineFilteringDialog.h>
#include <QtGui/QLineEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QPushButton>

using namespace std;

namespace OpenMS
{

	BaselineFilteringDialog::BaselineFilteringDialog(QWidget* parent )
		: QDialog(parent)
	{
		setupUi(this);
		
	  struc_elem_line_edit->setText(QString().setNum(5));
	  spacing_line_edit->setText(QString().setNum(0.2));
	  resampling_check_box->setChecked(false);
	}
	
	BaselineFilteringDialog::~BaselineFilteringDialog()
	{}
	
	void BaselineFilteringDialog::setStrucElemWidth(float kw)
	{
	  // update view
	  struc_elem_line_edit->setText(QString().setNum(kw));
	}
	
	void BaselineFilteringDialog::setSpacing(float r)
	{
	  // update view
	  spacing_line_edit->setText(QString().setNum(r));
	}
	
	void BaselineFilteringDialog::setResampling(bool r)
	{
	  resampling_check_box->setChecked(r);
	}
	
	float BaselineFilteringDialog::getSpacing()
	{
	  return QString(spacing_line_edit->text()).toFloat();
	}
	
	bool BaselineFilteringDialog::getResampling()
	{
	  return resampling_check_box->isChecked();
	}
	
	float BaselineFilteringDialog::getStrucElemWidth()
	{
	  return QString(struc_elem_line_edit->text()).toFloat();;
	}
	
	void BaselineFilteringDialog::startButton_clicked()
	{
	  // conversion succeded
	  done(QDialog::Accepted);
	}
	
	void BaselineFilteringDialog::resetButton_clicked()
	{
	  struc_elem_line_edit->setText(QString().setNum(5));
	  spacing_line_edit->setText(QString().setNum(0.2));
	  resampling_check_box->setChecked(false);
	}

} //namespace OpenMS

