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
#include <OpenMS/VISUAL/DIALOGS/SmoothingDialog.h>
#include <QtGui/QLineEdit>
#include <QtGui/QRadioButton>
#include <QtGui/QCheckBox>
#include <QtGui/QPushButton>

#include <iostream> 
using namespace std;

namespace OpenMS
{

	SmoothingDialog::SmoothingDialog(QWidget* parent )
	    : QDialog(parent)
	{
		setupUi(this);
		
	  sgolay_order_line_edit->setText(QString().setNum(4));
	  kernel_widh_line_edit->setText(QString().setNum(1.0));
	  spacing_line_edit->setText(QString().setNum(0.2));
	  resampling_check_box->setChecked(false);
	  gaussian_radio_button->setChecked(false);
	  sgolay_radio_button->setChecked(true);
	}
	
	SmoothingDialog::~SmoothingDialog()
	{}
	
	void SmoothingDialog::setKernelWidth(float kw)
	{
	  // update view
	  kernel_widh_line_edit->setText(QString().setNum(kw));
	}
	
	void SmoothingDialog::setSGolayOrder(unsigned int o)
	{
	  // update view
	  sgolay_order_line_edit->setText(QString().setNum(o));
	}
	
	void SmoothingDialog::setSpacing(float r)
	{
	  // update view
	  spacing_line_edit->setText(QString().setNum(r));
	}
	
	void SmoothingDialog::setResampling(bool r)
	{
	  resampling_check_box->setChecked(r);
	}
	
	void SmoothingDialog::setGaussian(bool g)
	{
	  gaussian_radio_button->setChecked(g);
	}
	
	void SmoothingDialog::setSGolay(bool s)
	{
	  sgolay_radio_button->setChecked(s);
	}
	
	
	float SmoothingDialog::getKernelWidth()
	{
	  return QString(kernel_widh_line_edit->text()).toFloat();;
	}
	
	unsigned int SmoothingDialog::getSGolayOrder()
	{
	  return QString(sgolay_order_line_edit->text()).toUInt();
	}
	
	float SmoothingDialog::getSpacing()
	{
	  return QString(spacing_line_edit->text()).toFloat();
	}
	
	bool SmoothingDialog::getResampling()
	{
	  return resampling_check_box->isChecked();
	}
	
	bool SmoothingDialog::getGaussian()
	{
	  return gaussian_radio_button->isChecked();
	}
	
	bool SmoothingDialog::getSGolay()
	{
	  return sgolay_radio_button->isChecked();
	}
	
	void SmoothingDialog::startButton_clicked()
	{
	  // conversion succeded
	  done(QDialog::Accepted);
	}
	
	void SmoothingDialog::resetButton_clicked()
	{
	  sgolay_order_line_edit->setText(QString().setNum(4));
	  kernel_widh_line_edit->setText(QString().setNum(1.0));
	  spacing_line_edit->setText(QString().setNum(.2));
	  resampling_check_box->setChecked(false);
	  gaussian_radio_button->setChecked(false);
	  sgolay_radio_button->setChecked(true);
	}

} //namespace OpenMS
