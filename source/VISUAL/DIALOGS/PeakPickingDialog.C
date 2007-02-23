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
#include <OpenMS/VISUAL/DIALOGS/PeakPickingDialog.h>
#include <QtGui/QLineEdit>
#include <QtGui/QCheckBox>
#include <QtGui/QPushButton>

using namespace std;

namespace OpenMS
{

	PeakPickingDialog::PeakPickingDialog(QWidget* parent )
	    : QDialog(parent)
	{
		setupUi(this);
		
	  height_line_edit->setText(QString().setNum(200));
	  height_ms2_line_edit->setText(QString().setNum(50));
	  signal_to_noise_line_edit->setText(QString().setNum(5));
	  fwhm_line_edit->setText(QString().setNum(0.2));
	  opt_check_box->setChecked(false);
	}
	
	PeakPickingDialog::~PeakPickingDialog()
	{}
	
	
	void PeakPickingDialog::setPeakHeight(float height)
	{
	  // update view
	  height_line_edit->setText(QString().setNum(height));
	}
	
	void PeakPickingDialog::setPeakHeightMs2(float height)
	{
	  // update view
	  height_ms2_line_edit->setText(QString().setNum(height));
	}
	
	
	void PeakPickingDialog::setSignalToNoise(float sn)
	{
	  // update view
	  signal_to_noise_line_edit->setText(QString().setNum(sn));
	}
	
	void PeakPickingDialog::setFwhm(float fwhm)
	{
	  // update view
	  fwhm_line_edit->setText(QString().setNum(fwhm));
	}
	
	void PeakPickingDialog::setOptimization(bool opt)
	{
	  // update view
	  opt_check_box->setChecked(opt);
	}
	
	
	float PeakPickingDialog::getPeakHeight()
	{
	  return  QString(height_line_edit->text()).toFloat();
	}
	
	float PeakPickingDialog::getPeakHeightMs2()
	{
	  return QString(height_ms2_line_edit->text()).toFloat();
	}
	
	
	float PeakPickingDialog::getSignalToNoise()
	{
	  return QString(signal_to_noise_line_edit->text()).toFloat();
	}
	
	float PeakPickingDialog::getFwhm()
	{
	  return QString(fwhm_line_edit->text()).toFloat();
	}
	
	bool PeakPickingDialog::getOptimization()
	{
	  return opt_check_box->isChecked();
	}
	
	
	void PeakPickingDialog::startButton_clicked()
	{
	  // conversion succeded
	  done(QDialog::Accepted);
	}
	
	void PeakPickingDialog::resetButton_clicked()
	{
	  height_line_edit->setText(QString().setNum(200));
	  height_ms2_line_edit->setText(QString().setNum(50));
	  signal_to_noise_line_edit->setText(QString().setNum(5));
	  fwhm_line_edit->setText(QString().setNum(0.2));
	  opt_check_box->setChecked(false);
	}

}//namespace OpenMS

