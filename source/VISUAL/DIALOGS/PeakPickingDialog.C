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
// $Maintainer: Timo Sachsenberg$
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/PeakPickingDialog.h>
#include <qlineedit.h>
#include <qradiobutton.h>
#include <qpushbutton.h>

using namespace std;
using namespace OpenMS;

PeakPickingDialog::PeakPickingDialog(QWidget* parent, const char* name , WFlags fl )
    : PeakPickingDialogTemplate(parent,name,fl),
    optimization_(false),
    height_(200),
    height_ms2_(50),
    signal_to_noise_(5),
    fwhm_(0.25)
{
  heightLineEdit->setText(QString().setNum(height_));
  heightMs2LineEdit->setText(QString().setNum(height_ms2_));
  signalToNoiseLineEdit->setText(QString().setNum(signal_to_noise_));
  fwhmLineEdit->setText(QString().setNum(fwhm_));
  optimizationRadioButton->setChecked(optimization_);
}

PeakPickingDialog::~PeakPickingDialog()
{}


void PeakPickingDialog::setPeakHeight(float height)
{
  // update model
  height_ = height;
  // update view
  heightLineEdit->setText(QString().setNum(height_));
}

void PeakPickingDialog::setPeakHeightMs2(float height)
{
  // update model
  height_ms2_ = height;
  // update view
  heightLineEdit->setText(QString().setNum(height_ms2_));
}


void PeakPickingDialog::setSignalToNoise(float sn)
{
  // update model
  signal_to_noise_ = sn;
  // update view
  signalToNoiseLineEdit->setText(QString().setNum(signal_to_noise_));
}

void PeakPickingDialog::setFwhm(float fwhm)
{
  // update model
  fwhm_ = fwhm;
  // update view
  fwhmLineEdit->setText(QString().setNum(fwhm_));
}

void PeakPickingDialog::setOptimization(bool opt)
{
  // update model
  optimization_ = opt;
  // update view
  optimizationRadioButton->setChecked(optimization_);
}


float PeakPickingDialog::getPeakHeight()
{
  return height_;
}

float PeakPickingDialog::getPeakHeightMs2()
{
  return height_ms2_;
}


float PeakPickingDialog::getSignalToNoise()
{
  return signal_to_noise_;
}

float PeakPickingDialog::getFwhm()
{
  return fwhm_;
}

bool PeakPickingDialog::getOptimization()
{
  return optimization_;
}


void PeakPickingDialog::startButton_clicked()
{
  // conversion succeded
  done(QDialog::Accepted);
}

void PeakPickingDialog::resetButton_clicked()
{
  height_ = 200;
  height_ms2_ = 50;
  signal_to_noise_ = 5;
  fwhm_ = 0.25;
  optimization_ = false;

  heightLineEdit->setText(QString().setNum(height_));
  heightMs2LineEdit->setText(QString().setNum(height_ms2_));
  signalToNoiseLineEdit->setText(QString().setNum(signal_to_noise_));
  fwhmLineEdit->setText(QString().setNum(fwhm_));
  optimizationRadioButton->setChecked(optimization_);
}

