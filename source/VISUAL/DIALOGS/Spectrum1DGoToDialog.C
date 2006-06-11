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
// $Id: Spectrum1DGoToDialog.C,v 1.2 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DGoToDialog.h>

using namespace std;
using namespace OpenMS;

Spectrum1DGoToDialog::Spectrum1DGoToDialog( QWidget * parent, const char * name, WFlags fl)
	:	Spectrum1DGoToDialogTemplate(parent,name,fl),
  minPos(0.0f),
  maxPos(1000.0f),
  centerPos((minPos+maxPos)/2.0f)
{
  minPositionLineEdit->setText(QString().setNum(minPos));
  maxPositionLineEdit->setText(QString().setNum(maxPos));
  centerPositionLineEdit->setText(QString().setNum(centerPos));
}

Spectrum1DGoToDialog::~Spectrum1DGoToDialog()
{

}

void Spectrum1DGoToDialog::setMinPosition(float min)
{
  // update model
  minPos = min; 
  centerPos = (minPos+maxPos)/2.0f;
  // update view
  minPositionLineEdit->setText(QString().setNum(min)); 
  centerPositionLineEdit->setText(QString().setNum(centerPos)); 
}

void Spectrum1DGoToDialog::setMaxPosition(float max)
{
  // update model
  maxPos = max;
  centerPos = (minPos+maxPos)/2.0f;
  // update view
  maxPositionLineEdit->setText(QString().setNum(max));
  centerPositionLineEdit->setText(QString().setNum(centerPos)); 
}

float Spectrum1DGoToDialog::getMinPosition()
{
  return minPos;
}

float Spectrum1DGoToDialog::getMaxPosition()
{
  return maxPos;
}


void Spectrum1DGoToDialog::gotoButton_clicked()
{
  // calculate translation of the center
  float newCenter = centerPositionLineEdit->text().toFloat();
  float translation = newCenter - centerPos;
  // update model
  minPos += translation;
  maxPos += translation;
  // conversion succeded
  done(QDialog::Accepted);
}

void Spectrum1DGoToDialog::setVisibleAreaButton_clicked()
{
  // update model
  minPos = minPositionLineEdit->text().toFloat();
  maxPos = maxPositionLineEdit->text().toFloat();
  // conversion succeded
  done(QDialog::Accepted);
}


