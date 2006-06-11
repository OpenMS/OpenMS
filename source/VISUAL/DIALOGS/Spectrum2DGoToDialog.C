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
// $Id: Spectrum2DGoToDialog.C,v 1.2 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Timo Sachsenberg$
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>
#include <qlineedit.h>

using namespace std;
using namespace OpenMS;

Spectrum2DGoToDialog::Spectrum2DGoToDialog( QWidget * parent, const char * name, WFlags fl)
	: Spectrum2DGoToDialogTemplate(parent,name,fl),
  area_(DRange<2>(0.0f, 0.0f, 1000.0f, 1000.0f)),
  centerX((area_.minX() + area_.maxX())/2.0f),
  centerY((area_.minY() + area_.maxY())/2.0f)
{
  centerPositionXLineEdit->setText(QString().setNum(centerX));
  centerPositionYLineEdit->setText(QString().setNum(centerY));
  minXLineEdit->setText(QString().setNum(area_.minX()));
  maxXLineEdit->setText(QString().setNum(area_.maxX()));
  minYLineEdit->setText(QString().setNum(area_.minY()));
  maxYLineEdit->setText(QString().setNum(area_.maxY()));
}

Spectrum2DGoToDialog::~Spectrum2DGoToDialog()
{

}


void Spectrum2DGoToDialog::setMinX(float minX)
{
  // update model
  area_.setMinX(minX);
  centerX = (area_.minX() + area_.maxX())/2.0f;
  // update view
  minXLineEdit->setText(QString().setNum(area_.minX()));
  centerPositionXLineEdit->setText(QString().setNum(centerX));
}

void Spectrum2DGoToDialog::setMaxX(float maxX)
{
  // update model
  area_.setMaxX(maxX);
  centerX = (area_.minX() + area_.maxX())/2.0f;
  // update view
  maxXLineEdit->setText(QString().setNum(area_.maxX()));
  centerPositionXLineEdit->setText(QString().setNum(centerX));  
}

float Spectrum2DGoToDialog::getMinX()
{
  return area_.minX();    
}

float Spectrum2DGoToDialog::getMaxX()
{
  return area_.maxX();    
}

void Spectrum2DGoToDialog::setMinY(float minY)
{
  // update model
  area_.setMinY(minY);
  centerY = (area_.minY() + area_.maxY())/2.0f;
  // update view
  minYLineEdit->setText(QString().setNum(area_.minY()));
  centerPositionYLineEdit->setText(QString().setNum(centerY));
  
}

void Spectrum2DGoToDialog::setMaxY(float maxY)
{
  // update model
  area_.setMaxY(maxY);
  centerY = (area_.minY() + area_.maxY())/2.0f;
  // update view
  maxYLineEdit->setText(QString().setNum(area_.maxY()));
  centerPositionYLineEdit->setText(QString().setNum(centerY));  
}

float Spectrum2DGoToDialog::getMinY()
{
  return area_.minY();  
}

float Spectrum2DGoToDialog::getMaxY()
{
  return area_.maxY();  
}

void Spectrum2DGoToDialog::gotoButton_clicked()
{
  // calculate translation of the center
  float newCenterX = centerPositionXLineEdit->text().toFloat();
  float translationX = newCenterX - centerX;
  
  float newCenterY = centerPositionYLineEdit->text().toFloat();
  float translationY = newCenterY - centerY;

  // update model
  area_.setMinX(area_.minX() + translationX);
  area_.setMaxX(area_.maxX() + translationX);

  area_.setMinY(area_.minY() + translationY);
  area_.setMaxY(area_.maxY() + translationY);
  
  // conversion succeded
  done(QDialog::Accepted);
}

void Spectrum2DGoToDialog::setVisibleAreaButton_clicked()
{
  // update model
  area_.setMinX(minXLineEdit->text().toFloat());
  area_.setMaxX(maxXLineEdit->text().toFloat());
  
  area_.setMinY(minYLineEdit->text().toFloat());
  area_.setMaxY(maxYLineEdit->text().toFloat()); 
  // conversion succeded
  done(QDialog::Accepted);
}

