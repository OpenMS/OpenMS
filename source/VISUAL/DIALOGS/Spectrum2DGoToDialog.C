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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>
#include <qlineedit.h>

using namespace std;

namespace OpenMS
{

	Spectrum2DGoToDialog::Spectrum2DGoToDialog( QWidget * parent, const char * name, WFlags fl)
		: Spectrum2DGoToDialogTemplate(parent,name,fl),
	  area_(DRange<2>(0.0f, 0.0f, 1000.0f, 1000.0f)),
	  center_x_((area_.minX() + area_.maxX())/2.0f),
	  center_y_((area_.minY() + area_.maxY())/2.0f)
	{
	  centerPositionXLineEdit->setText(QString().setNum(center_x_));
	  centerPositionYLineEdit->setText(QString().setNum(center_y_));
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
	  center_x_ = (area_.minX() + area_.maxX())/2.0f;
	  // update view
	  minXLineEdit->setText(QString().setNum(area_.minX()));
	  centerPositionXLineEdit->setText(QString().setNum(center_x_));
	}
	
	void Spectrum2DGoToDialog::setMaxX(float maxX)
	{
	  // update model
	  area_.setMaxX(maxX);
	  center_x_ = (area_.minX() + area_.maxX())/2.0f;
	  // update view
	  maxXLineEdit->setText(QString().setNum(area_.maxX()));
	  centerPositionXLineEdit->setText(QString().setNum(center_x_));  
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
	  center_y_ = (area_.minY() + area_.maxY())/2.0f;
	  // update view
	  minYLineEdit->setText(QString().setNum(area_.minY()));
	  centerPositionYLineEdit->setText(QString().setNum(center_y_));
	  
	}
	
	void Spectrum2DGoToDialog::setMaxY(float maxY)
	{
	  // update model
	  area_.setMaxY(maxY);
	  center_y_ = (area_.minY() + area_.maxY())/2.0f;
	  // update view
	  maxYLineEdit->setText(QString().setNum(area_.maxY()));
	  centerPositionYLineEdit->setText(QString().setNum(center_y_));  
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
	  float new_center_x = centerPositionXLineEdit->text().toFloat();
	  float translationX = new_center_x - center_x_;
	  
	  float new_center_y = centerPositionYLineEdit->text().toFloat();
	  float translationY = new_center_y - center_y_;
	
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

} //namespace OpenMS
