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
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DGoToDialog.h>

#include <qlineedit.h>

using namespace std;

namespace OpenMS
{

	Spectrum1DGoToDialog::Spectrum1DGoToDialog( QWidget * parent, const char * name, WFlags fl)
		:	Spectrum1DGoToDialogTemplate(parent,name,fl),
	  min_pos_(0.0f),
	  max_pos_(1000.0f),
	  center_pos_((min_pos_+max_pos_)/2.0f)
	{
	  minPositionLineEdit->setText(QString().setNum(min_pos_));
	  maxPositionLineEdit->setText(QString().setNum(max_pos_));
	  centerPositionLineEdit->setText(QString().setNum(center_pos_));
	}
	
	Spectrum1DGoToDialog::~Spectrum1DGoToDialog()
	{
	
	}
	
	void Spectrum1DGoToDialog::setMinPosition(float min)
	{
	  // update model
	  min_pos_ = min; 
	  center_pos_ = (min_pos_+max_pos_)/2.0f;
	  // update view
	  minPositionLineEdit->setText(QString().setNum(min)); 
	  centerPositionLineEdit->setText(QString().setNum(center_pos_)); 
	}
	
	void Spectrum1DGoToDialog::setMaxPosition(float max)
	{
	  // update model
	  max_pos_ = max;
	  center_pos_ = (min_pos_+max_pos_)/2.0f;
	  // update view
	  maxPositionLineEdit->setText(QString().setNum(max));
	  centerPositionLineEdit->setText(QString().setNum(center_pos_)); 
	}
	
	float Spectrum1DGoToDialog::getMinPosition()
	{
	  return min_pos_;
	}
	
	float Spectrum1DGoToDialog::getMaxPosition()
	{
	  return max_pos_;
	}
	
	
	void Spectrum1DGoToDialog::gotoButton_clicked()
	{
	  // calculate translation of the center
	  float newCenter = centerPositionLineEdit->text().toFloat();
	  float translation = newCenter - center_pos_;
	  // update model
	  min_pos_ += translation;
	  max_pos_ += translation;
	  // conversion succeded
	  done(QDialog::Accepted);
	}
	
	void Spectrum1DGoToDialog::setVisibleAreaButton_clicked()
	{
	  // update model
	  min_pos_ = minPositionLineEdit->text().toFloat();
	  max_pos_ = maxPositionLineEdit->text().toFloat();
	  // conversion succeded
	  done(QDialog::Accepted);
	}

}//namespace OpenMS
