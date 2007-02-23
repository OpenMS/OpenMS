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

#include <QtGui/QLineEdit>

using namespace std;

namespace OpenMS
{

	Spectrum1DGoToDialog::Spectrum1DGoToDialog( QWidget * parent)
		:	QDialog(parent)
	{
		setupUi(this);
	  min_->setText(QString::number(1.0));
	  max_->setText(QString::number(100.0));
	}
	
	Spectrum1DGoToDialog::~Spectrum1DGoToDialog()
	{
	
	}
	
	void Spectrum1DGoToDialog::setMin(double value)
	{
	  min_->setText(QString::number(value));
	}
	
	void Spectrum1DGoToDialog::setMax(double value)
	{
	  max_->setText(QString::number(value));
	}
	
	float Spectrum1DGoToDialog::getMin()
	{
	  return min_->text().toFloat();
	}
	
	float Spectrum1DGoToDialog::getMax()
	{
	  return max_->text().toFloat();
	}

}//namespace OpenMS
