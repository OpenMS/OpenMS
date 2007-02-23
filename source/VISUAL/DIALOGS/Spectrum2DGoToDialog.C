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
#include <QtGui/QLineEdit>

using namespace std;

namespace OpenMS
{

	Spectrum2DGoToDialog::Spectrum2DGoToDialog( QWidget * parent)
		: QDialog(parent)
	{
		setupUi(this);

	  min_mz_->setText(QString::number(1.0));
	  max_mz_->setText(QString::number(100.0));
	  min_rt_->setText(QString::number(1.0));
	  max_rt_->setText(QString::number(100.0));
	}
	
	Spectrum2DGoToDialog::~Spectrum2DGoToDialog()
	{
	
	}
	
	
	void Spectrum2DGoToDialog::setMinRT(double value)
	{
	  min_rt_->setText(QString::number(value));
	}
	
	void Spectrum2DGoToDialog::setMaxRT(double value)
	{
	  max_rt_->setText(QString::number(value));
	}
	
	float Spectrum2DGoToDialog::getMinRT()
	{
	  return min_rt_->text().toFloat();
	}
	
	float Spectrum2DGoToDialog::getMaxRT()
	{
	  return max_rt_->text().toFloat();  
	}
	
	void Spectrum2DGoToDialog::setMinMZ(double value)
	{
	  min_mz_->setText(QString::number(value));
	}
	
	void Spectrum2DGoToDialog::setMaxMZ(double value)
	{
	  max_mz_->setText(QString::number(value));
	}
	
	float Spectrum2DGoToDialog::getMinMZ()
	{
	  return min_mz_->text().toFloat();  
	}
	
	float Spectrum2DGoToDialog::getMaxMZ()
	{
	  return max_mz_->text().toFloat();  
	}

} //namespace OpenMS
