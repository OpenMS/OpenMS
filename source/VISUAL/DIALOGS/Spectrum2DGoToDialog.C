// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QLineEdit>

using namespace std;

namespace OpenMS
{

	Spectrum2DGoToDialog::Spectrum2DGoToDialog(QWidget* parent)
		: QDialog(parent)
	{
		setupUi(this);
	}
	
	Spectrum2DGoToDialog::~Spectrum2DGoToDialog()
	{
	}
	
	void Spectrum2DGoToDialog::setRange(Real min_rt, Real max_rt, Real min_mz, Real max_mz)
	{
		min_rt_->setText(QString::number(min_rt));
		max_rt_->setText(QString::number(max_rt));
		min_mz_->setText(QString::number(min_mz));
		max_mz_->setText(QString::number(max_mz));
	}
	
	Real Spectrum2DGoToDialog::getMinRT() const
	{
	  return min_rt_->text().toFloat();
	}
	
	Real Spectrum2DGoToDialog::getMaxRT() const
	{
	  return max_rt_->text().toFloat();  
	}
	
	Real Spectrum2DGoToDialog::getMinMZ() const
	{
	  return min_mz_->text().toFloat();  
	}
	
	Real Spectrum2DGoToDialog::getMaxMZ() const
	{
	  return max_mz_->text().toFloat();  
	}

	void Spectrum2DGoToDialog::enableFeatureNumber(bool enabled)
	{
		feature_label_->setEnabled(enabled);
		nr_->setEnabled(enabled);
		feature_number_->setEnabled(enabled);
		//Reorder tab order
		if (enabled)
		{
			setTabOrder(feature_number_, ok_button_);
			setTabOrder(ok_button_, cancel_button_);
			setTabOrder(cancel_button_, min_mz_);
			setTabOrder(min_mz_, max_mz_);
			setTabOrder(max_mz_, min_rt_);
			setTabOrder(min_rt_, max_rt_);
		}
		else
		{
			setTabOrder(min_mz_, max_mz_);
			setTabOrder(max_mz_, min_rt_);
			setTabOrder(min_rt_, max_rt_);
			setTabOrder(max_rt_, ok_button_);
			setTabOrder(ok_button_, cancel_button_);
		}
	}
	
	String Spectrum2DGoToDialog::getFeatureNumber() const
	{
		return feature_number_->text(); 
	}

	bool Spectrum2DGoToDialog::showRange() const
	{
		if (feature_number_->text().trimmed()!="")
		{
			return false;
		}
		return true;
	}
	
} //namespace OpenMS
