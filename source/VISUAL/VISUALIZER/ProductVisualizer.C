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
//  this library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/ProductVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	ProductVisualizer::ProductVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<Product>()
	{
		addLabel_("Modify processing method information.");	
		
		addSeparator_();  

		addDoubleLineEdit_(product_mz_, "m/z");
		addDoubleLineEdit_(product_window_low_, "Lower offset from target m/z");
		addDoubleLineEdit_(product_window_up_, "Upper offset from target m/z");
			
		finishAdding_();
	}
	
	void ProductVisualizer::update_()
	{
		product_mz_->setText(String( temp_.getMZ() ).c_str() );
		product_window_low_->setText(String( temp_.getIsolationWindowLowerOffset() ).c_str() );
		product_window_up_->setText(String( temp_.getIsolationWindowUpperOffset() ).c_str() );
	}
	
	void ProductVisualizer::store()
	{
		ptr_->setMZ(product_mz_->text().toFloat());
		ptr_->setIsolationWindowLowerOffset(product_window_low_->text().toFloat());		
		ptr_->setIsolationWindowUpperOffset(product_window_up_->text().toFloat());		
		
		temp_=(*ptr_);
	}
	
	void ProductVisualizer::undo_()
	{
		update_();
	}

}
