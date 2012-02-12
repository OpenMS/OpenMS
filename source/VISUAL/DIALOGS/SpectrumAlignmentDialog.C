// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/VISUAL/DIALOGS/SpectrumAlignmentDialog.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>

// QT includes
#include <QtGui/QButtonGroup>

#include <vector>

namespace OpenMS
{
	SpectrumAlignmentDialog::SpectrumAlignmentDialog(Spectrum1DWidget* parent)
		: layer_indices_1_(),
			layer_indices_2_()
	{
		setupUi(this);
		
		QButtonGroup* button_group = new QButtonGroup(this);
		button_group->addButton(ppm);
		button_group->addButton(da);
		da->setChecked(true);
		
		Spectrum1DCanvas* cc = parent->canvas();
		for (UInt i = 0; i < cc->getLayerCount(); ++i)
		{
			const LayerData& layer = cc->getLayer(i);
			if (layer.flipped)
			{
				layer_list_2->addItem(layer.name.toQString());
				layer_indices_2_.push_back(i);
			}
			else
			{
				layer_list_1->addItem(layer.name.toQString());
				layer_indices_1_.push_back(i);
			}
		}
		// select first item of each list
		if (layer_list_1->count() > 0)
		{
			layer_list_1->setCurrentRow(0);
		}
		if (layer_list_2->count() > 0)
		{
			layer_list_2->setCurrentRow(0);
		}
	}
	
	Int SpectrumAlignmentDialog::get1stLayerIndex()
	{
		if (layer_list_1->count() == 0 || layer_list_1->currentRow() == -1)
		{
			return -1;
		}
		if (layer_indices_1_.size() > (Size)(layer_list_1->currentRow()))
		{
			return layer_indices_1_[(Size)(layer_list_1->currentRow())];
		}
		return -1;
	}
	
	Int SpectrumAlignmentDialog::get2ndLayerIndex()
	{
		if (layer_list_2->count() == 0 || layer_list_2->currentRow() == -1)
		{
			return -1;
		}
		if (layer_indices_2_.size() > (Size)(layer_list_2->currentRow()))
		{
			return layer_indices_2_[(Size)(layer_list_2->currentRow())];
		}
		return -1;
	}
	
} // namespace
