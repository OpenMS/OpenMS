// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/DATASTRUCTURES/Param.h>


// QT includes
#include <QtGui/QButtonGroup>

// STL includes
#include <iostream>

using namespace std;

namespace OpenMS
{

	TOPPViewOpenDialog::TOPPViewOpenDialog(const String& data_name, bool as_window, bool as_2d, bool cutoff, QWidget * parent)
		: QDialog(parent),
			map_as_2d_disabled_(false)
	{
		setupUi(this);
		
		//init map view
		QButtonGroup* button_group = new QButtonGroup(this);
		button_group->addButton(d1_);
		button_group->addButton(d2_);
		button_group->addButton(d3_);
		if (!as_2d)
		{
			d3_->setChecked(true);
		}
		else
		{
			d2_->setChecked(true);
		}

		//init intensity cutoff
		button_group = new QButtonGroup(this);
		button_group->addButton(cutoff_);
		button_group->addButton(nocutoff_);
		if (!cutoff)
		{
			nocutoff_->setChecked(true);
		}
		else
		{
			cutoff_->setChecked(true);
		}
		
		//init open as
		button_group = new QButtonGroup(this);
		button_group->addButton(window_);
		button_group->addButton(layer_);
		button_group->addButton(merge_);
		connect(button_group,SIGNAL(buttonClicked(QAbstractButton*)),this,SLOT(updateViewMode_(QAbstractButton*)));
		if (!as_window)
		{
			layer_->setChecked(true);
		}
		else
		{
			window_->setChecked(true);
		}
		connect(merge_combo_,SIGNAL(activated(int)),merge_,SLOT(click()));
		
		//set title
		setWindowTitle((String("Open data options for ") + data_name).toQString());
	}
	
	TOPPViewOpenDialog::~TOPPViewOpenDialog()
	{
	}

	bool TOPPViewOpenDialog::viewMapAs2D() const
	{
		if (d2_->isChecked()) return true;
		return false;
	}
	
	bool TOPPViewOpenDialog::viewMapAs1D() const
	{
		if (d1_->isChecked()) return true;
		return false;
	}
	
	bool TOPPViewOpenDialog::isCutoffEnabled() const
	{
		if (cutoff_->isChecked()) return true;
		return false;
	}
	
	bool TOPPViewOpenDialog::openAsNewWindow() const
	{
		if (window_->isChecked()) return true;
		return false;	
	}

	void TOPPViewOpenDialog::disableDimension(bool as_2d)
	{
		d1_->setChecked(!as_2d);
		d1_->setEnabled(false);
		d2_->setChecked(as_2d);
		d2_->setEnabled(false);
		d3_->setEnabled(false);
		map_as_2d_disabled_ = true;
	}
	
	void TOPPViewOpenDialog::disableCutoff(bool cutoff_on)
	{
		cutoff_->setChecked(cutoff_on);
		cutoff_->setEnabled(false);
		nocutoff_->setEnabled(false);
	}
	
	void TOPPViewOpenDialog::disableLocation(bool as_window)
	{
		window_->setEnabled(false);
		layer_->setEnabled(false);
		merge_->setEnabled(false);
		merge_combo_->setEnabled(false);
		if (as_window)
		{
			window_->setChecked(true);
		}
		else
		{
			layer_->setChecked(true);
		}
	}

	void TOPPViewOpenDialog::updateViewMode_(QAbstractButton* button)
	{
		if (button==layer_ || button==merge_)
		{
			d1_->setEnabled(false);
			d2_->setEnabled(false);
			d3_->setEnabled(false);
		}
		else if (!map_as_2d_disabled_)
		{
			d1_->setEnabled(true);
			d2_->setEnabled(true);
			d3_->setEnabled(true);
		}
	}

	void TOPPViewOpenDialog::setMergeLayers(const Map<Size,String>& layers)
	{
		//remove all items
		merge_combo_->clear();
		
		if (layers.size()!=0)
		{
			merge_->setEnabled(true);
			merge_combo_->setEnabled(true);
			UInt i=0;
			for (Map<Size,String>::const_iterator it=layers.begin(); it!=layers.end(); ++it)
			{
				merge_combo_->insertItem(i++,it->second.toQString(), (int)(it->first));
			}
		}
		else
		{
			merge_->setEnabled(false);
			merge_combo_->setEnabled(false);
		}
	}

	Int TOPPViewOpenDialog::getMergeLayer() const
	{
		if (merge_->isChecked())
		{
			return merge_combo_->itemData(merge_combo_->currentIndex()).toInt();
		}
		
		return -1;
	}

}



