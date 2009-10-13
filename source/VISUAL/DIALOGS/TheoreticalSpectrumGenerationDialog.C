// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>

namespace OpenMS
{
	TheoreticalSpectrumGenerationDialog::TheoreticalSpectrumGenerationDialog()
	{
		setupUi(this);
		
		connect(list_widget, SIGNAL(itemChanged(QListWidgetItem*)), this, SLOT(itemChanged(QListWidgetItem*)));
		
		// select b- and y-ions as residue types by default
		list_widget->item(0)->setCheckState(Qt::Unchecked);
		list_widget->item(1)->setCheckState(Qt::Checked);
		list_widget->item(2)->setCheckState(Qt::Unchecked);
		list_widget->item(3)->setCheckState(Qt::Unchecked);
		list_widget->item(4)->setCheckState(Qt::Checked);
		list_widget->item(5)->setCheckState(Qt::Unchecked);
		list_widget->item(6)->setCheckState(Qt::Unchecked);
		list_widget->item(7)->setCheckState(Qt::Unchecked);
		list_widget->item(8)->setCheckState(Qt::Unchecked);
	}
	
	void TheoreticalSpectrumGenerationDialog::itemChanged(QListWidgetItem* item)
	{
		if (item->text() == "Isotope clusters")
		{
			if (item->checkState() == Qt::Checked)
			{
				max_iso_label->setEnabled(true);
				max_iso_spinbox->setEnabled(true);
			}
			else
			{
				max_iso_label->setEnabled(false);
				max_iso_spinbox->setEnabled(false);
			}
		}
	}

} // namespace
