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
// $Maintainer: Cornelia Friedle$
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/Spectrum3DWidgetPDP.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <qlayout.h>
#include <qgroupbox.h>
#include <qcheckbox.h>
#include <qcombobox.h>
#include<qlabel.h>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{
		
	  Spectrum3DWidgetPDP::Spectrum3DWidgetPDP( Spectrum3DWidget* manager, QWidget* parent, const char* name, WFlags f)
			: PreferencesDialogPage(manager,parent,name,f)
	  {
			help_ = "This is the preferences dialog of 3D spectrum!"
							"<br>";
		
			QGridLayout* grid;	
			grid = new QGridLayout(this,1,1);		
			canvas_ = manager->client("Canvas", this);
			grid->addMultiCellWidget(canvas_, 0,1,0,1);
			load();
	  }
		
	  Spectrum3DWidgetPDP::~Spectrum3DWidgetPDP()
	  {  	
	  }
		
	  void Spectrum3DWidgetPDP::load()
	  {	
			canvas_ -> load();
	  }
		
	  void Spectrum3DWidgetPDP::save()
	  {
	  		canvas_->save();
	  }

	} // namespace Internal

} //namespace


