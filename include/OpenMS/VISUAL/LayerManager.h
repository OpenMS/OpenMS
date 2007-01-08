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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_LAYERMANAGER_H
#define OPENMS_VISUAL_LAYERMANAGER_H

#include <OpenMS/VISUAL/UIC/LayerManagerTemplate.h>
#include <OpenMS/CONCEPT/Types.h>
#include <qlayout.h>

#include <vector>

namespace OpenMS 
{
	
	class LayerItem;
	
	/**
		@brief Displays a vertical list of LayerItem instances.
		
		The items consist of a checkbox, a label and can change the background color
		to show if they are activated (see example image).

		\image html LayerManager.png
		
		@ingroup Visual
	*/
	class LayerManager: public LayerManagerTemplate
	{
		Q_OBJECT
		
		public:
			/// Constructor
			LayerManager( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
			/// Destructor
			~LayerManager();
			/// Sets the checkbox of item @p i to @p b
			void setVisible(UnsignedInt index, bool b);
			/// Activates item @p index
			void activate(int index);
		
		public slots:
			/// Adds an item with label @p text
	  	virtual int addLayer( std::string text );
	  	/// Removes all items
	  	virtual void reset();
	  
		protected slots:
			virtual void itemVisibilityChanged(int index, bool b);
			virtual void itemActivated(int index);
			virtual void itemRemoveRequest(int index);
			virtual void itemPreferencesRequest(int index);
	
		protected:
			/// index of the active item
			SignedInt activated_item_;
			/// item count
			UnsignedInt count_;
	};

}
#endif // OPENMS_VISUAL_LAYERMANAGER_H

