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
// $Id: LayerManager.h,v 1.5 2006/03/03 17:55:35 marc_sturm Exp $
// $Author: marc_sturm $
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
		@brief Used to display and remove LayerItems.
		
		  
		@ingroup Visual
	*/
	class LayerManager: public LayerManagerTemplate
	{
		Q_OBJECT
		
		public:
			LayerManager( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
			~LayerManager();
			void setVisible(UnsignedInt i, bool b);
			void activate(int index);
		
		public slots:
	  	virtual int addLayer( std::string label );
	  	virtual void reset();
	  
		protected slots:
			virtual void itemVisibilityChanged(int index, bool b);
			virtual void itemActivated(int index);
			virtual void itemRemoveRequest(int index);
	
		protected:
			//layout where the layer items are added
			QVBoxLayout* main_layout_;
			QVBoxLayout* layout_;
			std::vector<LayerItem*> items_;
			SignedInt activated_item_;
	};

}
#endif // OPENMS_VISUAL_LAYERMANAGER_H

