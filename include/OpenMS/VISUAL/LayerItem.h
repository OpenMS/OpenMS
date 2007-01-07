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


#ifndef OPENMS_VISUAL_LAYERITEM_H
#define OPENMS_VISUAL_LAYERITEM_H

#include <OpenMS/VISUAL/UIC/LayerItemTemplate.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <qpopupmenu.h>

namespace OpenMS 
{

	/**
		@brief An item displayed in the LayerManager
		
		Visually it consists of a checkbox and a label.
		
		When the item is selected, the background color is blue.
	  
	  @ingroup Visual
	*/
	class LayerItem: public LayerItemTemplate
	{
		Q_OBJECT
		
		public:
			/// Constructor
			LayerItem( QWidget * parent = 0, const char * name = 0, WFlags fl = 0);
			/// Destructor
			~LayerItem();
			/// Sets the index associated with this item
			void setIndex(UnsignedInt index);
			/// Returns if this item of the LayerManager is activated
			bool isActivated();
			/// Returns the index associated with this item
			UnsignedInt getIndex() const;
			/// Return the label of this item
			String getLabel() const;
		
		public slots:
			/// Sets wether the checkbox is checked
	    virtual void changeState(bool state);
	    /// Changes the label
	    virtual void changeLabel(std::string l);
	    /// Activates the item
	  	void activate();
	  	/// Deactivates the item
	  	void deactivate();
	  		
	  protected slots:
	  	virtual void toggled(bool state);
	  
	  protected:
	  	UnsignedInt index_;
	  	bool activated_;
	  	String text_;
	  	virtual void mousePressEvent ( QMouseEvent * e );
			void contextMenuEvent( QContextMenuEvent * );
	};

}
#endif // OPENMS_VISUAL_LAYERITEM_H

