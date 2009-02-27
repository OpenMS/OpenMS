// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASVERTEX_H
#define OPENMS_VISUAL_TOPPASVERTEX_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QGraphicsItem>

namespace OpenMS
{
	class OPENMS_DLLAPI TOPPASVertex
		: public QGraphicsItem
	{		
		public:
			
			enum VertexType
			{
				VT_SOURCE,
				VT_TARGET,
				VT_TOOL
			};
			
			/// Constructor
			TOPPASVertex(const String& name, VertexType type);
			
			/// Returns the name of the tool
			const String& getName();
			
		protected:
			
			/// The name of the tool
			String name_;
			
			/// The type of this vertex
			VertexType vertex_type_;
			
	};
}

#endif
