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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASRESOURCES_H
#define OPENMS_VISUAL_TOPPASRESOURCES_H

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/TOPPASResource.h>

#include <QtCore/QString>
#include <QtCore/QObject>

namespace OpenMS
{
	/**
		@brief A dictionary mapping string keys to lists of TOPPASResource objects
		
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASResources
		: QObject
	{
		Q_OBJECT
		
		public:
			
			/// Constructor
			TOPPASResources();
			/// Copy constructor
			TOPPASResources(const TOPPASResources& rhs);
			/// Destructor
			virtual ~TOPPASResources();
			/// Assignment operator
			TOPPASResources& operator= (const TOPPASResources& rhs);
			/// Adds the (key,resource_list) pair to the dictionary
			void add(const QString& key, const QList<TOPPASResource>& resource_list);
			/// Returns the resource list that @p key is mapped to, or an empty list if @p key does not exist
			const QList<TOPPASResource>& get(const QString& key) const;
			/// Loads the dictionary from file @p file_name
			void load(const QString& file_name);
			/// Writes the dictionary to file @p file_name
			void store(const QString& file_name);
			/// Clears the dictionary
			void clear();
			
		protected:
			
			/// The dictionary
			Map<QString,QList<TOPPASResource> > map_;
			/// The empty list
			QList<TOPPASResource> empty_list_;
	};
}

#endif
