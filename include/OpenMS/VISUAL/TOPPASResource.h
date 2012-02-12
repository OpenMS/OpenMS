// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASRESOURCE_H
#define OPENMS_VISUAL_TOPPASRESOURCE_H

#include <OpenMS/config.h>

#include <QtCore/QString>
#include <QtCore/QStringList>
#include <QtCore/QUrl>
#include <QtCore/QObject>

namespace OpenMS
{
	/**
		@brief Represents a data resource for TOPPAS workflows.
		
		Currently, the only supported type of resource is local files.
		
		@ingroup TOPPAS_elements
	*/
	class OPENMS_GUI_DLLAPI TOPPASResource
		: QObject
	{
		Q_OBJECT
		
		public:
			
			/// Constructor
			TOPPASResource(const QString& file);
			/// Constructor from URL
			TOPPASResource(const QUrl& url);
			/// Copy constructor
			TOPPASResource(const TOPPASResource& rhs);
			/// Destructor
			virtual ~TOPPASResource();
			/// Assignment operator
			TOPPASResource& operator= (const TOPPASResource& rhs);
			/// Writes this resource to the local file @p file
			void writeToFile(const QString& file_name);
			/// Returns the file name of the local file, or "" if it has not been written yet
			const QString& getLocalFile() const;
			/// Returns the URL of this resource
			const QUrl& getURL() const;
			/// Sets the URL of this resource from @p file
			void fromLocalFile(const QString& file);
			
			/// Supported schemes
			static QStringList supported_schemes;
			
		protected:
			
			/// The URL of this resource
			QUrl url_;
			/// The name of the local file
			QString file_name_;
	};
}

#endif
