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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CVREFERENCE_H
#define OPENMS_DATASTRUCTURES_CVREFERENCE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief Controlled Vocabulary Reference
		
		Reference to a controlled vocabulary, defined in the first section of a mapping file.
		
		@ingroup Datastructures
	*/

			class OPENMS_DLLAPI CVReference
			{
				public:

				/// Default constructor
				CVReference();

				/// Copy constructor
				CVReference(const CVReference& rhs);

				/// Destructor
				virtual ~CVReference();

				/// Assignment operator
				CVReference& operator = (const CVReference& rhs);

				/** @name Accessors
				*/
				//@{
				/// sets the name of the CV reference
				void setName(const String& name);

				/// returns the name of the CV reference
				const String& getName() const;

				/// sets the CV identifier which is referenced
				void setIdentifier(const String& identifier);

				/// returns the CV identifier which is referenced
				const String& getIdentifier() const;
				//@}

				/** @name Predicates
				*/
				//@{
				/// equality operator 
				bool operator == (const CVReference& rhs) const;
				
				/// inequality operator
				bool operator != (const CVReference& rhs) const;
				//@}


				protected:

				String name_;

				String identifier_;
			};
			
			
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVREFERENCE_H
