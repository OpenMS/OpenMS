// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DFEATUREPAIRSFILE_H
#define OPENMS_FORMAT_DFEATUREPAIRSFILE_H

#include <OpenMS/FORMAT/SchemaFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/FORMAT/HANDLERS/DFeaturePairsHandler.h>
#include <OpenMS/KERNEL/DFeatureMap.h>

namespace OpenMS
{
	
	/**
		@brief This class provides Input/Output functionality for the class DFeaturePairVector.
		
		The features pairs are computed by an instance of DBaseFeatureMatcher during the
		matching of MS maps. The features pairs are stored in a pseudo XML format. No schema
		has been developed yet therefore no validation can be performed.
	
		@ingroup FileIO
	*/
	class DFeaturePairsFile 
		: public Internal::SchemaFile
	{
		public:
			/** @name Constructors and Destructor */
			//@{
			///Default constructor
			DFeaturePairsFile();
			///Destructor
			~DFeaturePairsFile();
			//@}

			/** @name Accessors */
			//@{
			/// loads the file with name @p filename into @p pairs.
			template<Size D> 
			void load(String filename, DFeaturePairVector<D>& pairs) throw (Exception::FileNotFound, Exception::ParseError)
			{
				Internal::DFeaturePairsHandler<D> handler(pairs,filename);
				parse_(filename, &handler);
			}

			/// stores the pair vector @p pairs in file with name @p filename. 
			template<Size D> 
			void store(String filename, const DFeaturePairVector<D>& pairs) const throw (Exception::UnableToCreateFile)
			{
				if (pairs.empty()) return;
				Internal::DFeaturePairsHandler<D> handler(pairs,filename);
				save_(filename, &handler);
			}

			/** 
				@brief Convert pair vector into feature map

			*/
			static void pairsToFeatures(const DFeaturePairVector<2>& pairs, DFeatureMap<2>& map);

			//@}
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_DFEATUREPAIRSFILE_H
