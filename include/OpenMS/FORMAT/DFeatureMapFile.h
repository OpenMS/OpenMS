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

#ifndef OPENMS_FORMAT_DFEATUREMAPFILE_H
#define OPENMS_FORMAT_DFEATUREMAPFILE_H

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/FORMAT/SchemaFile.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>

namespace OpenMS
{
	/**
  	@brief This class provides Input/Output functionality for DFeatureMaps
  		
  	@ingroup FileIO
  */
  class DFeatureMapFile : public Internal::SchemaFile
  {
	 public:
		/** @name Constructors and Destructor */
		//@{

		///Default constructor
		DFeatureMapFile();
		///Destructor
		~DFeatureMapFile();
		//@}

		/** @name Accessors */
		//@{
		/// loads the file with name @p filename into @p map. General case is not implemented!
		template<Size D> 
		void load(String filename, DFeatureMap<D>& map) throw (Exception::FileNotFound);
		
		/// loads the file with name @p filename into @p map.
		void load(String filename, DFeatureMap<2>& feature_map) throw (Exception::FileNotFound, Exception::ParseError);

		/// stores the map @p map in file with name @p filename. General case is not implemented!
		template<Size D> 
		void store(String filename, const DFeatureMap<D>& map) const throw (Exception::UnableToCreateFile);
					
		/// stores the map @p feature_map in file with name @p filename.
		void store(String filename, const DFeatureMap<2>& feature_map) const throw (Exception::UnableToCreateFile);
		//@}
		
		/// Mutable access to the options for loading/storing 
		PeakFileOptions& getOptions() { return options_; }

	protected:
		PeakFileOptions options_;
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_DEFEATUREMAPFILE_H
