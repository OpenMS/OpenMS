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
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm, Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHM_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHM_H

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

	/**
		@brief Base class for all map-alignment algorithms

		It takes two or more maps and corrects for retention time distortions.

		The input maps are transformed and the transformation description is returned.

		@improvement The maps should not be all loaded before the algorithm  - in order to save memory e.g. in the star-wise approach (Clemens)
	*/
	class OPENMS_DLLAPI MapAlignmentAlgorithm
		: public FactoryProduct,
			public ProgressLogger
	{
	 public:
		/// Default constructor
		MapAlignmentAlgorithm();

		/// Destructor
		virtual ~MapAlignmentAlgorithm();

		/**
		@brief Aligns peak maps

		@exception Exception::NotImplemented is thrown if an algorithm cannot align peak maps
		*/
		virtual void alignPeakMaps(std::vector< MSExperiment<> >&, std::vector<TransformationDescription>&);

		/**
		@brief Aligns feature maps

		@exception Exception::NotImplemented is thrown if an algorithm cannot align feature maps
		*/
		virtual void alignFeatureMaps(std::vector< FeatureMap<> >&, std::vector<TransformationDescription>&);

		/**
		@brief Aligns peptide identifications

		@exception Exception::NotImplemented is thrown if an algorithm cannot align peptide identifications
		*/
		virtual void alignPeptideIdentifications(std::vector< std::vector< PeptideIdentification > >&, std::vector<TransformationDescription>&);

		/// Register all derived classes in this method
		static void registerChildren();

	 private:
		/// Copy constructor is not implemented -> private
		MapAlignmentAlgorithm(const MapAlignmentAlgorithm& );
		/// Assignment operator is not implemented -> private
		MapAlignmentAlgorithm& operator=(const MapAlignmentAlgorithm& );

	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHM_H
