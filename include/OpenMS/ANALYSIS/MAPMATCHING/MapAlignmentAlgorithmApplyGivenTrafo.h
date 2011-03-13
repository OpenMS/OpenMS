// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMAPPLYGIVENTRAFO_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMAPPLYGIVENTRAFO_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

namespace OpenMS
{
	/**
		@brief This "map alignment" algorithm does not really perform an
		alignment.  Instead it simply applies a given transformation.
		
	  @htmlinclude OpenMS_MapAlignmentAlgorithmApplyGivenTrafo.parameters

		@ingroup MapAlignment
	*/
	class OPENMS_DLLAPI MapAlignmentAlgorithmApplyGivenTrafo
		: public MapAlignmentAlgorithm
	{
	 public:
		/// Default constructor
		MapAlignmentAlgorithmApplyGivenTrafo();

		/// Destructor
		virtual ~MapAlignmentAlgorithmApplyGivenTrafo();


		// Docu in base class
		virtual void alignPeakMaps( std::vector< MSExperiment<> >&, std::vector<TransformationDescription>& );

		// Docu in base class
		virtual void alignFeatureMaps( std::vector< FeatureMap<> >&, std::vector<TransformationDescription>& );

		// Docu in base class
		virtual void alignConsensusMaps( std::vector< ConsensusMap >&, std::vector<TransformationDescription>& );

		// Docu in base class
		virtual void alignPeptideIdentifications( std::vector< std::vector< PeptideIdentification > >&, std::vector<TransformationDescription>& );


    /**@brief Reads the "given transformations" from files that were set in the parameters.

		Note that this bypasses the factory concept, since this method is not
		defined in the base class.  It is provided here for "direct use" of this
		class.
		*/
		void readGivenTrafos(std::vector<TransformationDescription>&);

		/// Gets the default model (incl. parameters) for the alignment algorithm
		void getDefaultModel(String& model_type, Param& params);

		/// Creates a new instance of this class (for Factory)
		static MapAlignmentAlgorithm* create()
		{
			return new MapAlignmentAlgorithmApplyGivenTrafo();
		}
			
		/// Returns the product name (for the Factory)
		static String getProductName()
		{
			return "apply_given_trafo";
		}
		
	 private:

		/// Copy constructor intentionally not implemented -> private
		MapAlignmentAlgorithmApplyGivenTrafo(const MapAlignmentAlgorithmApplyGivenTrafo& );
		///Assignment operator intentionally not implemented -> private
		MapAlignmentAlgorithmApplyGivenTrafo& operator=(const MapAlignmentAlgorithmApplyGivenTrafo& );

	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMAPPLYGIVENTRAFO_H
